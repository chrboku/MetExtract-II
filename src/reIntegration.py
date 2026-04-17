#    MetExtract II
#    Copyright (C) 2015
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

from __future__ import print_function, division, absolute_import
import gc
import logging
import time
import traceback
from copy import copy
from multiprocessing import Pool, Manager

from .Chromatogram import Chromatogram
from .chromPeakPicking.MassSpecWavelet import MassSpecWavelet
from .utils import Bunch, get_main_dir, getDBSuffix, getDBFormat
from .PolarsDB import PolarsDB
import polars as pl


# Constants for peak abundance calculation
peakAbundanceUseSignalsSides = 2


def smoothDataSeries(times, eic, windowLen, window, polynom):
    """
    Smooth a data series using various smoothing methods.

    Args:
        times: Time points
        eic: Extracted ion chromatogram data
        windowLen: Window length for smoothing
        window: Type of smoothing window
        polynom: Polynomial order for Savitzky-Golay filter

    Returns:
        Smoothed EIC data
    """
    if window == "None" or windowLen == 0:
        return eic

    # Import smoothing functions based on window type
    if window == "SavitzkyGolay":
        from .SGR import savitzkyGolay

        return savitzkyGolay(eic, window_size=windowLen, order=polynom)
    elif window == "Gaussian":
        from scipy.ndimage import gaussian_filter1d

        return gaussian_filter1d(eic, sigma=windowLen / 4.0)
    elif window == "MovingAverage":
        from numpy import convolve, ones

        kernel = ones(windowLen) / windowLen
        return convolve(eic, kernel, mode="same")
    else:
        return eic


class ReIntegrationProcessor:
    """
    Handles re-integration of chromatographic peaks for a single LC-HRMS file.
    """

    def __init__(self, forFile, chromPeakFile=None):
        """
        Initialize re-integration processor for one LC-HRMS file.

        Args:
            forFile: Path to the LC-HRMS data file
            chromPeakFile: Kept for API compatibility (no longer used)
        """
        self.chromPeakFile = chromPeakFile  # kept for backward compatibility
        self.forFile = forFile
        self.queue = None
        self.pID = None
        self.results = []

    def findPeakFor(
        self,
        mz,
        rt,
        ppm,
        scanEvent,
        maxRTShift,
        smoothingWindow,
        smoothingWindowSize,
        smoothingWindowPolynom,
    ):
        """
        Find chromatographic peak for one detected feature.

        Args:
            mz: m/z value
            rt: retention time
            ppm: mass accuracy in ppm
            scanEvent: scan event filter
            maxRTShift: maximum RT shift allowed
            smoothingWindow: smoothing method
            smoothingWindowSize: window size for smoothing
            smoothingWindowPolynom: polynomial order for smoothing

        Returns:
            Tuple of (area, abundance, SNR) or None if no peak found
        """
        eic, times, scanids, mzs = self.chromatogram.getEIC(mz, ppm, filterLine=scanEvent)
        eicSmoothed = smoothDataSeries(
            times,
            eic,
            windowLen=smoothingWindowSize,
            window=smoothingWindow,
            polynom=smoothingWindowPolynom,
        )
        ret = self.peakPicker.getPeaksFor(times, eicSmoothed)

        best = None

        for ri, r in enumerate(ret):
            if abs(times[r.peakIndex] / 60.0 - rt) <= maxRTShift and (best is None or abs(times[r.peakIndex] / 60.0 - rt) < abs(times[ret[best].peakIndex] / 60.0 - rt)):
                best = ri

        if best is not None:
            lb = int(
                min(
                    ret[best].peakIndex - peakAbundanceUseSignalsSides,
                    ret[best].peakIndex - peakAbundanceUseSignalsSides,
                )
            )
            rb = (
                int(
                    max(
                        ret[best].peakIndex + peakAbundanceUseSignalsSides,
                        ret[best].peakIndex + peakAbundanceUseSignalsSides,
                    )
                )
                + 1
            )
            peak = eic[lb:rb]

            peakAbundance = max(peak)
            return ret[best].peakArea, peakAbundance, ret[best].peakSNR

        return None

    def processFeaturePair(
        self,
        rowNum,
        mz,
        rt,
        lmz,
        ionMode,
        fileName,
        positiveScanEvent,
        negativeScanEvent,
        ppm,
        maxRTShift,
        smoothingWindow,
        smoothingWindowSize,
        smoothingWindowPolynom,
        addPeakArea,
        addPeakAbundance,
        addPeakSNR,
    ):
        """
        Re-integrate one detected feature pair.

        Returns:
            Dictionary with re-integration results or None
        """
        scanEvent = ""
        if ionMode == "+":
            scanEvent = positiveScanEvent
        elif ionMode == "-":
            scanEvent = negativeScanEvent
        else:
            logging.error(f"Undefined scan event: {ionMode}")
            return None

        result = {
            "Num": rowNum,
            "fileName": fileName,
        }

        # Re-integrate N (natural/unlabeled) peak
        r = None
        try:
            r = self.findPeakFor(
                mz,
                rt,
                ppm,
                scanEvent,
                maxRTShift,
                smoothingWindow,
                smoothingWindowSize,
                smoothingWindowPolynom,
            )
        except Exception as exc:
            logging.error(f"   - Reintegration failed for feature pair (N) {self.forFile} ({mz} {rt}) [{exc}]")

        nFound = False
        areaN = 0.0
        if r is not None:
            area, abundance, snr = r
            areaN = area
            result[f"{fileName}_Found"] = "Reintegrated"
            result[f"{fileName}_Area_N"] = area
            result[f"{fileName}_Abundance_N"] = abundance
            result[f"{fileName}_SNR_N"] = snr
            nFound = True

        # Re-integrate L (labeled) peak
        r = None
        try:
            r = self.findPeakFor(
                lmz,
                rt,
                ppm,
                scanEvent,
                maxRTShift,
                smoothingWindow,
                smoothingWindowSize,
                smoothingWindowPolynom,
            )
        except Exception as exc:
            logging.error(f"   - Reintegration failed for feature pair (L) {self.forFile} ({lmz} {rt}) [{exc}]")

        lFound = False
        areaL = 0.0
        if r is not None:
            area, abundance, snr = r
            areaL = area
            result[f"{fileName}_Found"] = "Reintegrated"
            result[f"{fileName}_Area_L"] = area
            result[f"{fileName}_Abundance_L"] = abundance
            result[f"{fileName}_SNR_L"] = snr
            lFound = True

        # Track for database if found
        if nFound or lFound:
            result["found"] = True
        else:
            result["found"] = False

        return result

    def processFile(self, params):
        """
        Re-integrate all detected feature pairs in this LC-HRMS data file.

        Args:
            params: Bunch object containing all processing parameters

        Returns:
            List of re-integration results
        """
        startProc = time.time()
        logging.info(f"   Reintegration started for file {self.forFile}")

        if self.queue is not None:
            self.queue.put(Bunch(pid=self.pID, mes="start"))
            self.queue.put(Bunch(pid=self.pID, mes="text", val=str(self.forFile)))

        try:
            # Import LC-HRMS data file
            logging.info(f"   Reading file {self.forFile}")
            self.chromatogram = Chromatogram()
            self.chromatogram.parse_file(self.forFile, intensityCutoff=params.reintegrateIntensityCutoff)

            # Initialize peak picking algorithm (using native Python wavelet implementation)
            self.peakPicker = MassSpecWavelet(scales=params.scales, snrTh=params.snrTH, minScans=1)

            scanEventsPerPolarity = self.chromatogram.getFilterLinesPerPolarity()

            if self.queue is not None:
                self.queue.put(Bunch(pid=self.pID, mes="max", val=len(params.features)))

            # Process each feature pair
            for idx, feature in enumerate(params.features):
                if self.queue is not None:
                    self.queue.put(Bunch(pid=self.pID, mes="value", val=idx))

                try:
                    # Check if scan event is available for this ionization mode
                    shouldProcess = False
                    if feature["ionMode"] == "+" and params.positiveScanEvent in scanEventsPerPolarity.get("+", []):
                        shouldProcess = True
                    elif feature["ionMode"] == "-" and params.negativeScanEvent in scanEventsPerPolarity.get("-", []):
                        shouldProcess = True

                    if shouldProcess:
                        result = self.processFeaturePair(
                            rowNum=feature["Num"],
                            mz=feature["MZ"],
                            rt=feature["RT"],
                            lmz=feature["L_MZ"],
                            ionMode=feature["ionMode"],
                            fileName=params.fileName,
                            positiveScanEvent=params.positiveScanEvent,
                            negativeScanEvent=params.negativeScanEvent,
                            ppm=params.ppm,
                            maxRTShift=params.maxRTShift,
                            smoothingWindow=params.smoothingWindow,
                            smoothingWindowSize=params.smoothingWindowSize,
                            smoothingWindowPolynom=params.smoothingWindowPolynom,
                            addPeakArea=params.addPeakArea,
                            addPeakAbundance=params.addPeakAbundance,
                            addPeakSNR=params.addPeakSNR,
                        )
                        if result and result.get("found", False):
                            del result["found"]
                            self.results.append(result)

                except Exception as ex:
                    logging.error(f"Error processing feature {feature.get('Num', 'unknown')}: {ex}")
                    traceback.print_exc()

            # Clean up
            self.chromatogram.freeMe()
            gc.collect()

            logging.info(f"   Reintegration finished for file {self.forFile} ({(time.time() - startProc) / 60.0:.1f} minutes)")

        except Exception as e:
            logging.error(f"Error during reintegration of {self.forFile}: {e}")
            traceback.print_exc()

        finally:
            if self.queue is not None:
                self.queue.put(Bunch(pid=self.pID, mes="end"))

        return self.forFile, self.results

    def setMultiProcessingQueue(self, qu, pID):
        """Set multiprocessing queue for progress reporting."""
        self.queue = qu
        self.pID = pID


def processReIntegrationFile(processor):
    """
    Helper function for multiprocessing.

    Args:
        processor: ReIntegrationProcessor instance with params set

    Returns:
        Tuple of (file_path, results_list)
    """
    return processor.processFile(processor.params)


def reIntegrateResultsFile(
    file,
    sheet_name,
    new_sheet_name,
    fDict,
    addPeakArea=True,
    addPeakAbundance=True,
    addPeakSNR=True,
    ppm=5.0,
    maxRTShift=0.25,
    scales=[3, 19],
    reintegrateIntensityCutoff=0,
    snrTH=1,
    smoothingWindow="None",
    smoothingWindowSize=0,
    smoothingWindowPolynom=0,
    positiveScanEvent="None",
    negativeScanEvent="None",
    pw=None,
    selfObj=None,
    cpus=1,
    start=0,
):
    """
    Re-integrate all LC-HRMS data files with the grouped feature pairs results using PolarsDB.

    This function:
    1. Loads the last sheet from the Excel results file using PolarsDB
    2. Creates re-integration tasks for each file
    3. Processes files in parallel using multiprocessing
    4. Collects results and updates the dataframe
    5. Saves the updated results back to Excel

    Args:
        file: Path to the results file (without extension)
        fDict: Dictionary mapping file paths to file names
        addPeakArea: Whether to add peak area columns
        addPeakAbundance: Whether to add peak abundance columns
        addPeakSNR: Whether to add peak SNR columns
        ppm: Mass accuracy in ppm
        maxRTShift: Maximum RT shift allowed (in minutes)
        scales: Wavelet scales for peak picking
        reintegrateIntensityCutoff: Intensity cutoff for reintegration
        snrTH: SNR threshold
        smoothingWindow: Smoothing method
        smoothingWindowSize: Window size for smoothing
        smoothingWindowPolynom: Polynomial order for smoothing
        positiveScanEvent: Positive scan event name
        negativeScanEvent: Negative scan event name
        pw: Progress wrapper for UI updates
        selfObj: Reference to main object (for job termination)
        cpus: Number of CPU cores to use
        start: Start time for elapsed time calculation
    """
    logging.info("Starting re-integration of missed peaks using PolarsDB")
    excel_file = file.replace(".xlsx", ".tsv").replace(".tsv", ".txt").replace(".txt", "") + ".xlsx"

    # Load the PolarsDB once
    try:
        plDB = PolarsDB(excel_file, format="xlsx", load_all_tables=True)
    except Exception as e:
        logging.error(f"Failed to load results file: {e}")
        return

    # Load the dataframe
    results_df = plDB.get_table(sheet_name).clone()

    # Check required columns
    required_cols = ["Num", "MZ", "RT", "L_MZ", "Ionisation_Mode", "Xn", "Charge"]
    missing_cols = [col for col in required_cols if col not in results_df.columns]
    if missing_cols:
        logging.error(f"Missing required columns: {missing_cols}")
        plDB.close()
        return

    # Create re-integration processors for each file
    processors = []
    for fileIdx, (filePath, fileName) in enumerate(fDict.items()):
        processor = ReIntegrationProcessor(filePath)

        # Prepare feature data for processing
        # select those features that need reintegration
        features_list = []

        for row in results_df.iter_rows(named=True):
            if row[f"{fileName}_Found"] is not None:
                pass
            else:
                features_list.append(
                    {
                        "Num": row["Num"],
                        "MZ": row["MZ"],
                        "RT": row["RT"],
                        "L_MZ": row["L_MZ"],
                        "ionMode": row["Ionisation_Mode"],
                    }
                )

        # Create params object
        params = Bunch(
            features=features_list,
            fileName=fileName,
            ppm=ppm,
            maxRTShift=maxRTShift,
            scales=scales,
            reintegrateIntensityCutoff=reintegrateIntensityCutoff,
            snrTH=snrTH,
            smoothingWindow=smoothingWindow,
            smoothingWindowSize=smoothingWindowSize,
            smoothingWindowPolynom=smoothingWindowPolynom,
            positiveScanEvent=positiveScanEvent,
            negativeScanEvent=negativeScanEvent,
            addPeakArea=addPeakArea,
            addPeakAbundance=addPeakAbundance,
            addPeakSNR=addPeakSNR,
        )
        processor.params = params

        processors.append(processor)

    # Set up multiprocessing
    p = Pool(processes=min(len(fDict), cpus), maxtasksperchild=1)
    manager = Manager()
    queue = manager.Queue()

    # Set queue for each processor
    for idx, processor in enumerate(processors):
        processor.setMultiProcessingQueue(queue, idx + 1)

    if pw is not None and selfObj is not None:
        from .utils import CallBackMethod

        pw.setCloseCallback(closeCallBack=CallBackMethod(_target=interruptReIntegrationProcessing, pool=p, selfObj=selfObj).getRunMethod())

    # Start the multiprocessing pool
    imap_results = p.imap_unordered(processReIntegrationFile, processors)

    # Monitor progress
    loop = True
    freeSlots = list(range(min(len(fDict), cpus)))
    assignedThreads = {}
    completed = 0

    while loop and (selfObj is None or not selfObj.terminateJobs):
        completed_now = imap_results._index
        if completed_now == len(fDict):
            loop = False
        else:
            # Process queue messages
            mess = {}
            while not queue.empty():
                try:
                    mes = queue.get(block=False, timeout=1)
                    if mes.pid not in mess:
                        mess[mes.pid] = {}
                    mess[mes.pid][mes.mes] = mes
                except:
                    break

            # Handle thread assignment
            for v in mess.values():
                if "end" in v.keys():
                    mes = v["end"]
                    if mes.pid in assignedThreads:
                        freeS = assignedThreads[mes.pid]
                        if freeS != -1 and pw is not None:
                            pw.getCallingFunction(freeS + 1)("text")("")
                            pw.getCallingFunction(freeS + 1)("value")(0)
                        assignedThreads[mes.pid] = -1
                        if freeS != -1:
                            freeSlots.append(freeS)

            for v in mess.values():
                if "start" in v.keys():
                    mes = v["start"]
                    if len(freeSlots) > 0:
                        w = freeSlots.pop()
                        assignedThreads[mes.pid] = w

            # Update progress bars
            if pw is not None:
                for v in mess.values():
                    for mes in v.values():
                        if mes.mes in ["log", "max", "value", "text"]:
                            if mes.pid in assignedThreads and assignedThreads[mes.pid] != -1:
                                pw.getCallingFunction(assignedThreads[mes.pid] + 1)(mes.mes)(mes.val)

            # Update overall progress
            elapsed = (time.time() - start) / 60.0
            hours = ""
            if elapsed >= 60.0:
                hours = f"{int(elapsed // 60)} hours "
            mins = f"{elapsed % 60.0:.2f} mins"

            if pw is not None:
                pw.getCallingFunction()("value")(completed_now)
                pw.getCallingFunction()("text")(f"<p align='right' >{hours}{mins} elapsed</p>\n\n{completed_now} / {len(fDict)} files done ({min(cpus, len(fDict))} parallel)")

            # Allow UI updates
            try:
                from PySide6 import QtWidgets

                QtWidgets.QApplication.processEvents()
            except:
                pass

            time.sleep(0.5)

    elapsed = (time.time() - start) / 60.0
    hours = ""
    if elapsed >= 60.0:
        hours = f"{int(elapsed // 60)} hours "
    mins = f"{elapsed % 60.0:.2f} mins"
    logging.info(f"Re-integrating finished after {hours}{mins}. Collecting results...")

    # Collect results from all processors
    all_results = {}
    for filePath, results in imap_results:
        for result in results:
            rowNum = result["Num"]
            if rowNum not in all_results:
                all_results[rowNum] = {}
            all_results[rowNum].update(result)

    # Update the dataframe with new columns and values
    new_columns = {}

    # Initialize new columns
    for filePath, fileName in fDict.items():
        if addPeakArea:
            if f"{fileName}_Area_N" not in results_df.columns:
                new_columns[f"{fileName}_Area_N"] = [None] * len(results_df)
            if f"{fileName}_Area_L" not in results_df.columns:
                new_columns[f"{fileName}_Area_L"] = [None] * len(results_df)

        if addPeakAbundance:
            if f"{fileName}_Abundance_N" not in results_df.columns:
                new_columns[f"{fileName}_Abundance_N"] = [None] * len(results_df)
            if f"{fileName}_Abundance_L" not in results_df.columns:
                new_columns[f"{fileName}_Abundance_L"] = [None] * len(results_df)

        if addPeakSNR:
            if f"{fileName}_SNR_N" not in results_df.columns:
                new_columns[f"{fileName}_SNR_N"] = [None] * len(results_df)
            if f"{fileName}_SNR_L" not in results_df.columns:
                new_columns[f"{fileName}_SNR_L"] = [None] * len(results_df)

    # Add new columns to dataframe
    for col_name, col_data in new_columns.items():
        results_df = results_df.with_columns(pl.Series(col_name, col_data))

    # Update values from results
    print(f"Updating results in all_results table")

    for row_idx, row in enumerate(results_df.iter_rows()):
        num = int(results_df[row_idx, "Num"])
        if num not in all_results.keys():
            pass
        else:
            for key, value in all_results[num].items():
                if key not in ["Num", "fileName", "areaN", "areaL", "found"] and value is not None:
                    results_df[row_idx, key] = value

    # Save the updated dataframe back to PolarsDB
    logging.info(f"Saving re-integrated results to sheet: {new_sheet_name}")
    plDB.insert_table(new_sheet_name, results_df)

    # Save configuration parameters
    if not plDB.has_table("Parameters"):
        plDB.create_table("Parameters", {"Parameter": pl.Utf8, "Value": pl.Utf8})

    plDB.insert_row("Parameters", {"Parameter": "FPREINT_ppm", "Value": str(ppm)})
    plDB.insert_row("Parameters", {"Parameter": "FPREINT_maxRTShift", "Value": str(maxRTShift)})
    plDB.insert_row("Parameters", {"Parameter": "FPREINT_scales", "Value": str(scales)})
    plDB.insert_row("Parameters", {"Parameter": "FPREINT_reintegrateIntensityCutoff", "Value": str(reintegrateIntensityCutoff)})
    plDB.insert_row("Parameters", {"Parameter": "FPREINT_positiveScanEvent", "Value": str(positiveScanEvent)})
    plDB.insert_row("Parameters", {"Parameter": "FPREINT_negativeScanEvent", "Value": str(negativeScanEvent)})

    # Commit and close
    plDB.commit()
    plDB.close()

    elapsed = (time.time() - start) / 60.0
    hours = ""
    if elapsed >= 60.0:
        hours = f"{int(elapsed // 60)} hours "
    mins = f"{elapsed % 60.0:.2f} mins"
    logging.info(f"Results saved after {hours}{mins}. Done")

    # Clean up multiprocessing
    p.close()
    try:
        p.terminate()
    except:
        pass
    p.join()

    if selfObj is not None and selfObj.terminateJobs:
        return


def interruptReIntegrationProcessing(pool, selfObj):
    """
    Interrupt handler for re-integration processing.

    Args:
        pool: Multiprocessing pool
        selfObj: Reference to main object

    Returns:
        True if processing should be interrupted, False otherwise
    """
    try:
        from PySide6 import QtWidgets, QtCore

        if (
            QtWidgets.QMessageBox.question(
                selfObj,
                "MetExtract",
                "Are you sure you want to cancel?",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            )
            == QtWidgets.QMessageBox.Yes
        ):
            pool.close()
            pool.terminate()
            pool.join()

            selfObj.terminateJobs = True
            logging.info("Re-integration processing stopped by user")

            return True
        else:
            return False
    except:
        # If GUI not available, just return False
        return False
