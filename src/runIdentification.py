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
import base64
import functools
import logging
import os
import platform
import pprint
import re
import time
from copy import copy
from math import floor
from pickle import loads, dumps
from sqlite3 import *

try:
    import rpy2

    R_AVAILABLE = True
except ImportError:
    R_AVAILABLE = False
    rpy2 = None
try:
    import rpy2.robjects as ro

    R_AVAILABLE = True
except ImportError:
    R_AVAILABLE = False
    ro = None

import numpy as np

from . import HCA_general
from .utils import USEGRADIENTDESCENDPEAKPICKING
from .utils import getNormRatio, CallBackMethod, get_main_dir

pp = pprint.PrettyPrinter(indent=1)

from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.charts.lineplots import ScatterPlot
from reportlab.lib import pagesizes
from reportlab.lib import colors
from reportlab.lib.colors import Color
from reportlab.lib.units import mm
from reportlab.platypus import Paragraph, Table
from reportlab.platypus.flowables import Flowable
from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF
from reportlab.lib.styles import getSampleStyleSheet


# helper function for ReportLab
def coord(x, y):
    return x, y


# helper function for ReportLab
def noLabel(a):
    return ""


# helper class for ReportLab
class TTR(Flowable):  # TableTextRotate
    """Rotates a tex in a table cell."""

    def __init__(self, text):
        Flowable.__init__(self)
        self.text = text

    def draw(self):
        canvas = self.canv
        canvas.rotate(90)
        canvas.drawString(0, -1, self.text)


from .Chromatogram import Chromatogram
from .runIdentification_matchPartners import matchPartners
from .formulaTools import formulaTools, getIsotopeMass
from .utils import (
    getAtomAdd,
    mean,
    weightedMean,
    sd,
    weightedSd,
    smoothDataSeries,
    SQLInsert,
    SQLSelectAsObject,
)
from .MZHCA import HierarchicalClustering, cutTreeSized
from .chromPeakPicking.MassSpecWavelet import MassSpecWavelet
from . import Baseline
from .utils import corr, getSubGraphs, ChromPeakPair, Bunch
from .SGR import SGRGenerator
from .mePyGuis.TracerEdit import ConfiguredTracer
from . import exportAsFeatureML
from .chromPeakPicking.CommonChromPeak import getPeakStats
import numpy as np
import scipy


def getDBSuffix():
    return ".identified.sqlite"


# returns the abbreviation for a given element
# e.g. element="Carbon" --> return="C"
def getShortName(element):
    fT = formulaTools()
    for i in fT.elemDetails:
        d = fT.elemDetails[i]
        if d[0] == element:
            return d[1]
    return ""


# counts how ofter an entry i is listed in the vector x and returns it as a dictionary
# e.g. x=[1,2,1,3,4,4,5,4] --> return={1:2, 2:1, 3:1, 4:3, 5:1}
def countEntries(x):
    ret = {}
    for i in x:
        if not (i in ret):
            ret[i] = 0
        ret[i] = ret[i] + 1
    return ret


def getCorrelationShifted(xic1, xic2, peakStartMin, peakEndMin, rtShiftsMin=None):
    if rtShiftsMin is None:
        raise RuntimeError("parameter rtShiftsMin must be an array of retention time shifts in minutes, e.g., [-3.0, -2.6, -2.2, ... 2.2, 2.6, 3.0]")

    if type(rtShiftsMin) != list and type(rtShiftsMin) != np.ndarray:
        raise RuntimeError("parameter rtShiftsMin must be an array/numpy ndarray of retention time shifts in minutes, e.g., [-3.0, -2.6, -2.2, ... 2.2, 2.6, 3.0]")

    assert xic1.shape[0] == xic2.shape[0]

    corrs = []
    nonShiftCor = correlationFor(xic1, xic2, peakStartMin, peakEndMin)

    for rtShiftMini, rtShiftMin in enumerate(rtShiftsMin):
        ## shift second EIC by current rtShiftMin
        xic2_ = np.copy(xic2)
        xic2_[:, 0] = xic2_[:, 0] + rtShiftMin

        ## get new rts
        newRTs = np.sort(np.concatenate((xic1[:, 0], xic2_[:, 0])))

        xic1_ = interpolateToNewRTs(np.copy(xic1), np.arange(np.min(newRTs), np.max(newRTs) + 0.0001, 0.01))
        xic2_ = interpolateToNewRTs(xic2_, np.arange(np.min(newRTs), np.max(newRTs) + 0.0001, 0.01))
        ## calculate correlation using peak 1
        c = correlationFor(xic1_, xic2_, peakStartMin, peakEndMin)
        corrs.append(c)

    return nonShiftCor, np.max(corrs), rtShiftsMin[np.argmax(corrs)]


def interpolateToNewRTs(xic, newRTs):
    newRTs = np.unique(newRTs)

    nxic = np.zeros((len(newRTs), 2), dtype=xic.dtype)
    for curPosN, rt in enumerate(newRTs):
        nxic[curPosN, 0] = rt

        if rt <= xic[0, 0]:
            nxic[curPosN, 1] = xic[0, 1]

        elif rt >= xic[xic.shape[0] - 1, 0]:
            nxic[curPosN, 1] = xic[xic.shape[0] - 1, 1]
        else:
            bestNextRtInd = np.argmin(np.abs(rt - xic[:, 0]))

            if xic[bestNextRtInd, 0] == rt:
                nxic[curPosN, :] = xic[bestNextRtInd, :]

            else:
                leftInd = None
                rightInd = None

                if xic[bestNextRtInd, 0] < rt:
                    leftInd = bestNextRtInd
                    rightInd = bestNextRtInd + 1

                elif xic[bestNextRtInd, 0] > rt:
                    leftInd = bestNextRtInd - 1
                    rightInd = bestNextRtInd

                assert leftInd is not None, "Branch should not exists..."

                nxic[curPosN, 1] = xic[leftInd, 1] + (xic[rightInd, 1] - xic[leftInd, 1]) / (xic[rightInd, 0] - xic[leftInd, 0]) * (rt - xic[leftInd, 0])

    return np.copy(nxic)


def smooth(abundances):
    algorithm = "rollingaverage"
    if algorithm.lower() == "rollingaverage":
        kernel = np.ones(kernel_size) / kernel_size
        data_convolved = np.convolve(abundances, kernel, mode="same")
        return data_convolved

    elif algorithm.lower() == "savgol":
        kernel_size = self.smoothingkernel_size
        window_length = self.smoothingwindow_length
        polynom_degree = self.smoothingpolynom_degree
        return scipy.signal.savgol_filter(abundances, window_length, polynom_degree)

    else:
        raise RuntimeError("Unknown algorithm '%s' for smoothing" % (algorithm))


def correlationFor(xic1, xic2, peakStartMin, peakEndMin):
    algorithm = "pearson"

    leftInd = np.argmin(np.abs(xic1[:, 0] - peakStartMin))
    rightInd = np.argmin(np.abs(xic1[:, 0] - peakEndMin))

    if algorithm == "pearson":
        try:
            return corr(xic1[leftInd:rightInd, 1], xic2[leftInd:rightInd, 1])
        except Exception as ex:
            print("There was an error calculating the pearson correlation between the two provided XICs. these are")
            print("XIC 1")
            print(xic1)
            print("")
            print("XIC 2")
            print(xic2)
            print("")
            print("Peak start and end are", peakStartMin, peakEndMin)
            raise RuntimeError("Could not calcualte correlation: '%s'" % (ex))

    else:
        raise RuntimeError("Unknown algorithm '%s' for smoothing" % (algorithm))


peakAbundanceUseSignals = 5
peakAbundanceUseSignalsSides = int((peakAbundanceUseSignals - 1) / 2)


# This class is used as a Command for each LC-HRMS file and is called by the multiprocessing module in MExtract.py
# during the processing of each individual LC-HRMS data files
class RunIdentification:
    # Constructor of the class, which stores the processing parameters
    # writeMZXML: 0001: 12C  0010: 12C-Iso  0100: 13C-Iso  1000: 13C
    def __init__(
        self,
        file,
        exOperator="",
        exExperimentID="",
        exComments="",
        exExperimentName="",
        writePDF=False,
        writeFeatureML=False,
        writeTSV=False,
        writeMZXML=9,
        metabolisationExperiment=False,
        labellingisotopeA="12C",
        labellingisotopeB="13C",
        useCIsotopePatternValidation=False,
        xOffset=1.00335,
        minRatio=0,
        maxRatio=9999999,
        useRatio=False,
        configuredTracer=None,
        intensityThreshold=0,
        intensityCutoff=0,
        maxLoading=1,
        xCounts="",
        ppm=2.0,
        isotopicPatternCountLeft=2,
        isotopicPatternCountRight=2,
        lowAbundanceIsotopeCutoff=True,
        intensityThresholdIsotopologs=1000,
        intensityErrorN=0.25,
        intensityErrorL=0.25,
        purityN=0.99,
        purityL=0.99,
        minSpectraCount=1,
        clustPPM=8.0,
        chromPeakPPM=5.0,
        snrTh=1.0,
        scales=[1, 35],
        peakCenterError=5,
        peakScaleError=3,
        minPeakCorr=0.85,
        checkPeaksRatio=False,
        minPeaksRatio=0,
        maxPeaksRatio=99999999,
        calcIsoRatioNative=1,
        calcIsoRatioLabelled=-1,
        calcIsoRatioMoiety=1,
        startTime=2,
        stopTime=37,
        positiveScanEvent="None",
        negativeScanEvent="None",
        eicSmoothingWindow="None",
        eicSmoothingWindowSize=0,
        eicSmoothingPolynom=0,
        artificialMPshift_start=0,
        artificialMPshift_stop=0,
        correctCCount=True,
        minCorrelation=0.85,
        minCorrelationConnections=0.4,
        hAIntensityError=5.0,
        hAMinScans=3,
        adducts=[],
        elements=[],
        heteroAtoms=[],
        simplifyInSourceFragments=True,
        chromPeakFile=None,
        lock=None,
        queue=None,
        pID=-1,
        rVersion="NA",
        meVersion="NA",
    ):
        # Check if configuredTracer is None when metabolisationExperiment is True
        if metabolisationExperiment and configuredTracer is None:
            self.printMessage("Metabolisation experiment requires a configured tracer, but none was provided.", type="error")
            raise Exception("Metabolisation experiment requires a configured tracer, but none was provided.")

        ma, ea = getIsotopeMass(configuredTracer.isotopeA if metabolisationExperiment else labellingisotopeA)
        mb, eb = getIsotopeMass(configuredTracer.isotopeB if metabolisationExperiment else labellingisotopeB)
        xOffset = mb - ma if metabolisationExperiment else xOffset

        if (not metabolisationExperiment and len(ea) > 0 and ea == eb) or metabolisationExperiment or (ea == "Hydrogen" and eb == "Deuterium"):
            pass
        else:
            self.printMessage("Labelling specifications are invalid..", type="warning")
            raise Exception("Labelling specifications are invalid..")

        if positiveScanEvent == "None" and negativeScanEvent == "None":
            self.printMessage("No scan event was specified..", type="warning")
            raise Exception("No scan event was specified..")

        self.file = file
        self.writePDF = writePDF
        self.writeFeatureML = writeFeatureML
        self.writeTSV = writeTSV
        self.writeMZXML = writeMZXML

        # E. Experiment
        self.experimentOperator = exOperator
        self.experimentID = exExperimentID
        self.experimentComments = exComments
        self.experimentName = exExperimentName

        # 0. General
        self.startTime = startTime
        self.stopTime = stopTime
        self.positiveScanEvent = positiveScanEvent
        self.negativeScanEvent = negativeScanEvent

        self.metabolisationExperiment = metabolisationExperiment
        self.labellingElement = getShortName(ea)
        self.isotopeA = ma
        self.isotopeB = mb
        self.useCIsotopePatternValidation = useCIsotopePatternValidation
        self.minRatio = minRatio
        self.maxRatio = maxRatio
        self.useRatio = useRatio

        self.configuredTracer = configuredTracer
        if not self.metabolisationExperiment:
            self.configuredTracer = ConfiguredTracer(name="FML", id=0)

        # 1. Mass picking
        self.intensityThreshold = intensityThreshold
        self.intensityCutoff = intensityCutoff
        self.maxLoading = maxLoading

        self.xCounts = []
        self.xCountsString = xCounts
        a = xCounts.replace(" ", "").replace(";", ",").split(",")
        for j in a:
            if "-" in j:
                self.xCounts.extend(range(int(j[0 : j.find("-")]), int(j[j.find("-") + 1 : len(j)]) + 1))
            elif ":" in j:
                self.xCounts.extend(range(int(j[0 : j.find(":")]), int(j[j.find(":") + 1 : len(j)]) + 1))
            elif j != "":
                self.xCounts.append(int(j))
        self.xCounts = sorted(list(set(self.xCounts)))

        self.xOffset = xOffset
        self.ppm = ppm
        self.isotopicPatternCountLeft = isotopicPatternCountLeft
        self.isotopicPatternCountRight = isotopicPatternCountRight
        self.lowAbundanceIsotopeCutoff = lowAbundanceIsotopeCutoff
        self.intensityThresholdIsotopologs = intensityThresholdIsotopologs
        self.intensityErrorN = intensityErrorN
        self.intensityErrorL = intensityErrorL
        self.purityN = purityN
        self.purityL = purityL

        # 2. Results clustering
        self.minSpectraCount = max(1, minSpectraCount)
        self.clustPPM = clustPPM

        # 3. Peak detection
        self.chromPeakPPM = chromPeakPPM

        if eicSmoothingWindow == None or eicSmoothingWindow == "" or eicSmoothingWindow == "none":
            eicSmoothingWindow = "None"
        self.eicSmoothingWindow = eicSmoothingWindow
        self.eicSmoothingWindowSize = eicSmoothingWindowSize
        self.eicSmoothingPolynom = eicSmoothingPolynom
        self.artificialMPshift_start = artificialMPshift_start
        self.artificialMPshift_stop = artificialMPshift_stop
        self.snrTh = snrTh
        self.scales = scales
        self.peakCenterError = peakCenterError
        self.peakScaleError = peakScaleError
        self.minPeakCorr = minPeakCorr
        self.checkPeaksRatio = checkPeaksRatio
        self.minPeaksRatio = minPeaksRatio
        self.maxPeaksRatio = maxPeaksRatio

        self.calcIsoRatioNative = calcIsoRatioNative
        self.calcIsoRatioLabelled = calcIsoRatioLabelled
        self.calcIsoRatioMoiety = calcIsoRatioMoiety

        self.chromPeakFile = chromPeakFile
        if self.chromPeakFile is None:
            cpf = get_main_dir() + "/src/chromPeakPicking/MassSpecWaveletIdentification.r"
            self.chromPeakFile = cpf

        self.performCorrectCCount = correctCCount

        self.minCorrelation = minCorrelation
        self.minCorrelationConnections = minCorrelationConnections

        self.hAIntensityError = hAIntensityError
        self.hAMinScans = hAMinScans
        self.adducts = adducts
        self.elements = {}
        for elem in elements:
            self.elements[elem.name] = elem
        self.heteroAtoms = {}
        for heteroAtom in heteroAtoms:
            self.heteroAtoms[heteroAtom.name] = heteroAtom
        self.simplifyInSourceFragments = simplifyInSourceFragments

        # System
        self.lock = lock
        self.queue = queue
        self.pID = pID
        self.rVersion = rVersion
        self.meVersion = meVersion

    # Thread safe printing function
    def printMessage(self, message, type="info"):
        if self.lock is not None:
            self.lock.acquire()
            if type.lower() == "info":
                logging.info("   %d: %s" % (self.pID, message))
            elif type.lower() == "warning":
                logging.warning("   %d: %s" % (self.pID, message))
            elif type.lower() == "error":
                logging.error("   %d: %s" % (self.pID, message))
            else:
                logging.debug("   %d: %s" % (self.pID, message))
            self.lock.release()

    # helper function used to update the status of the current processing in the Process Dialog
    def postMessageToProgressWrapper(self, mes, val=""):
        if self.pID != -1 and self.queue is not None:
            if mes.lower() == "text":
                self.queue.put(Bunch(pid=self.pID, mes="text", val="%d: %s" % (self.pID, val)))
            elif mes == "value" or mes == "max":
                self.queue.put(Bunch(pid=self.pID, mes=mes, val=val))
            elif mes == "start" or mes == "end" or mes == "failed":
                self.queue.put(Bunch(pid=self.pID, mes=mes))

    def getMostLikelyHeteroIsotope(self, foundIsotopes):
        if len(foundIsotopes) == 0:
            return

        for iso in foundIsotopes:
            isoD = foundIsotopes[iso]
            useCn = -1
            useCnRatio = 0.0
            useCnRatioTheo = 0.0

            for cn in isoD:
                cnF = isoD[cn]

                cnRatio = mean([x[1] for x in cnF])
                cnRatioTheo = mean([x[2] for x in cnF])
                cnScans = len(cnF)

                if useCn == -1 or (abs(cnRatioTheo - cnRatio) < abs(useCnRatio - useCnRatioTheo)):
                    useCn = cn
                    useCnRatio = cnRatio
                    useCnRatioTheo = cnRatioTheo

            foundIsotopes[iso] = {useCn: isoD[useCn]}

    # store configuration used for processing the LC-HRMS file into the database
    def writeConfigurationToDB(self, conn, curs):
        curs.execute("DROP TABLE IF EXISTS config")
        curs.execute("CREATE TABLE config (id INTEGER PRIMARY KEY AUTOINCREMENT, key TEXT, value TEXT)")

        curs.execute("DROP TABLE IF EXISTS tracerConfiguration")
        curs.execute(
            "create table tracerConfiguration(id INTEGER PRIMARY KEY, name TEXT, elementCount INTEGER, natural TEXT, labelling TEXT, deltaMZ REAL, purityN REAL, purityL REAL, amountN REAL, amountL REAL, monoisotopicRatio REAL, checkRatio TEXT, lowerError REAL, higherError REAL, tracertype TEXT)"
        )

        curs.execute("DROP TABLE IF EXISTS MZs")
        curs.execute("create table MZs(id INTEGER PRIMARY KEY, tracer INTEGER, mz REAL, lmz REAL, tmz REAL, xcount INTEGER, scanid INTEGER, scantime REAL, loading INTEGER, intensity FLOAT, intensityL FLOAT, ionMode TEXT)")

        curs.execute("DROP TABLE IF EXISTS MZBins")
        curs.execute("CREATE TABLE MZBins(id INTEGER PRIMARY KEY, mz REAL, ionMode TEXT)")

        curs.execute("DROP TABLE IF EXISTS MZBinsKids")
        curs.execute("CREATE TABLE MZBinsKids(mzbinID INTEGER, mzID INTEGER)")

        curs.execute("DROP TABLE IF EXISTS XICs")
        curs.execute(
            "CREATE TABLE XICs(id INTEGER PRIMARY KEY, avgmz REAL, xcount INTEGER, loading INTEGER, polarity TEXT, xic TEXT, xic_smoothed TEXT, xic_baseline TEXT, xicL TEXT, xicL_smoothed TEXT, xicL_baseline TEXT, xicfirstiso TEXT, xicLfirstiso TEXT, xicLfirstisoconjugate TEXT, mzs TEXT, mzsL TEXT, mzsfirstiso TEXT, mzsLfirstiso TEXT, mzsLfirstisoconjugate TEXT, times TEXT, scanCount INTEGER, allPeaks TEXT)"
        )

        curs.execute("DROP TABLE IF EXISTS tics")
        curs.execute("CREATE TABLE tics(id INTEGER PRIMARY KEY, loading INTEGER, scanevent TEXT, times TEXT, intensities TEXT)")

        curs.execute("DROP TABLE IF EXISTS chromPeaks")
        curs.execute(
            "CREATE TABLE chromPeaks(id INTEGER PRIMARY KEY, tracer INTEGER, eicID INTEGER, NPeakCenter INTEGER, NPeakCenterMin REAL, NPeakScale FLOAT, NSNR REAL, NPeakArea REAL, NPeakAbundance REAL, mz REAL, lmz REAL, tmz REAL, xcount INTEGER, xcountId INTEGER, LPeakCenter INTEGER, LPeakCenterMin REAL, LPeakScale FLOAT, LSNR REAL, LPeakArea REAL, LPeakAbundance REAL, Loading INTEGER, peaksCorr FLOAT, heteroAtoms TEXT, NBorderLeft INTEGER, NBorderRight INTEGER, LBorderLeft INTEGER, LBorderRight INTEGER, adducts TEXT, heteroAtomsFeaturePairs TEXT, massSpectrumID INTEGER, ionMode TEXT, assignedMZs INTEGER, fDesc TEXT, peaksRatio FLOAT, peaksRatioMp1 FLOAT, peaksRatioMPm1 FLOAT, isotopesRatios TEXT, mzDiffErrors TEXT, peakType TEXT, assignedName TEXT, correlationsToOthers TEXT, comments TEXT, artificialEICLShift INTEGER, new_FWHM_M FLOAT, new_FWHM_Mp FLOAT, new_Area_M FLOAT, new_Area_Mp FLOAT, new_SNR_M FLOAT, new_SNR_Mp FLOAT)"
        )

        curs.execute("drop table if exists allChromPeaks")
        curs.execute(
            "CREATE TABLE allChromPeaks(id INTEGER PRIMARY KEY, tracer INTEGER, eicID INTEGER, NPeakCenter INTEGER, NPeakCenterMin REAL, NPeakScale FLOAT, NSNR REAL, NPeakArea REAL, NPeakAbundance REAL, mz REAL, lmz REAL, tmz REAL, xcount INTEGER, xcountId INTEGER, LPeakCenter INTEGER, LPeakCenterMin REAL, LPeakScale FLOAT, LSNR REAL, LPeakArea REAL, LPeakAbundance REAL, Loading INTEGER, peaksCorr FLOAT, heteroAtoms TEXT, NBorderLeft INTEGER, NBorderRight INTEGER, LBorderLeft INTEGER, LBorderRight INTEGER, adducts TEXT, heteroAtomsFeaturePairs TEXT, ionMode TEXT, assignedMZs INTEGER, fDesc TEXT, peaksRatio FLOAT, peaksRatioMp1 FLOAT, peaksRatioMPm1 FLOAT, isotopesRatios TEXT, mzDiffErrors TEXT, peakType TEXT, assignedName TEXT, comment TEXT, comments TEXT, artificialEICLShift INTEGER, new_FWHM_M FLOAT, new_FWHM_Mp FLOAT, new_Area_M FLOAT, new_Area_Mp FLOAT, new_SNR_M FLOAT, new_SNR_Mp FLOAT)"
        )

        curs.execute("DROP TABLE IF EXISTS featureGroups")
        curs.execute("CREATE TABLE featureGroups (id INTEGER PRIMARY KEY, featureName TEXT, tracer INTEGER)")

        curs.execute("DROP TABLE IF EXISTS featureGroupFeatures")
        curs.execute("CREATE TABLE featureGroupFeatures (id INTEGER PRIMARY KEY, fID INTEGER, fDesc TEXT, fGroupID INTEGER)")

        curs.execute("DROP TABLE IF EXISTS featurefeatures")
        curs.execute("CREATE TABLE featurefeatures (fID1 INTEGER, fID2 INTEGER, corr FLOAT, silRatioValue FLOAT, desc1 TEXT, desc2 TEXT, add1 TEXT, add2 TEXT)")

        curs.execute("DROP TABLE IF EXISTS massspectrum")
        curs.execute("CREATE TABLE massspectrum (mID INTEGER, fgID INTEGER, time FLOAT, mzs TEXT, intensities TEXT, ionMode TEXT)")

        curs.execute("DROP TABLE IF EXISTS stats")
        curs.execute("CREATE TABLE stats (key TEXT, value TEXT)")

        SQLInsert(curs, "config", key="MetExtractVersion", value=self.meVersion)
        SQLInsert(curs, "config", key="RVersion", value=self.rVersion)

        SQLInsert(curs, "config", key="ExperimentName", value=self.experimentName)
        SQLInsert(curs, "config", key="ExperimentOperator", value=self.experimentOperator)
        SQLInsert(curs, "config", key="ExperimentID", value=self.experimentID)
        SQLInsert(curs, "config", key="ExperimentComments", value=self.experimentComments)

        SQLInsert(curs, "config", key="labellingElement", value=self.labellingElement)
        SQLInsert(curs, "config", key="isotopeA", value=self.isotopeA)
        SQLInsert(curs, "config", key="isotopeB", value=self.isotopeB)
        SQLInsert(
            curs,
            "config",
            key="useCValidation",
            value=self.useCIsotopePatternValidation,
        )
        SQLInsert(curs, "config", key="minRatio", value=self.minRatio)
        SQLInsert(curs, "config", key="maxRatio", value=self.maxRatio)
        SQLInsert(curs, "config", key="useRatio", value=str(self.useRatio))
        SQLInsert(
            curs,
            "config",
            key="metabolisationExperiment",
            value=str(self.metabolisationExperiment),
        )
        SQLInsert(
            curs,
            "config",
            key="configuredTracer",
            value=base64.b64encode(dumps(self.configuredTracer)).decode("utf-8"),
        )
        SQLInsert(curs, "config", key="startTime", value=self.startTime)
        SQLInsert(curs, "config", key="stopTime", value=self.stopTime)
        SQLInsert(curs, "config", key="positiveScanEvent", value=self.positiveScanEvent)
        SQLInsert(curs, "config", key="negativeScanEvent", value=self.negativeScanEvent)
        SQLInsert(curs, "config", key="intensityThreshold", value=self.intensityThreshold)
        SQLInsert(curs, "config", key="intensityCutoff", value=self.intensityCutoff)
        SQLInsert(curs, "config", key="maxLoading", value=self.maxLoading)
        SQLInsert(curs, "config", key="xCounts", value=self.xCountsString)
        SQLInsert(curs, "config", key="xOffset", value=self.xOffset)
        SQLInsert(curs, "config", key="ppm", value=self.ppm)
        SQLInsert(
            curs,
            "config",
            key="isotopicPatternCountLeft",
            value=self.isotopicPatternCountLeft,
        )
        SQLInsert(
            curs,
            "config",
            key="isotopicPatternCountRight",
            value=self.isotopicPatternCountRight,
        )
        SQLInsert(
            curs,
            "config",
            key="lowAbundanceIsotopeCutoff",
            value=str(self.lowAbundanceIsotopeCutoff),
        )
        SQLInsert(
            curs,
            "config",
            key="intensityThresholdIsotopologs",
            value=str(self.intensityThresholdIsotopologs),
        )
        SQLInsert(curs, "config", key="intensityErrorN", value=self.intensityErrorN)
        SQLInsert(curs, "config", key="intensityErrorL", value=self.intensityErrorL)
        SQLInsert(curs, "config", key="purityN", value=self.purityN)
        SQLInsert(curs, "config", key="purityL", value=self.purityL)
        SQLInsert(curs, "config", key="clustPPM", value=self.clustPPM)
        SQLInsert(curs, "config", key="chromPeakPPM", value=self.chromPeakPPM)
        SQLInsert(curs, "config", key="eicSmoothing", value=self.eicSmoothingWindow)
        SQLInsert(
            curs,
            "config",
            key="eicSmoothingWindowSize",
            value=self.eicSmoothingWindowSize,
        )
        SQLInsert(curs, "config", key="eicSmoothingPolynom", value=self.eicSmoothingPolynom)
        SQLInsert(
            curs,
            "config",
            key="artificialMPshift_start",
            value=self.artificialMPshift_start,
        )
        SQLInsert(
            curs,
            "config",
            key="artificialMPshift_stop",
            value=self.artificialMPshift_stop,
        )
        SQLInsert(curs, "config", key="snrTh", value=self.snrTh)
        SQLInsert(
            curs,
            "config",
            key="scales",
            value=base64.b64encode(dumps(self.scales)).decode("utf-8"),
        )
        SQLInsert(
            curs,
            "config",
            key="peakAbundanceCriteria",
            value="Center +- %d signals (%d total)" % (peakAbundanceUseSignalsSides, peakAbundanceUseSignals),
        )
        SQLInsert(curs, "config", key="peakCenterError", value=self.peakCenterError)
        SQLInsert(curs, "config", key="peakScaleError", value=self.peakScaleError)
        SQLInsert(curs, "config", key="minPeakCorr", value=self.minPeakCorr)
        SQLInsert(curs, "config", key="checkPeaksRatio", value=str(self.checkPeaksRatio))
        SQLInsert(curs, "config", key="minPeaksRatio", value=self.minPeaksRatio)
        SQLInsert(curs, "config", key="maxPeaksRatio", value=self.maxPeaksRatio)
        SQLInsert(curs, "config", key="calcIsoRatioNative", value=self.calcIsoRatioNative)
        SQLInsert(curs, "config", key="calcIsoRatioLabelled", value=self.calcIsoRatioLabelled)
        SQLInsert(curs, "config", key="calcIsoRatioMoiety", value=self.calcIsoRatioMoiety)
        SQLInsert(curs, "config", key="minSpectraCount", value=self.minSpectraCount)
        SQLInsert(
            curs,
            "config",
            key="configuredHeteroAtoms",
            value=base64.b64encode(dumps(self.heteroAtoms)).decode("utf-8"),
        )
        SQLInsert(curs, "config", key="haIntensityError", value=self.hAIntensityError)
        SQLInsert(curs, "config", key="haMinScans", value=self.hAMinScans)
        SQLInsert(curs, "config", key="minCorrelation", value=self.minCorrelation)
        SQLInsert(
            curs,
            "config",
            key="minCorrelationConnections",
            value=self.minCorrelationConnections,
        )
        SQLInsert(
            curs,
            "config",
            key="adducts",
            value=base64.b64encode(dumps(self.adducts)).decode("utf-8"),
        )
        SQLInsert(
            curs,
            "config",
            key="elements",
            value=base64.b64encode(dumps(self.elements)).decode("utf-8"),
        )
        SQLInsert(
            curs,
            "config",
            key="simplifyInSourceFragments",
            value=str(self.simplifyInSourceFragments),
        )

        import uuid
        import platform
        import datetime

        self.processingUUID = "%s_%s_%s" % (
            str(uuid.uuid1()),
            str(platform.node()),
            str(datetime.datetime.now()),
        )
        SQLInsert(
            curs,
            "config",
            key="processingUUID_ext",
            value=base64.b64encode(str(self.processingUUID).encode("utf-8")),
        )

        i = 1
        if self.metabolisationExperiment:
            tracer = self.configuredTracer
            tracer.id = i
            SQLInsert(
                curs,
                "tracerConfiguration",
                id=tracer.id,
                name=tracer.name,
                elementCount=tracer.elementCount,
                natural=tracer.isotopeA,
                labelling=tracer.isotopeB,
                deltaMZ=getIsotopeMass(tracer.isotopeB)[0] - getIsotopeMass(tracer.isotopeA)[0],
                purityN=tracer.enrichmentA,
                purityL=tracer.enrichmentB,
                amountN=tracer.amountA,
                amountL=tracer.amountB,
                monoisotopicRatio=tracer.monoisotopicRatio,
                checkRatio=str(tracer.checkRatio),
                lowerError=tracer.maxRelNegBias,
                higherError=tracer.maxRelPosBias,
                tracertype=tracer.tracerType,
            )

        else:
            # ConfiguredTracer(name="Full metabolome labeling experiment", id=0)
            SQLInsert(curs, "tracerConfiguration", id=0, name="FLE")
        conn.commit()

    def parseMzXMLFile(self):
        mzxml = Chromatogram()
        mzxml.parse_file(self.file, intensityCutoff=self.intensityCutoff)
        return mzxml

    # creates a new PDF page which contains the used data processing parameters
    def writeSettingsToPDF(self, pdf):
        currentHeight = 800

        pdf.drawString(50, currentHeight, "Experiment name")
        pdf.drawString(240, currentHeight, self.experimentName)
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2)
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "Operator")
        pdf.drawString(240, currentHeight, self.experimentOperator)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "ID")
        pdf.drawString(240, currentHeight, self.experimentID)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Comments")
        currentHeight -= 15

        p = Paragraph(self.experimentComments, style=getSampleStyleSheet()["Normal"])
        w, h = p.wrap(450, 120)
        p.wrapOn(pdf, 450, 120)
        p.drawOn(pdf, 60, currentHeight - h + 10)
        pdf.showPage()

        currentHeight = 800

        pdf.drawString(50, currentHeight, "Parameters")
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2)
        currentHeight -= 20

        if self.metabolisationExperiment:
            pdf.drawString(70, currentHeight, "Metabolisation Experiment")
            pdf.drawString(240, currentHeight, "1 tracer")
            currentHeight -= 25

        else:
            pdf.drawString(50, currentHeight, "Full Metabolome Experiment")
            currentHeight -= 15
            pdf.drawString(50, currentHeight, "Labelling")

            pdf.drawString(
                240,
                currentHeight,
                "%s (%s %s)" % (self.labellingElement, self.isotopeA, self.isotopeB),
            )
            currentHeight -= 15

            pdf.drawString(
                50,
                currentHeight,
                "Abundance %d%s" % (self.isotopeA, str(self.labellingElement)),
            )
            pdf.drawString(240, currentHeight, "%.2f%%" % self.purityN)
            currentHeight -= 15

            pdf.drawString(
                70,
                currentHeight,
                "Abundance %d%s" % (self.isotopeB, str(self.labellingElement)),
            )
            pdf.drawString(240, currentHeight, "%.2f%%" % self.purityL)
            currentHeight -= 35

        pdf.drawString(50, currentHeight, "M/Z Picking")
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2)
        currentHeight -= 15
        pdf.drawString(70, currentHeight, "Scan event(s)")
        currentHeight -= 15

        if self.positiveScanEvent != "None":
            pdf.drawString(90, currentHeight, "Positive mode")
            pdf.drawString(240, currentHeight, "%s" % self.positiveScanEvent)
            currentHeight -= 15
        if self.negativeScanEvent != "None":
            pdf.drawString(90, currentHeight, "Negative mode")
            pdf.drawString(220, currentHeight, "%s" % self.negativeScanEvent)
            currentHeight -= 15

        pdf.drawString(70, currentHeight, "Intensity threshold")
        pdf.drawString(
            240,
            currentHeight,
            "%d%s"
            % (
                self.intensityThreshold,
                " (Low abundance cutoff for isotopologues: %.0f)" % self.intensityThresholdIsotopologs if self.lowAbundanceIsotopeCutoff else "",
            ),
        )
        currentHeight -= 15
        pdf.drawString(70, currentHeight, "Intensity cutoff")
        pdf.drawString(240, currentHeight, "%d" % (self.intensityCutoff))
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Time")
        pdf.drawString(240, currentHeight, "%.1f-%.1f min" % (self.startTime, self.stopTime))
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Atom count")
        pdf.drawString(240, currentHeight, "%s" % (",".join(str(x) for x in self.xCountsString)))
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Max. charge")
        pdf.drawString(240, currentHeight, "%d" % self.maxLoading)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Max. mass deviation ")
        pdf.drawString(240, currentHeight, "%.2f" % self.ppm)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Isotopic pattern count")
        pdf.drawString(
            240,
            currentHeight,
            "Native: %d Labelled: %d" % (self.isotopicPatternCountLeft, self.isotopicPatternCountRight),
        )
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Intensity abundance error")
        pdf.drawString(
            240,
            currentHeight,
            "Native: %.1f%% Labelled: %.1f%%" % (self.intensityErrorN * 100.0, self.intensityErrorL * 100.0),
        )
        currentHeight -= 35

        pdf.drawString(50, currentHeight, "Post Processing")
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2)
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "Clustering PPM")
        pdf.drawString(240, currentHeight, "%.2f" % self.clustPPM)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "EIC ppm")
        pdf.drawString(240, currentHeight, "%.2f" % self.chromPeakPPM)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "EIC Smoothing Window")
        pdf.drawString(240, currentHeight, self.eicSmoothingWindow)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "M' artificial shift")
        pdf.drawString(
            240,
            currentHeight,
            "%d - %d" % (self.artificialMPshift_start, self.artificialMPshift_stop),
        )
        currentHeight -= 15

        if self.eicSmoothingWindow.lower() != "none":
            pdf.drawString(70, currentHeight, "EIC Smoothing Window Size")
            pdf.drawString(240, currentHeight, str(self.eicSmoothingWindowSize))
            currentHeight -= 15

        if self.eicSmoothingWindow.lower() == "SavitzkyGolay":
            pdf.drawString(70, currentHeight, "EIC Smoothing Polynom")
            pdf.drawString(240, currentHeight, str(self.eicSmoothingPolynom))
            currentHeight -= 15

        pdf.drawString(70, currentHeight, "Scales")
        pdf.drawString(240, currentHeight, "%d - %d" % (self.scales[0], self.scales[1]))
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Peak matching")
        # pdf.drawString(240, currentHeight, "Center: %d Scales: %d Min. Corr: %.2f" % (
        #    self.peakCenterError, self.peakScaleError, self.minPeakCorr));
        pdf.drawString(
            240,
            currentHeight,
            "Center: %d, Min. Corr: %.2f" % (self.peakCenterError, self.minPeakCorr),
        )
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Min. Corr:")
        pdf.drawString(240, currentHeight, "%.2f" % self.minPeakCorr)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Min. scans")
        pdf.drawString(240, currentHeight, "%d" % self.minSpectraCount)
        currentHeight -= 35

        pdf.drawString(50, currentHeight, "Hetero atom annotation")
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2)
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "Max. intensity error")
        pdf.drawString(240, currentHeight, "%.1f%%" % (self.hAIntensityError * 100.0))
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Min. scans")
        pdf.drawString(240, currentHeight, "%d" % self.hAMinScans)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Hetero atom isotopes")
        p = Paragraph(
            str(
                ", ".join(
                    [
                        "%s (m: %.4f, rel. ab. %.1f%%, min: %d, max: %d)"
                        % (
                            pIso,
                            self.heteroAtoms[pIso].mzOffset,
                            self.heteroAtoms[pIso].relativeAbundance * 100.0,
                            self.heteroAtoms[pIso].minCount,
                            self.heteroAtoms[pIso].maxCount,
                        )
                        for pIso in self.heteroAtoms
                    ]
                )
            ),
            style=getSampleStyleSheet()["Normal"],
        )
        w, h = p.wrap(300, 60)
        p.wrapOn(pdf, 300, 60)
        p.drawOn(pdf, 240, currentHeight - h + 10)
        currentHeight -= max(35, h + 10)

        pdf.drawString(50, currentHeight, "Non-targeted feature grouping")
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2)
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "Min. correlation")
        pdf.drawString(240, currentHeight, "%.2f" % self.minCorrelation)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Min. correlation connections")
        pdf.drawString(240, currentHeight, "%.2f" % self.minCorrelationConnections)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Adducts")
        currentHeight -= 15
        p = Paragraph(
            str(", ".join(["%s (m/z: %.4f, z:%d%s)" % (ad.name, ad.mzoffset, ad.charge, ad.polarity) for ad in self.adducts])),
            style=getSampleStyleSheet()["Normal"],
        )
        w, h = p.wrap(460, 60)
        p.wrapOn(pdf, 460, 60)
        p.drawOn(pdf, 80, currentHeight - h + 10)
        currentHeight -= max(20, h + 10)

        pdf.drawString(70, currentHeight, "Elements")
        currentHeight -= 15
        p = Paragraph(
            str(", ".join(["%s (m: %.4f)" % (el, self.elements[el].weight) for el in self.elements.keys()])),
            style=getSampleStyleSheet()["Normal"],
        )
        w, h = p.wrap(460, 60)
        p.wrapOn(pdf, 460, 60)
        p.drawOn(pdf, 80, currentHeight - h + 10)
        currentHeight -= max(35, h + 10)

        pdf.drawString(
            70,
            currentHeight,
            "Simplify in-source fragments: " + str(self.simplifyInSourceFragments),
        )
        currentHeight -= 15

        pdf.drawString(50, currentHeight, "Software")
        pdf.line(50, currentHeight - 2, 540, currentHeight - 2)
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "MetExtract version")
        pdf.drawString(240, currentHeight, str(self.meVersion))
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "R-Version")
        pdf.drawString(240, currentHeight, str(self.rVersion))
        currentHeight -= 15

        p = Paragraph("UUID_ext: " + self.processingUUID, style=getSampleStyleSheet()["Normal"])
        w, h = p.wrap(450, 120)
        p.wrapOn(pdf, 450, 120)
        p.drawOn(pdf, 60, currentHeight - h + 10)
        currentHeight -= 15

        pdf.showPage()

    # creates a new PDF page which contains the tracer used in this experiment
    def writeCurrentTracerToPDF(self, pdf, tracer):
        currentHeight = 800
        pdf.drawString(50, currentHeight, "Tracer: %s" % tracer.name)
        pdf.line(50, currentHeight - 2, 540, currentHeight - 1)
        currentHeight -= 20

        pdf.drawString(70, currentHeight, "Labelling")
        pdf.drawString(240, currentHeight, "%s, %s" % (tracer.isotopeA, tracer.isotopeB))
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Delta m/z")
        pdf.drawString(240, currentHeight, "%.5f" % self.xOffset)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Element count")
        pdf.drawString(240, currentHeight, "%d" % tracer.elementCount)
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Purity %s" % tracer.isotopeA)
        pdf.drawString(240, currentHeight, "%.2f%%" % (tracer.enrichmentA * 100.0))
        currentHeight -= 15

        pdf.drawString(70, currentHeight, "Purity %s" % tracer.isotopeB)
        pdf.drawString(240, currentHeight, "%.2f%%" % (tracer.enrichmentB * 100.0))
        currentHeight -= 15

        if tracer.checkRatio:
            pdf.drawString(
                70,
                currentHeight,
                "Monoisotopic ratio (%s:%s)" % (tracer.isotopeA, tracer.isotopeB),
            )
            pdf.drawString(
                240,
                currentHeight,
                "%.5f (%.2f:%.2f [v/v])" % (tracer.monoisotopicRatio, tracer.amountA, tracer.amountB),
            )
            currentHeight -= 15

            pdf.drawString(70, currentHeight, "Max. intensity error")
            pdf.drawString(
                240,
                currentHeight,
                "%.1f%%, %.1f%% (%.3f - %.3f)"
                % (
                    tracer.maxRelNegBias * 100.0,
                    tracer.maxRelPosBias * 100.0,
                    tracer.monoisotopicRatio * tracer.maxRelNegBias,
                    tracer.monoisotopicRatio * tracer.maxRelPosBias,
                ),
            )
            currentHeight -= 15

        else:
            pdf.drawString(70, currentHeight, "Ratio was not checked")
            currentHeight -= 15

        pdf.showPage()

    # data processing step 1: searches each mass spectrum for isotope patterns of native and highly isotope enriched
    # metabolite ions. The actual calculation and processing of the data is performed in the file runIdentification_matchPartners.py.
    # The positive and negative ionisation modes are processed separately.
    def findSignalPairs(self, curProgress, mzxml, tracer, reportFunction=None):
        mzs = []
        posFound = 0
        negFound = 0

        checkRatio = False
        minRatio = 0.0
        maxRatio = 0.0

        if self.metabolisationExperiment:
            checkRatio = tracer.checkRatio
            minRatio = tracer.monoisotopicRatio * tracer.maxRelNegBias
            maxRatio = tracer.monoisotopicRatio * tracer.maxRelPosBias
        else:
            checkRatio = self.useRatio
            minRatio = self.minRatio
            maxRatio = self.maxRatio

        def reportFunctionHelper(curVal, text):
            reportFunction(curVal, text)

        if self.positiveScanEvent != "None":
            if self.negativeScanEvent != "None":

                def reportFunctionHelper(curVal, text):
                    reportFunction(curVal / 2, text)

            p = matchPartners(
                mzXMLData=mzxml,
                forFile=self.file,
                labellingElement=self.labellingElement,
                useCIsotopePatternValidation=self.useCIsotopePatternValidation,
                intensityThres=self.intensityThreshold,
                isotopologIntensityThres=self.intensityThresholdIsotopologs,
                maxLoading=self.maxLoading,
                xCounts=self.xCounts,
                xOffset=self.xOffset,
                ppm=self.ppm,
                intensityErrorN=self.intensityErrorN,
                intensityErrorL=self.intensityErrorL,
                purityN=tracer.enrichmentA if self.metabolisationExperiment else self.purityN,
                purityL=tracer.enrichmentB if self.metabolisationExperiment else self.purityL,
                startTime=self.startTime,
                stopTime=self.stopTime,
                filterLine=self.positiveScanEvent,
                ionMode="+",
                peakCountLeft=self.isotopicPatternCountLeft,
                peakCountRight=self.isotopicPatternCountRight,
                lowAbundanceIsotopeCutoff=self.lowAbundanceIsotopeCutoff,
                metabolisationExperiment=self.metabolisationExperiment,
                checkRatio=checkRatio,
                minRatio=minRatio,
                maxRatio=maxRatio,
                reportFunction=reportFunctionHelper,
            )
            posFound = len(p)
            mzs.extend(p)

        def reportFunctionHelper(curVal, text):
            reportFunction(curVal, text)

        if self.negativeScanEvent != "None":
            if self.positiveScanEvent != "None":

                def reportFunctionHelper(curVal, text):
                    reportFunction(0.5 + curVal / 2, text)

            n = matchPartners(
                mzXMLData=mzxml,
                forFile=self.file,
                labellingElement=self.labellingElement,
                useCIsotopePatternValidation=self.useCIsotopePatternValidation,
                intensityThres=self.intensityThreshold,
                isotopologIntensityThres=self.intensityThresholdIsotopologs,
                maxLoading=self.maxLoading,
                xCounts=self.xCounts,
                xOffset=self.xOffset,
                ppm=self.ppm,
                intensityErrorN=self.intensityErrorN,
                intensityErrorL=self.intensityErrorL,
                purityN=tracer.enrichmentA if self.metabolisationExperiment > 1 else self.purityN,
                purityL=tracer.enrichmentB if self.metabolisationExperiment > 1 else self.purityL,
                startTime=self.startTime,
                stopTime=self.stopTime,
                filterLine=self.negativeScanEvent,
                ionMode="-",
                peakCountLeft=self.isotopicPatternCountLeft,
                peakCountRight=self.isotopicPatternCountRight,
                lowAbundanceIsotopeCutoff=self.lowAbundanceIsotopeCutoff,
                metabolisationExperiment=self.metabolisationExperiment,
                checkRatio=checkRatio,
                minRatio=minRatio,
                maxRatio=maxRatio,
                reportFunction=reportFunctionHelper,
            )
            negFound = len(n)
            mzs.extend(n)

        return mzs, negFound, posFound

    # store detected signal pairs (1st data processing step) in the database
    def writeSignalPairsToDB(self, mzs, mzxml, tracerID):
        conn = connect(self.file + getDBSuffix(), isolation_level="DEFERRED")
        conn.execute("""PRAGMA synchronous = OFF""")
        conn.execute("""PRAGMA journal_mode = OFF""")
        curs = conn.cursor()

        for mz in mzs:
            mz.id = self.curMZId
            mz.tid = tracerID

            scanEvent = ""
            if mz.ionMode == "+":
                scanEvent = self.positiveScanEvent
            elif mz.ionMode == "-":
                scanEvent = self.negativeScanEvent

            SQLInsert(
                curs,
                "MZs",
                id=mz.id,
                tracer=mz.tid,
                mz=mz.mz,
                lmz=mz.lmz,
                tmz=mz.tmz,
                xcount=mz.xCount,
                scanid=mzxml.getIthMS1Scan(mz.scanIndex, scanEvent).id,
                scanTime=mzxml.getIthMS1Scan(mz.scanIndex, scanEvent).retention_time,
                loading=mz.loading,
                intensity=mz.nIntensity,
                intensityL=mz.lIntensity,
                ionMode=mz.ionMode,
            )

            self.curMZId = self.curMZId + 1

        conn.commit()
        curs.close()
        conn.close()

    # data processing step 2: cluster detected signal pairs with HCA
    def clusterFeaturePairs(self, mzs, reportFunction=None):
        mzbins = {}
        mzbins["+"] = []
        mzbins["-"] = []

        uniquexCounts = list(set([mz.xCount for mz in mzs]))

        # cluster each detected number of carbon atoms separately
        donei = 0
        for xCount in uniquexCounts:
            if reportFunction is not None:
                reportFunction(donei / len(uniquexCounts), "Current Xn: %s" % xCount)
            # cluster each detected number of loadings separately
            for loading in range(self.maxLoading, 0, -1):
                for ionMode in ["+", "-"]:
                    xAct = sorted(
                        [mz for mz in mzs if mz.ionMode == ionMode and mz.xCount == xCount and mz.loading == loading],
                        key=lambda x: x.mz,
                    )

                    doClusterings = []

                    if len(xAct) > 0:
                        lastUnused = 0
                        usedCount = 0

                        # non-consecutive mz regions are prematurely separated without HCA to reduce the number of
                        # signal pairs for HCA (faster processing)
                        for i in range(1, len(xAct)):
                            ppmDiff = (xAct[i].mz - xAct[i - 1].mz) * 1000000 / xAct[i - 1].mz
                            if ppmDiff > self.clustPPM:
                                doClusterings.append(xAct[lastUnused:i])

                                usedCount = usedCount + i - lastUnused
                                lastUnused = i

                        # cluster remainign signal pairs
                        if lastUnused < len(xAct):
                            doClusterings.append(xAct[lastUnused : len(xAct)])
                            usedCount = usedCount + len(xAct) - lastUnused

                        assert usedCount == len(xAct)

                        temp = doClusterings
                        doClusterings = []
                        for t in temp:
                            xAct = sorted([mz for mz in t], key=lambda x: x.scanIndex)

                            lastUnused = 0
                            usedCount = 0

                            for i in range(1, len(xAct)):
                                scanDiff = xAct[i].scanIndex - xAct[i - 1].scanIndex
                                if scanDiff > 20:
                                    doClusterings.append(xAct[lastUnused:i])

                                    usedCount = usedCount + i - lastUnused
                                    lastUnused = i

                            # cluster remainign signal pairs
                            if lastUnused < len(xAct):
                                doClusterings.append(xAct[lastUnused : len(xAct)])
                                usedCount = usedCount + len(xAct) - lastUnused

                            assert usedCount == len(xAct)

                        for doClustering in doClusterings:
                            hc = HierarchicalClustering(
                                doClustering,
                                dist=lambda x, y: x.getValue() - y.getValue(),
                                val=lambda x: x.mz,
                                mean=lambda x, y: x / y,
                                add=lambda x, y: x + y,
                            )

                            for n in cutTreeSized(hc.getTree(), self.clustPPM):
                                mzbins[ionMode].append(n)

        return mzbins

    # store signal pair clusters (2nd data processing step) in the database
    def writeFeaturePairClustersToDB(self, mzbins):
        conn = connect(self.file + getDBSuffix(), isolation_level="DEFERRED")
        conn.execute("""PRAGMA synchronous = OFF""")
        conn.execute("""PRAGMA journal_mode = OFF""")
        curs = conn.cursor()

        for ionMode in ["+", "-"]:
            for mzbin in mzbins[ionMode]:
                SQLInsert(
                    curs,
                    "MZBins",
                    id=self.curMZBinId,
                    mz=mzbin.getValue(),
                    ionMode=ionMode,
                )
                for kid in mzbin.getKids():
                    SQLInsert(
                        curs,
                        "MZBinsKids",
                        mzbinID=self.curMZBinId,
                        mzID=kid.getObject().id,
                    )
                self.curMZBinId = self.curMZBinId + 1

        conn.commit()
        curs.close()
        conn.close()

    def removeImpossibleFeaturePairClusters(self, mzbins):
        mzBinsNew = {}
        for ionMode in mzbins.keys():
            mzBinsNew[ionMode] = [mzbin for mzbin in mzbins[ionMode] if len(mzbin.getKids()) >= self.minSpectraCount]
        return mzBinsNew

    # EXPERIMENTAL: calulcate the mean intensity ratio of two chromatographic peaks at the same retention time
    # in two different EICs. Thus, the internal standardisation is not calculated using the peak area
    # but rather the ratio of each MS peak that contributes to the chromatographic peak.
    def getMeanRatioOfScans(self, eicA, eicB, lib, rib, perfWeighted=True, minInt=1000, minRatiosNecessary=3):
        try:
            if perfWeighted:
                sumEIC = sum(eicA[lib:rib])
                sumEICL = sum(eicB[lib:rib])

                normEIC = eicA
                if sumEICL > sumEIC:
                    normEIC = eicB

                os = [o for o in range(lib, rib) if eicA[o] >= minInt and eicB[o] >= minInt and eicB[o] > 0]
                normSum = sum([normEIC[o] for o in os])
                weights = [1.0 * normEIC[o] / normSum for o in os]
                ratios = [1.0 * eicA[o] / eicB[o] for o in os if eicB[o] > minInt and eicB[o] > 0]

                if len(ratios) >= minRatiosNecessary:
                    assert len(ratios) == len(weights)
                    ratio = sum([ratios[o] * weights[o] for o in range(len(ratios))])
                    deviation = sum([ratios[o] * weights[o] for o in range(len(ratios))])
                    return ratio
                else:
                    return -1
            else:
                ratios = [(eicA[o] / eicB[o]) for o in range(lib, rib) if eicA[o] >= minInt and eicB[o] >= minInt]
                if len(ratios) >= minRatiosNecessary:
                    return mean(ratios)
                else:
                    return -1

        except IndexError:
            return -1

    # data processing step 3: for each signal pair cluster extract the EICs, detect chromatographic peaks
    # present in both EICs at approximately the same retention time and very their chromatographic peak shapes.
    # if all criteria are passed, write this detected feature pair to the database
    def findChromatographicPeaksAndWriteToDB(self, mzbins, mzxml, tracerID, reportFunction=None):
        conn = connect(self.file + getDBSuffix(), isolation_level="DEFERRED")
        conn.execute("""PRAGMA synchronous = OFF""")
        conn.execute("""PRAGMA journal_mode = OFF""")
        curs = conn.cursor()
        chromPeaks = []

        totalBins = len(mzbins["+"]) + len(mzbins["-"])

        # process each signal pair cluster
        for ionMode in ["+", "-"]:
            if ionMode == "+":
                scanEvent = self.positiveScanEvent
            elif ionMode == "-":
                scanEvent = self.negativeScanEvent
            else:
                raise Exception("Unknown Ion Mode")

            if scanEvent != "None":
                for mzbin, i in zip(mzbins[ionMode], range(0, len(mzbins[ionMode]))):
                    if reportFunction is not None:
                        doneBins = i
                        if ionMode == "-":
                            doneBins += len(mzbins["+"])

                        reportFunction(
                            1.0 * doneBins / totalBins,
                            "%d mzbins remaining" % (totalBins - doneBins),
                        )

                    kids = mzbin.getKids()
                    if len(kids) < self.minSpectraCount:
                        continue

                    # calulcate mean mz value for this signal pair cluster
                    meanmz = weightedMean(
                        [kid.getObject().mz for kid in kids],
                        [kid.getObject().nIntensity for kid in kids],
                    )
                    meanmzLabelled = weightedMean(
                        [kid.getObject().lmz for kid in kids],
                        [kid.getObject().lIntensity for kid in kids],
                    )
                    meantmz = weightedMean(
                        [kid.getObject().tmz for kid in kids],
                        [kid.getObject().lIntensity for kid in kids],
                    )

                    xcount = kids[0].getObject().xCount
                    assert all([kid.getObject().xCount == xcount for kid in kids])

                    loading = kids[0].getObject().loading
                    assert all([kid.getObject().loading == loading for kid in kids])

                    # extract the EIC of the native ion and detect its chromatographic peaks
                    # optionally: smoothing
                    eic, times, scanIds, mzs = mzxml.getEIC(meanmz, self.chromPeakPPM, filterLine=scanEvent)
                    eicBaseline = self.BL.getBaseline(copy(eic), times)
                    eicSmoothed = smoothDataSeries(
                        times,
                        copy(eic),
                        windowLen=self.eicSmoothingWindowSize,
                        window=self.eicSmoothingWindow,
                        polynom=self.eicSmoothingPolynom,
                    )
                    # extract the EIC of the labelled ion and detect its chromatographic peaks
                    # optionally: smoothing
                    eicL, times, scanIds, mzsL = mzxml.getEIC(meanmzLabelled, self.chromPeakPPM, filterLine=scanEvent)
                    eicLBaseline = self.BL.getBaseline(copy(eicL), times)
                    eicLSmoothed = smoothDataSeries(
                        times,
                        copy(eicL),
                        windowLen=self.eicSmoothingWindowSize,
                        window=self.eicSmoothingWindow,
                        polynom=self.eicSmoothingPolynom,
                    )

                    # determine boundaries for chromatographic peak picking
                    minInd = min([kid.getObject().scanIndex for kid in kids])
                    maxInd = max([kid.getObject().scanIndex for kid in kids])

                    ## TODO needs to be optimized
                    # startIndex=max(0, minInd-int(ceil(self.scales[1]*10)))
                    # endIndex=min(len(eic)-1, int(ceil(maxInd+self.scales[1]*10)))
                    startIndex = 0
                    endIndex = len(eic) - 1

                    peaksN = []
                    try:
                        peaksN = self.CP.getPeaksFor(times, eicSmoothed, startIndex=startIndex, endIndex=endIndex)
                    except Exception as ex:
                        self.printMessage("Errora: %s" % str(ex), type="error")

                    peaksL = []
                    try:
                        peaksL = self.CP.getPeaksFor(
                            times,
                            eicLSmoothed,
                            startIndex=startIndex,
                            endIndex=endIndex,
                        )
                    except Exception as ex:
                        self.printMessage("Errorb: %s" % str(ex), type="error")

                    # get EICs of M+1, M'-1 and M'+1 for the database
                    eicfirstiso, timesL, scanIdsL, mzsfirstiso = mzxml.getEIC(
                        meanmz + 1.00335484 / loading,
                        self.chromPeakPPM,
                        filterLine=scanEvent,
                    )
                    # eicfirstisoSmoothed = smoothDataSeries(times, eicfirstiso, windowLen=self.eicSmoothingWindowSize, window=self.eicSmoothingWindow, polynom=self.eicSmoothingPolynom)

                    eicLfirstiso, timesL, scanIdsL, mzsLfirstiso = mzxml.getEIC(
                        meanmzLabelled - 1.00335484 / loading,
                        self.chromPeakPPM,
                        filterLine=scanEvent,
                    )
                    # eicLfirstisoSmoothed = smoothDataSeries(times, eicLfirstiso, windowLen=self.eicSmoothingWindowSize, window=self.eicSmoothingWindow, polynom=self.eicSmoothingPolynom)

                    eicLfirstisoconjugate, timesL, scanIdsL, mzsLfirstisoconjugate = mzxml.getEIC(
                        meanmzLabelled + 1.00335484 / loading,
                        self.chromPeakPPM,
                        filterLine=scanEvent,
                    )
                    # eicLfirstisoconjugateSmoothed = smoothDataSeries(times, eicLfirstisoconjugate, windowLen=self.eicSmoothingWindowSize, window=self.eicSmoothingWindow, polynom=self.eicSmoothingPolynom)

                    peaksBoth = []

                    # match detected chromatographic peaks from both EICs
                    for peakN in peaksN:
                        closestMatch = -1
                        closestOffset = 10000000000

                        # search closest peak pair
                        for li, peakL in enumerate(peaksL):
                            ##print(peakN.peakIndex, peakL.peakIndex, meanmz, meanmzLabelled)
                            ##if abs(peakN[0] - peakL[0]) <= self.peakCenterError:
                            if abs(peakN.peakIndex - peakL.peakIndex) <= self.peakCenterError:
                                if closestMatch == -1 or closestOffset > abs(peakN.peakIndex - peakL.peakIndex):
                                    closestMatch = li
                                    closestOffset = abs(peakN.peakIndex - peakL.peakIndex)

                        if closestMatch != -1:
                            peakL = peaksL[closestMatch]

                            # Calculate new peak stats
                            fwhm_M, area_M, snr_M = getPeakStats(np.vstack((times, eicSmoothed)), int(peakN.peakIndex - peakN.peakLeftFlank), peakN.peakIndex, int(peakN.peakIndex + peakN.peakRightFlank))
                            fwhm_Mp, area_Mp, snr_Mp = getPeakStats(np.vstack((times, eicLSmoothed)), int(peakL.peakIndex - peakL.peakLeftFlank), peakL.peakIndex, int(peakL.peakIndex + peakL.peakRightFlank))

                            # print("M: ", int(peakN.peakIndex - peakN.peakLeftFlank), peakN.peakIndex, int(peakN.peakIndex + peakN.peakRightFlank))
                            # print(f"   FWHM: {fwhm_M}, area: {area_M}, SNR: {snr_M}")
                            # print("Mp: ", int(peakL.peakIndex - peakL.peakLeftFlank), peakL.peakIndex, int(peakL.peakIndex + peakL.peakRightFlank))
                            # print(f"   FWHM: {fwhm_Mp}, area: {area_Mp}, SNR: {snr_Mp}")

                            peak = ChromPeakPair(
                                mz=meanmz,
                                lmz=meanmzLabelled,
                                tmz=meantmz,
                                xCount=xcount,
                                loading=loading,
                                ionMode=ionMode,
                                NPeakCenter=peakN.peakIndex,
                                NPeakCenterMin=times[peakN.peakIndex],
                                NPeakScale=peakN.peakScale,
                                NSNR=peakN.peakSNR,
                                NPeakArea=peakN.peakArea,
                                LPeakCenter=peakL.peakIndex,
                                LPeakCenterMin=times[peakL.peakIndex],
                                LPeakScale=peakL.peakScale,
                                LSNR=peakL.peakSNR,
                                LPeakArea=peakL.peakArea,
                                NXIC=eic,
                                LXIC=eicL,
                                times=times,
                                fDesc=[],
                                adducts=[],
                                heteroAtomsFeaturePairs=[],
                                NXICSmoothed=eicSmoothed,
                                LXICSmoothed=eicLSmoothed,
                                NBorderLeft=peakN.peakLeftFlank,
                                NBorderRight=peakN.peakRightFlank,
                                LBorderLeft=peakL.peakLeftFlank,
                                LBorderRight=peakL.peakRightFlank,
                                isotopeRatios=[],
                                mzDiffErrors=Bunch(),
                                comments=[],
                                artificialEICLShift=0,
                                new_FWHM_M=fwhm_M,
                                new_Area_M=area_M,
                                new_SNR_M=snr_M,
                                new_FWHM_Mp=fwhm_Mp,
                                new_Area_Mp=area_Mp,
                                new_SNR_Mp=snr_Mp,
                            )

                            lb = int(
                                max(
                                    0,
                                    min(
                                        peak.NPeakCenter - peak.NBorderLeft,
                                        peak.LPeakCenter - peak.LBorderLeft,
                                    ),
                                )
                            )
                            rb = int(
                                min(
                                    max(
                                        peak.NPeakCenter + peak.NBorderRight,
                                        peak.LPeakCenter + peak.LBorderRight,
                                    )
                                    + 1,
                                    len(eic) - 1,
                                )
                            )
                            peakN = eicSmoothed[lb:rb]
                            peakL = eicLSmoothed[lb:rb]

                            def findBestArtificialShift(eicN, eicL, lb, rb, shiftFrom=0, shiftTo=0):
                                correlations = []

                                for artShift in range(shiftFrom, shiftTo + 1):
                                    peakN = eicN[lb:rb]
                                    peakL = eicL[(lb + artShift) : (rb + artShift)]
                                    silRatios = [
                                        peakN[i] / peakL[i]
                                        for i in range(
                                            int(len(peakN) * 0.25),
                                            int(len(peakN) * 0.75) + 1,
                                        )
                                        if peakL[i] > 0 and peakN[i] > 0
                                    ]
                                    correlations.append(
                                        Bunch(
                                            correlation=corr(peakN, peakL),
                                            artificialShift=artShift,
                                            silRatios=silRatios,
                                            peakNInts=[
                                                peakN[i]
                                                for i in range(
                                                    int(len(peakN) * 0.25),
                                                    int(len(peakN) * 0.75) + 1,
                                                )
                                                if peakL[i] > 0 and peakN[i] > 0
                                            ],
                                            peakLInts=[
                                                peakL[i]
                                                for i in range(
                                                    int(len(peakN) * 0.25),
                                                    int(len(peakN) * 0.75) + 1,
                                                )
                                                if peakL[i] > 0 and peakN[i] > 0
                                            ],
                                        )
                                    )
                                bestFit = max(correlations, key=lambda x: x.correlation)

                                return bestFit

                            co = findBestArtificialShift(
                                eicSmoothed,
                                eicLSmoothed,
                                lb,
                                rb,
                                self.artificialMPshift_start,
                                self.artificialMPshift_stop,
                            )

                            # check peak shape (Pearson correlation)

                            if co.correlation >= self.minPeakCorr and ((not self.checkPeaksRatio) or self.minPeaksRatio <= (peak.NPeakArea / peak.LPeakArea) <= self.maxPeaksRatio):
                                peak.peaksCorr = co.correlation
                                peak.silRatios = Bunch(
                                    silRatios=co.silRatios,
                                    peakNInts=co.peakNInts,
                                    peakLInts=co.peakLInts,
                                )
                                if co.artificialShift != 0:
                                    peak.artificialEICLShift = co.artificialShift

                                peaksBoth.append(peak)

                                libs = int(
                                    max(
                                        peak.NPeakCenter - floor(peak.NBorderLeft * 0.33),
                                        peak.LPeakCenter - floor(peak.LBorderLeft * 0.33),
                                    )
                                )
                                ribs = (
                                    int(
                                        min(
                                            peak.NPeakCenter + floor(peak.NBorderRight * 0.33),
                                            peak.LPeakCenter + floor(peak.LBorderRight * 0.33),
                                        )
                                    )
                                    + 1
                                )

                                peak.peaksRatio = self.getMeanRatioOfScans(
                                    eic,
                                    eicL,
                                    libs,
                                    ribs,
                                    minInt=self.intensityThreshold,
                                )
                                peak.peaksRatioMp1 = self.getMeanRatioOfScans(
                                    eicfirstiso,
                                    eic,
                                    libs,
                                    ribs,
                                    minInt=self.intensityThreshold,
                                )
                                peak.peaksRatioMPm1 = self.getMeanRatioOfScans(
                                    eicLfirstiso,
                                    eicL,
                                    libs,
                                    ribs,
                                    minInt=self.intensityThreshold,
                                )

                    # check, if enough MS scans fulfill the requirements (isotope patterns)
                    for kid in kids:
                        kido = kid.getObject()

                        closestPeak = -1
                        distance = len(eic)
                        i = 0
                        for peak in peaksBoth:
                            if abs(kido.scanIndex - peak.NPeakCenter) < distance:
                                closestPeak = i
                                distance = abs(kido.scanIndex - peak.NPeakCenter)
                            i = i + 1
                        if closestPeak != -1 and distance < (
                            min(
                                peaksBoth[closestPeak].NBorderLeft,
                                peaksBoth[closestPeak].NBorderRight,
                            )
                            * 2
                        ):
                            peaksBoth[closestPeak].assignedMZs.append(kido)

                    curChromPeaks = []
                    for peak in peaksBoth:
                        if len(peak.assignedMZs) >= self.minSpectraCount:
                            curChromPeaks.append(peak)

                    for peak in curChromPeaks:
                        lb = int(
                            max(
                                0,
                                min(
                                    peak.NPeakCenter - peakAbundanceUseSignalsSides,
                                    peak.LPeakCenter - peakAbundanceUseSignalsSides,
                                ),
                            )
                        )
                        rb = int(
                            min(
                                max(
                                    peak.NPeakCenter + peakAbundanceUseSignalsSides,
                                    peak.LPeakCenter + peakAbundanceUseSignalsSides,
                                )
                                + 1,
                                len(eic) - 1,
                            )
                        )
                        peakN = eic[lb:rb]
                        peakL = eicL[lb:rb]

                        peak.NPeakAbundance = max(peakN)  # mean(peakN)
                        peak.LPeakAbundance = max(peakL)  # mean(peakL)

                    # write detected feature pair to the database
                    if len(curChromPeaks) > 0:
                        assert len(eic) == len(eicL) == len(eicfirstiso) == len(eicLfirstiso) == len(eicLfirstisoconjugate) == len(times)

                        keepScans = set()
                        for cs in range(len(eic)):
                            if eic[cs] > 0 or eicL[cs] > 0 or eicfirstiso[cs] > 0 or eicLfirstiso[cs] > 0 or eicLfirstisoconjugate[cs] > 0:
                                keepScans.add(cs)
                        keepScans.add(0)
                        keepScans.add(len(eic) - 1)

                        # save the native and labelled EICs to the database
                        SQLInsert(
                            curs,
                            "XICs",
                            id=self.curEICId,
                            avgmz=meanmz,
                            xcount=xcount,
                            loading=kids[0].getObject().loading,
                            polarity=ionMode,
                            xic=";".join(["%f" % i for i in eic]),
                            xicL=";".join(["%f" % i for i in eicL]),
                            xicfirstiso=";".join(["%f" % i for i in eicfirstiso]),
                            xicLfirstiso=";".join(["%f" % i for i in eicLfirstiso]),
                            xicLfirstisoconjugate=";".join(["%f" % i for i in eicLfirstisoconjugate]),
                            xic_smoothed=";".join(["%f" % i for i in eicSmoothed]),
                            xicL_smoothed=";".join(["%f" % i for i in eicLSmoothed]),
                            xic_baseline=";".join(["%f" % i for i in eicBaseline]),
                            xicL_baseline=";".join(["%f" % i for i in eicLBaseline]),
                            mzs=";".join(["%f" % (mzs[i]) for i in range(0, len(eic))]),
                            mzsL=";".join(["%f" % (mzsL[i]) for i in range(0, len(eic))]),
                            mzsfirstiso=";".join(["%f" % (mzsfirstiso[i]) for i in range(0, len(eic))]),
                            mzsLfirstiso=";".join(["%f" % (mzsLfirstiso[i]) for i in range(0, len(eic))]),
                            mzsLfirstisoconjugate=";".join(["%f" % (mzsLfirstisoconjugate[i]) for i in range(0, len(eic))]),
                            times=";".join(["%f" % (times[i]) for i in range(0, len(times))]),
                            scanCount=len(eic),
                            allPeaks=base64.b64encode(dumps({"peaksN": peaksN, "peaksL": peaksL})).decode("utf-8"),
                        )

                        # save the detected feature pairs in these EICs to the database
                        for peak in curChromPeaks:
                            peak.id = self.curPeakId
                            peak.eicID = self.curEICId
                            adjcCount = peak.xCount
                            peak.correctedXCount = peak.xCount

                            if self.performCorrectCCount and False:
                                adjcCount = adjcCount + getAtomAdd(self.purityL, peak.xCount) + getAtomAdd(self.purityN, peak.xCount)
                                peak.correctedXCount = adjcCount

                            SQLInsert(
                                curs,
                                "chromPeaks",
                                id=peak.id,
                                tracer=tracerID,
                                eicID=peak.eicID,
                                mz=peak.mz,
                                lmz=peak.lmz,
                                tmz=peak.tmz,
                                xcount=peak.correctedXCount,
                                xcountId=peak.xCount,
                                Loading=peak.loading,
                                ionMode=ionMode,
                                NPeakCenter=peak.NPeakCenter,
                                NPeakCenterMin=peak.NPeakCenterMin,
                                NPeakScale=peak.NPeakScale,
                                NSNR=peak.NSNR,
                                NPeakArea=peak.NPeakArea,
                                NPeakAbundance=peak.NPeakAbundance,
                                NBorderLeft=peak.NBorderLeft,
                                NBorderRight=peak.NBorderRight,
                                LPeakCenter=peak.LPeakCenter,
                                LPeakCenterMin=peak.LPeakCenterMin,
                                LPeakScale=peak.LPeakScale,
                                LSNR=peak.LSNR,
                                LPeakArea=peak.LPeakArea,
                                LPeakAbundance=peak.LPeakAbundance,
                                LBorderLeft=peak.LBorderLeft,
                                LBorderRight=peak.LBorderRight,
                                peaksCorr=peak.peaksCorr,
                                heteroAtoms="",
                                adducts="",
                                heteroAtomsFeaturePairs="",
                                massSpectrumID=0,
                                assignedMZs=base64.b64encode(dumps(peak.assignedMZs)).decode("utf-8"),
                                fDesc=base64.b64encode(dumps([])).decode("utf-8"),
                                peaksRatio=peak.peaksRatio,
                                peaksRatioMp1=peak.peaksRatioMp1,
                                peaksRatioMPm1=peak.peaksRatioMPm1,
                                peakType="patternFound",
                                comments=base64.b64encode(dumps(peak.comments)).decode("utf-8"),
                                artificialEICLShift=peak.artificialEICLShift,
                                new_FWHM_M=peak.new_FWHM_M,
                                new_FWHM_Mp=peak.new_FWHM_Mp,
                                new_Area_M=peak.new_Area_M,
                                new_Area_Mp=peak.new_Area_Mp,
                                new_SNR_M=peak.new_SNR_M,
                                new_SNR_Mp=peak.new_SNR_Mp,
                            )

                            SQLInsert(
                                curs,
                                "allChromPeaks",
                                id=peak.id,
                                tracer=tracerID,
                                eicID=peak.eicID,
                                mz=peak.mz,
                                lmz=peak.lmz,
                                tmz=peak.tmz,
                                xcount=peak.correctedXCount,
                                xcountId=peak.xCount,
                                Loading=peak.loading,
                                ionMode=ionMode,
                                NPeakCenter=peak.NPeakCenter,
                                NPeakCenterMin=peak.NPeakCenterMin,
                                NPeakScale=peak.NPeakScale,
                                NSNR=peak.NSNR,
                                NPeakArea=peak.NPeakArea,
                                NPeakAbundance=peak.NPeakAbundance,
                                NBorderLeft=peak.NBorderLeft,
                                NBorderRight=peak.NBorderRight,
                                LPeakCenter=peak.LPeakCenter,
                                LPeakCenterMin=peak.LPeakCenterMin,
                                LPeakScale=peak.LPeakScale,
                                LSNR=peak.LSNR,
                                LPeakArea=peak.LPeakArea,
                                LPeakAbundance=peak.LPeakAbundance,
                                LBorderLeft=peak.LBorderLeft,
                                LBorderRight=peak.LBorderRight,
                                peaksCorr=peak.peaksCorr,
                                heteroAtoms="",
                                adducts="",
                                heteroAtomsFeaturePairs="",
                                assignedMZs=len(peak.assignedMZs),
                                fDesc=base64.b64encode(dumps([])).decode("utf-8"),
                                peaksRatio=peak.peaksRatio,
                                peaksRatioMp1=peak.peaksRatioMp1,
                                peaksRatioMPm1=peak.peaksRatioMPm1,
                                peakType="patternFound",
                                comments=base64.b64encode(dumps(peak.comments)).decode("utf-8"),
                                artificialEICLShift=peak.artificialEICLShift,
                                new_FWHM_M=peak.new_FWHM_M,
                                new_FWHM_Mp=peak.new_FWHM_Mp,
                                new_Area_M=peak.new_Area_M,
                                new_Area_Mp=peak.new_Area_Mp,
                                new_SNR_M=peak.new_SNR_M,
                                new_SNR_Mp=peak.new_SNR_Mp,
                            )

                            chromPeaks.append(peak)
                            self.curPeakId = self.curPeakId + 1

                    self.curEICId = self.curEICId + 1
        conn.commit()
        curs.close()
        conn.close()
        return chromPeaks

    # data processing step 4: remove those feature pairs, which represent incorrect pairings of
    # isotoplogs of either the native or the labelled or both ions. Such incorrect pairings always have
    # and increased mz value and/or a decreased number of labelled carbon atoms. Such identified
    # incorrect pairings are then removed from the database and thus the processing results
    def removeFalsePositiveFeaturePairsAndUpdateDB(self, chromPeaks, reportFunction=None):
        conn = connect(self.file + getDBSuffix(), isolation_level="DEFERRED")
        conn.execute("""PRAGMA synchronous = OFF""")
        conn.execute("""PRAGMA journal_mode = OFF""")
        curs = conn.cursor()

        todel = {}

        # iterate over all detected feature pairs and compare those
        # if the a) originate from the same ionisation mode
        #    and b) if they are not the same feature pair
        for a in range(0, len(chromPeaks)):
            if reportFunction is not None:
                reportFunction(
                    1.0 * a / len(chromPeaks),
                    "%d feature pairs remaining" % (len(chromPeaks) - a),
                )

            peakA = chromPeaks[a]
            for b in range(0, len(chromPeaks)):
                peakB = chromPeaks[b]
                if a != b and peakA.ionMode == peakB.ionMode and (peakA.mz < peakB.mz or abs(peakA.mz - peakB.mz) <= (peakA.mz * 2.5 * self.ppm / 1000000.0)):
                    # not same feature pair but same ionisation mode

                    if abs(peakA.NPeakCenter - peakB.NPeakCenter) <= self.peakCenterError:
                        # same retention time
                        if abs(peakA.mz - peakB.mz) <= (peakA.mz * 2.5 * self.ppm / 1000000.0):
                            # same mz value
                            if isinstance(peakA.xCount, float) and (peakB.xCount - peakA.xCount) == 2 and peakB.loading == peakA.loading and ((peakB.LPeakArea / peakA.LPeakArea) <= 0.1):
                                # incorrect matching of 18O atoms detected
                                # e.g. a and b have 869.4153; a has C39 and b has C41; peaksRatio(a)=1.46, peaksRatio(b)=43.48
                                # --> 1.46/43.48~0.0335 (which is in good agreement with On) and thus b has to be removed
                                if b not in todel.keys():
                                    todel[b] = []
                                todel[b].append("incorrectly matched hetero atoms (most likely O) with " + str(peakA.mz) + " " + str(peakA.xCount))
                            elif isinstance(peakA.xCount, float) and (peakB.xCount - peakA.xCount) in [1, 2, 3] and peakB.loading == peakA.loading:
                                # different number of 1 atoms (peakA has less than peakB)
                                if a not in todel.keys():
                                    todel[a] = []
                                todel[a].append("same mz but reduced number of carbon atoms with " + str(peakB.mz) + " " + str(peakB.xCount))

                            if isinstance(peakA.xCount, float) and (peakA.xCount * 2) == peakB.xCount and peakA.loading == peakB.loading:
                                # a is intermediate pairing of polymer
                                if a not in todel.keys():
                                    todel[a] = []
                                todel[a].append("polymer mismatch with half the number of carbon atoms with " + str(peakB.mz) + " " + str(peakB.xCount))

                            if isinstance(peakA.xCount, float) and (peakA.xCount * 2) == peakB.xCount and peakA.loading == 1 and peakB.loading == 2:
                                # a is intermediate pairing of polymer
                                if a not in todel.keys():
                                    todel[a] = []
                                todel[a].append("singly charged mismatch with half the number of carbon atoms with " + str(peakB.mz) + " " + str(peakB.xCount))

                            if isinstance(peakA.xCount, float) and 0 <= (peakB.xCount - peakA.xCount * 3) <= 2 and peakA.loading == 1 and peakB.loading == 3:
                                # a is intermediate pairing of polymer
                                if a not in todel.keys():
                                    todel[a] = []
                                todel[a].append("singly charged mismatch with third the number of carbon atoms with " + str(peakB.mz) + " " + str(peakB.xCount))

                        elif abs(peakB.mz - peakA.mz - 1.00335 / peakA.loading) <= (peakA.mz * 2.5 * self.ppm / 1000000.0) and peakB.xCount == peakA.xCount and peakB.NPeakArea < 2 * peakA.NPeakArea:
                            ## increased m/z value by one carbon atom, but the same number of labeling atoms ## happens quite often with 15N-labeled metabolites e.g. 415.21374 and 415.7154 or 414.6978 and 415.19966
                            if b not in todel.keys():
                                todel[b] = []
                            todel[b].append("increased m/z by one labeling atoms while the number of labeled atoms is identical." + str(peakA.mz))

                        else:
                            for i in [1, 2, 3]:
                                if peakA.loading == peakB.loading and abs(peakB.mz - peakA.mz - i * self.xOffset / peakA.loading) <= peakA.mz * 2.5 * self.ppm / 1000000.0 and isinstance(peakA.xCount, float) and (peakA.xCount - peakB.xCount) in [1, 2, 3]:
                                    # b has an mz offset and a reduced number of carbon atoms (by one labelling atom)
                                    if b not in todel.keys():
                                        todel[b] = []
                                    todel[b].append("increased mz (by %d) and reduced number of labeling atoms (by %d) with " % (i, peakA.xCount - peakB.xCount) + str(peakA.mz) + " " + str(peakA.xCount))
                                elif peakA.loading == peakB.loading and abs(peakB.mz - peakA.mz - i * 1.00335 / peakA.loading) <= peakA.mz * 2.5 * self.ppm / 1000000.0 and peakA.xCount == peakB.xCount and peakB.NPeakArea <= peakA.NPeakArea * 0.1:
                                    # b has an mz offset and a reduced number of carbon atoms (by one labelling atom)
                                    if b not in todel.keys():
                                        todel[b] = []
                                    todel[b].append("increased mz (by %d 13C) and same number of labeling atoms with " % (i) + str(peakA.mz) + " " + str(peakA.xCount))

        for a in range(0, len(chromPeaks)):
            for b in range(0, len(chromPeaks)):
                peakA = chromPeaks[a]
                peakB = chromPeaks[b]
                if a != b and peakA.ionMode == peakB.ionMode:
                    if abs(peakA.NPeakCenter - peakB.NPeakCenter) <= self.peakCenterError:
                        # same chrom peak
                        if abs(peakA.mz - peakB.mz - (peakA.lmz - peakA.mz)) <= peakA.mz * 2 * self.ppm / 1000000.0 and ((isinstance(peakA.xCount, float) and abs(peakA.xCount * 2 - peakB.xCount) <= 1) or (peakA.xCount == peakB.xCount)):
                            if a not in todel.keys():
                                todel[a] = []
                            todel[a].append("dimer (1)" + str(peakB.mz) + " " + str(peakB.xCount))
        delet = []
        for x in todel.keys():
            delet.append(x)
        delet.sort()
        delet.reverse()
        for dele in delet:
            cp = chromPeaks.pop(dele)
            curs.execute("DELETE FROM chromPeaks WHERE id=%d" % cp.id)
            curs.execute("UPDATE allChromPeaks SET comment='%s' WHERE id=%d" % (",".join(todel[dele]), cp.id))

        curs.execute("DELETE FROM XICs WHERE id NOT IN (SELECT eicID FROM chromPeaks)")

        conn.commit()
        curs.close()
        conn.close()

    # data processing step 5: in full metabolome labeling experiments hetero atoms (e.g. S, Cl) may
    # show characterisitc isotope patterns on the labeled metabolite ion side. There, these peaks may not be
    # dominated by the usually much more abundant carbon isotopes and can thus be easier seen. However, for
    # low abundant metabolite ions these isotope peaks may not be present at all. The search is performed
    # on a MS scan level and does not directly use the chromatographic information (no chromatographic
    # peak picking is performed)
    def annotateFeaturePairs(self, chromPeaks, mzxml, tracer, reportFunction=None):
        conn = connect(self.file + getDBSuffix(), isolation_level="DEFERRED")
        conn.execute("""PRAGMA synchronous = OFF""")
        conn.execute("""PRAGMA journal_mode = OFF""")
        curs = conn.cursor()

        self.postMessageToProgressWrapper("text", "%s: Annotating feature pairs" % tracer.name)
        for i in range(0, len(chromPeaks)):
            if reportFunction is not None:
                reportFunction(
                    1.0 * i / len(chromPeaks),
                    "%d features remaining" % (len(chromPeaks) - i),
                )

            peak = chromPeaks[i]

            ## Annotate hetero atoms
            for pIso in self.heteroAtoms:
                pIsotope = self.heteroAtoms[pIso]

                mz = 0
                if pIsotope.mzOffset < 0:
                    mz = peak.mz + pIsotope.mzOffset / peak.loading  # delta m/z is negative, therefore this decreases the search m/z
                else:
                    mz = peak.lmz + pIsotope.mzOffset / peak.loading  # delta m/z is positive

                refMz = 0
                if pIsotope.mzOffset < 0:
                    refMz = peak.mz
                else:
                    refMz = peak.lmz

                scanEvent = ""
                if peak.ionMode == "+":
                    scanEvent = self.positiveScanEvent
                elif peak.ionMode == "-":
                    scanEvent = self.negativeScanEvent

                for haCount in range(pIsotope.minCount, pIsotope.maxCount + 1):
                    if haCount == 0:
                        continue

                    for curScanNum in range(
                        int(
                            max(
                                peak.NPeakCenter - peak.NBorderLeft,
                                peak.LPeakCenter - peak.LBorderLeft,
                            )
                        ),
                        int(
                            min(
                                peak.NPeakCenter + peak.NBorderRight,
                                peak.LPeakCenter + peak.LBorderRight,
                            )
                        )
                        + 1,
                    ):
                        scan = mzxml.getIthMS1Scan(curScanNum, scanEvent)
                        if scan is not None:
                            mzBounds = scan.findMZ(refMz, self.ppm)
                            mzBounds = scan.getMostIntensePeak(mzBounds[0], mzBounds[1])
                            if mzBounds != -1:
                                peakIntensity = scan.intensity_list[mzBounds]

                                isoBounds = scan.findMZ(mz, self.ppm)
                                isoBounds = scan.getMostIntensePeak(isoBounds[0], isoBounds[1])
                                if isoBounds != -1:
                                    isoIntensity = scan.intensity_list[isoBounds]

                                    relativeIntensity = isoIntensity / peakIntensity
                                    if abs(relativeIntensity - pIsotope.relativeAbundance * haCount) < self.hAIntensityError:
                                        if not (pIso in peak.heteroIsotopologues):
                                            peak.heteroIsotopologues[pIso] = {}
                                        if not (haCount in peak.heteroIsotopologues[pIso]):
                                            peak.heteroIsotopologues[pIso][haCount] = []
                                        peak.heteroIsotopologues[pIso][haCount].append(
                                            (
                                                scan.id,
                                                relativeIntensity,
                                                pIsotope.relativeAbundance * haCount,
                                            )
                                        )

            rmHIso = []
            for hI in peak.heteroIsotopologues:
                for haCount in peak.heteroIsotopologues[hI]:
                    if len(peak.heteroIsotopologues[hI][haCount]) < self.hAMinScans:
                        rmHIso.append((hI, haCount))

            for rmHI, rmHACount in rmHIso:
                del peak.heteroIsotopologues[rmHI][rmHACount]
            rmHIso = []
            for hi in peak.heteroIsotopologues:
                if len(peak.heteroIsotopologues[hi]) == 0:
                    rmHIso.append(hi)
            for rmHI in rmHIso:
                del peak.heteroIsotopologues[rmHI]

            self.getMostLikelyHeteroIsotope(peak.heteroIsotopologues)

            ## Annotate isotopolog ratios

            def findRatiosForMZs(mzFrom, mzTo, fromScan, toScan, mzxml, scanEvent, ppm):
                isoRatios = []
                for curScanNum in range(fromScan, toScan):
                    scan = mzxml.getIthMS1Scan(curScanNum, scanEvent)
                    if scan is not None:
                        mzBounds = scan.findMZ(mzTo, ppm)
                        mzBounds = scan.getMostIntensePeak(mzBounds[0], mzBounds[1])
                        if mzBounds != -1:
                            peakIntensity = scan.intensity_list[mzBounds]

                            isoBounds = scan.findMZ(mzFrom, ppm)
                            isoBounds = scan.getMostIntensePeak(isoBounds[0], isoBounds[1])

                            if isoBounds != -1:
                                isoPeakIntensity = scan.intensity_list[isoBounds]

                                isoRatios.append(
                                    Bunch(
                                        ratio=isoPeakIntensity / peakIntensity,
                                        refInt=peakIntensity,
                                    )
                                )

                return isoRatios

            peak.isotopeRatios = []
            for i in range(1, self.calcIsoRatioNative + 1):
                isoRatios = findRatiosForMZs(
                    peak.mz + 1.00335484 * i / peak.loading,
                    peak.mz,
                    int(
                        max(
                            peak.NPeakCenter - peak.NBorderLeft,
                            peak.LPeakCenter - peak.LBorderLeft,
                        )
                    ),
                    int(
                        min(
                            peak.NPeakCenter + peak.NBorderRight,
                            peak.LPeakCenter + peak.LBorderRight,
                        )
                    )
                    + 1,
                    mzxml,
                    scanEvent,
                    self.ppm,
                )
                observedRatioMean = weightedMean([t.ratio for t in isoRatios], [t.refInt for t in isoRatios])
                observedRatioSD = weightedSd([t.ratio for t in isoRatios], [t.refInt for t in isoRatios])
                if isinstance(peak.xCount, str):
                    theoreticalRatio = -1
                else:
                    theoreticalRatio = getNormRatio(self.purityN, peak.xCount, i)
                peak.isotopeRatios.append(
                    Bunch(
                        type="native",
                        subs=i,
                        observedRatioMean=observedRatioMean,
                        observedRatioSD=observedRatioSD,
                        theoreticalRatio=theoreticalRatio,
                    )
                )

            for i in range(-1, self.calcIsoRatioLabelled - 1, -1):
                isoRatios = findRatiosForMZs(
                    peak.lmz + 1.00335484 * i / peak.loading,
                    peak.lmz,
                    int(
                        max(
                            peak.NPeakCenter - peak.NBorderLeft,
                            peak.LPeakCenter - peak.LBorderLeft,
                        )
                    ),
                    int(
                        min(
                            peak.NPeakCenter + peak.NBorderRight,
                            peak.LPeakCenter + peak.LBorderRight,
                        )
                    )
                    + 1,
                    mzxml,
                    scanEvent,
                    self.ppm,
                )
                observedRatioMean = weightedMean([t.ratio for t in isoRatios], [t.refInt for t in isoRatios])
                observedRatioSD = weightedSd([t.ratio for t in isoRatios], [t.refInt for t in isoRatios])
                if isinstance(peak.xCount, str):
                    theoreticalRatio = -1
                else:
                    theoreticalRatio = getNormRatio(self.purityL, peak.xCount, abs(i))
                peak.isotopeRatios.append(
                    Bunch(
                        type="labelled",
                        subs=abs(i),
                        observedRatioMean=observedRatioMean,
                        observedRatioSD=observedRatioSD,
                        theoreticalRatio=theoreticalRatio,
                    )
                )
            if self.metabolisationExperiment:
                for i in range(1, self.calcIsoRatioMoiety + 1):
                    isoRatios = findRatiosForMZs(
                        peak.lmz + i * 1.00335484 / peak.loading,
                        peak.lmz,
                        int(
                            max(
                                peak.NPeakCenter - peak.NBorderLeft,
                                peak.LPeakCenter - peak.LBorderLeft,
                            )
                        ),
                        int(
                            min(
                                peak.NPeakCenter + peak.NBorderRight,
                                peak.LPeakCenter + peak.LBorderRight,
                            )
                        )
                        + 1,
                        mzxml,
                        scanEvent,
                        self.ppm,
                    )
                    observedRatioMean = weightedMean([t.ratio for t in isoRatios], [t.refInt for t in isoRatios])
                    observedRatioSD = weightedSd([t.ratio for t in isoRatios], [t.refInt for t in isoRatios])
                    peak.isotopeRatios.append(
                        Bunch(
                            type="moiety",
                            subs=i,
                            observedRatioMean=observedRatioMean,
                            observedRatioSD=observedRatioSD,
                            theoreticalRatio=0,
                        )
                    )

            def findMZDifferenceRelativeToXnForMZs(mzFrom, mzTo, tmz, fromScan, toScan, mzxml, scanEvent, ppm):
                diffs = []
                for curScanNum in range(fromScan, toScan):
                    scan = mzxml.getIthMS1Scan(curScanNum, scanEvent)
                    if scan is not None:
                        mzBounds = scan.findMZ(mzTo, ppm)
                        mzBounds = scan.getMostIntensePeak(mzBounds[0], mzBounds[1])
                        if mzBounds != -1:
                            peakMZ = scan.mz_list[mzBounds]

                            refBounds = scan.findMZ(mzFrom, ppm)
                            refBounds = scan.getMostIntensePeak(refBounds[0], refBounds[1])

                            if refBounds != -1:
                                isoPeakMZ = scan.mz_list[refBounds]

                                diffs.append((abs(isoPeakMZ - peakMZ) - tmz) * 1000000.0 / mzFrom)

                return diffs

            diffs = findMZDifferenceRelativeToXnForMZs(
                peak.mz,
                peak.lmz,
                peak.lmz - peak.mz,
                int(
                    max(
                        peak.NPeakCenter - peak.NBorderLeft,
                        peak.LPeakCenter - peak.LBorderLeft,
                    )
                ),
                int(
                    min(
                        peak.NPeakCenter + peak.NBorderRight,
                        peak.LPeakCenter + peak.LBorderRight,
                    )
                )
                + 1,
                mzxml,
                scanEvent,
                self.ppm,
            )
            peak.mzDiffErrors = Bunch(mean=mean(diffs), sd=sd(diffs), vals=diffs)

        ## TODO this is really, really slow as
        for i in range(len(chromPeaks)):
            peak = chromPeaks[i]
            curs.execute(
                "UPDATE chromPeaks SET heteroAtoms=? WHERE id=?",
                (
                    base64.b64encode(dumps(peak.heteroIsotopologues)).decode("utf-8"),
                    peak.id,
                ),
            )
            curs.execute(
                "UPDATE allChromPeaks SET heteroAtoms=? WHERE id=?",
                (
                    base64.b64encode(dumps(peak.heteroIsotopologues)).decode("utf-8"),
                    peak.id,
                ),
            )

            curs.execute(
                "UPDATE chromPeaks SET isotopesRatios=? WHERE id=?",
                (base64.b64encode(dumps(peak.isotopeRatios)).decode("utf-8"), peak.id),
            )
            curs.execute(
                "UPDATE allChromPeaks SET isotopesRatios=? WHERE id=?",
                (base64.b64encode(dumps(peak.isotopeRatios)).decode("utf-8"), peak.id),
            )

            curs.execute(
                "UPDATE chromPeaks SET mzDiffErrors=? WHERE id=?",
                (base64.b64encode(dumps(peak.mzDiffErrors)).decode("utf-8"), peak.id),
            )
            curs.execute(
                "UPDATE allChromPeaks SET mzDiffErrors=? WHERE id=?",
                (base64.b64encode(dumps(peak.mzDiffErrors)).decode("utf-8"), peak.id),
            )

        self.printMessage("%s: Annotating feature pairs done." % tracer.name, type="info")

        conn.commit()
        curs.close()
        conn.close()

    # annotate a metabolite group (consisting of ChromPeakPair instances) with the defined
    # hetero atoms based on detected feature pairs
    def annotateFeaturePairsWithHeteroAtoms(self, group, peaksInGroup):
        # iterate all peaks pairwise to find different adducts of the metabolite
        for pa in group:
            peakA = peaksInGroup[pa]

            for pb in group:
                peakB = peaksInGroup[pb]

                # peakA shall always have the lower mass
                if peakA.mz < peakB.mz and peakA.xCount == peakB.xCount:
                    ## check, if it could be a hetero atom
                    for pIso in self.heteroAtoms:
                        pIsotope = self.heteroAtoms[pIso]

                        bestFit = None
                        bestFitRatio = 100000
                        bestFitID = None

                        if abs(peakB.mz - peakA.mz - pIsotope.mzOffset) <= (max(peakA.mz, peakB.mz) * 2.5 * self.ppm / 1000000.0):
                            if pIsotope.mzOffset > 0:
                                ratio = peakB.LPeakArea / peakA.LPeakArea
                            else:
                                ratio = peakB.NPeakArea / peakA.NPeakArea

                            for nHAtoms in range(pIsotope.minCount, pIsotope.maxCount + 1):
                                if abs(ratio - nHAtoms * pIsotope.relativeAbundance) <= self.hAIntensityError:
                                    if abs(ratio - nHAtoms * pIsotope.relativeAbundance) < bestFitRatio:
                                        bestFitRatio = abs(ratio - nHAtoms * pIsotope.relativeAbundance)
                                        bestFit = Bunch(pIso=pIso, nHAtoms=nHAtoms)
                                        bestFitID = pb

                        if bestFit is not None:
                            peakB = peaksInGroup[bestFitID]
                            # peakA.heteroAtomsFeaturePairs.append(Bunch(name="M_%s%d"%(pIso, bestFit),  partnerAdd="_%s%d"%(pIso, bestFit),  toPeak=pb))
                            peakB.heteroAtomsFeaturePairs.append(
                                Bunch(
                                    name="_%s%d" % (bestFit.pIso, bestFit.nHAtoms),
                                    partnerAdd="M_%s%d" % (bestFit.pIso, bestFit.nHAtoms),
                                    toPeak=pa,
                                )
                            )

    # annotate a metabolite group (consisting of ChromPeakPair instances) with the defined
    # adducts and generate possible in-source fragments. Remove inconsistencies
    # in the form of impossible adduct combinations (e.g. [M+H]+ and [M+Na]+ for the same ion)
    def annotateChromPeaks(self, group, peaksInGroup):
        for pe in group:
            peak = peaksInGroup[pe]

            if not hasattr(peak, "fDesc"):
                setattr(peak, "fDesc", [])
            if not hasattr(peak, "adducts"):
                setattr(peak, "adducts", [])
            if not hasattr(peak, "Ms"):
                setattr(peak, "Ms", [])

        if len(group) <= 40:
            self.annotateFeaturePairsWithHeteroAtoms(group, peaksInGroup)
            for pa in group:
                peakA = peaksInGroup[pa]
                peakA.adducts.extend(peakA.heteroAtomsFeaturePairs)

            fT = formulaTools()

            # prepare adducts list
            addPol = {}
            addPol["+"] = []
            addPol["-"] = []
            adductsDict = {}
            for adduct in self.adducts:
                adductsDict[adduct.name] = adduct
                if adduct.polarity != "":
                    addPol[adduct.polarity].append(adduct)
            adducts = addPol
            adducts["+"] = sorted(adducts["+"], key=lambda x: x.mzoffset)
            adducts["-"] = sorted(adducts["-"], key=lambda x: x.mzoffset)

            # 1. iterate all peaks pairwise to find different adducts of the metabolite
            for pa in group:
                peakA = peaksInGroup[pa]
                for pb in group:
                    peakB = peaksInGroup[pb]

                    # peakA shall always have the lower mass
                    if peakA.mz < peakB.mz and (len(peakB.heteroAtomsFeaturePairs) == 0 or not all([s.name.startswith("_") for s in peakB.heteroAtomsFeaturePairs])):
                        # search for different adduct combinations
                        # different adducts must have the same number of labelled atoms
                        if peakA.xCount == peakB.xCount:
                            # check, if it could be some kind of adduct combination
                            for adA in adducts[peakA.ionMode]:
                                for adB in adducts[peakB.ionMode]:
                                    if adA.mCount == 1 and adB.mCount == 1 and adA.mzoffset < adB.mzoffset:
                                        if (
                                            peakA.ionMode == "-"
                                            and peakB.ionMode == "+"
                                            and peakA.loading == peakB.loading
                                            and abs(abs(adB.mzoffset - adA.mzoffset) - 2 * 1.007276) <= (max(peakA.mz, peakB.mz) * 2.5 * self.ppm / 1000000.0)
                                            and abs(peakB.mz - peakA.mz - 2 * 1.007276 / peakA.loading) <= (max(peakA.mz, peakB.mz) * 2.5 * self.ppm / 1000000.0)
                                        ):
                                            peakA.adducts.append(
                                                Bunch(
                                                    name="[M-H]-",
                                                    partnerAdd="[M+H]+",
                                                    toPeak=pb,
                                                )
                                            )
                                            peakB.adducts.append(
                                                Bunch(
                                                    name="[M+H]+",
                                                    partnerAdd="[M-H]-",
                                                    toPeak=pa,
                                                )
                                            )

                                        elif peakA.loading == adA.charge and peakB.loading == adB.charge and abs((peakA.mz - adA.mzoffset) / adA.mCount * peakA.loading - (peakB.mz - adB.mzoffset) / adB.mCount * peakB.loading) <= (max(peakA.mz, peakB.mz) * 2.5 * self.ppm / 1000000.0):
                                            peakA.adducts.append(
                                                Bunch(
                                                    name=adA.name,
                                                    partnerAdd=adB.name,
                                                    toPeak=pb,
                                                )
                                            )
                                            peakB.adducts.append(
                                                Bunch(
                                                    name=adB.name,
                                                    partnerAdd=adA.name,
                                                    toPeak=pa,
                                                )
                                            )

                        # search for adducts of the form [2M+XX]+-
                        elif peakA.xCount * 2 == peakB.xCount:
                            for adA in adducts[peakA.ionMode]:
                                for adB in adducts[peakB.ionMode]:
                                    if adA.mCount == 1 and adB.mCount == 2:
                                        if peakA.loading == adA.charge and peakB.loading == adB.charge and abs((peakA.mz - adA.mzoffset) / adA.mCount * peakA.loading - (peakB.mz - adB.mzoffset) / adB.mCount * peakB.loading) <= (max(peakA.mz, peakB.mz) * 2.5 * self.ppm / 1000000.0):
                                            peakA.adducts.append(
                                                Bunch(
                                                    name=adA.name,
                                                    partnerAdd=adB.name,
                                                    toPeak=pb,
                                                )
                                            )
                                            peakB.adducts.append(
                                                Bunch(
                                                    name=adB.name,
                                                    partnerAdd=adA.name,
                                                    toPeak=pa,
                                                )
                                            )

            def isAdductPairPresent(pA, adductAName, adductsA, pB, adductBName, adductsB, checkPartners=True):
                aFound = False
                for add in adductsA:
                    if add.name == adductAName and (not checkPartners or (add.partnerAdd == adductBName and add.toPeak == pB)):
                        aFound = True
                bFound = False
                for add in adductsB:
                    if add.name == adductBName and (not checkPartners or (add.partnerAdd == adductAName and add.toPeak == pA)):
                        bFound = True
                return aFound and bFound

            def removeAdductFromFeaturePair(pA, adductAName, adductsA):
                aFound = []
                for i, add in enumerate(adductsA):
                    if add.name == adductAName:
                        aFound.append(i)

                aFound = reversed(sorted(list(set(aFound))))
                for i in aFound:
                    del adductsA[i]

            def removeAdductPair(pA, adductAName, adductsA, pB, adductBName, adductsB):
                aFound = []
                for i, add in enumerate(adductsA):
                    if add.name == adductAName and add.partnerAdd == adductBName and add.toPeak == pB:
                        aFound.append(i)
                bFound = []
                for i, add in enumerate(adductsB):
                    if add.name == adductBName and add.partnerAdd == adductAName and add.toPeak == pA:
                        bFound.append(i)

                aFound = reversed(sorted(list(set(aFound))))
                for i in aFound:
                    del adductsA[i]
                bFound = reversed(sorted(list(set(bFound))))
                for i in bFound:
                    del adductsB[i]

            # 2. remove incorrect annotations from feature pairs e.g. A: [M+H]+ and [M+Na]+ with B: [M+Na]+ and [M+2Na-H]+
            for p in group:
                peak = peaksInGroup[p]
                peak.adducts = list(set(peak.adducts))

            for pa in group:
                peakA = peaksInGroup[pa]
                if len(peakA.adducts) > 0:
                    for pb in group:
                        peakB = peaksInGroup[pb]
                        if len(peakB.adducts) > 0:
                            if peakA.mz < peakB.mz:
                                for pc in group:
                                    peakC = peaksInGroup[pc]
                                    if len(peakC.adducts) > 0:
                                        if (
                                            isAdductPairPresent(
                                                pa,
                                                "[M+H]+",
                                                peakA.adducts,
                                                pb,
                                                "[M+Na]+",
                                                peakB.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+Na]+",
                                                peakA.adducts,
                                                pb,
                                                "[M+2Na-H]+",
                                                peakB.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+H]+",
                                                peakA.adducts,
                                                pc,
                                                "[M+2Na-H]+",
                                                peakC.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pb,
                                                "[M+H]+",
                                                peakB.adducts,
                                                pc,
                                                "[M+Na]+",
                                                peakC.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pb,
                                                "[M+Na]+",
                                                peakB.adducts,
                                                pc,
                                                "[M+2Na-H]+",
                                                peakC.adducts,
                                            )
                                            and peakA.mz < peakB.mz < peakC.mz
                                        ):
                                            removeAdductPair(
                                                pa,
                                                "[M+Na]+",
                                                peakA.adducts,
                                                pb,
                                                "[M+2Na-H]+",
                                                peakB.adducts,
                                            )
                                            removeAdductPair(
                                                pb,
                                                "[M+H]+",
                                                peakB.adducts,
                                                pc,
                                                "[M+Na]+",
                                                peakC.adducts,
                                            )

                                        if (
                                            isAdductPairPresent(
                                                pa,
                                                "[M+H]+",
                                                peakA.adducts,
                                                pb,
                                                "[M+K]+",
                                                peakB.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+K]+",
                                                peakA.adducts,
                                                pb,
                                                "[M+2K-H]+",
                                                peakB.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+H]+",
                                                peakA.adducts,
                                                pc,
                                                "[M+2K-H]+",
                                                peakC.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pb,
                                                "[M+H]+",
                                                peakB.adducts,
                                                pc,
                                                "[M+K]+",
                                                peakC.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pb,
                                                "[M+K]+",
                                                peakB.adducts,
                                                pc,
                                                "[M+2K-H]+",
                                                peakC.adducts,
                                            )
                                            and peakA.mz < peakB.mz < peakC.mz
                                        ):
                                            removeAdductPair(
                                                pa,
                                                "[M+K]+",
                                                peakA.adducts,
                                                pb,
                                                "[M+2K-H]+",
                                                peakB.adducts,
                                            )
                                            removeAdductPair(
                                                pb,
                                                "[M+H]+",
                                                peakB.adducts,
                                                pc,
                                                "[M+K]+",
                                                peakC.adducts,
                                            )

                                        if (
                                            isAdductPairPresent(
                                                pa,
                                                "[M+2H]++",
                                                peakA.adducts,
                                                pb,
                                                "[M+H+Na]++",
                                                peakB.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+H+Na]++",
                                                peakA.adducts,
                                                pb,
                                                "[M+2Na]++",
                                                peakB.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+H]++",
                                                peakA.adducts,
                                                pc,
                                                "[M+2Na]++",
                                                peakC.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pb,
                                                "[M+2H]++",
                                                peakB.adducts,
                                                pc,
                                                "[M+H+Na]++",
                                                peakC.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pb,
                                                "[M+H+Na]++",
                                                peakB.adducts,
                                                pc,
                                                "[M+2Na]++",
                                                peakC.adducts,
                                            )
                                            and peakA.mz < peakB.mz < peakC.mz
                                        ):
                                            removeAdductPair(
                                                pa,
                                                "[M+H+Na]++",
                                                peakA.adducts,
                                                pb,
                                                "[M+2Na]++",
                                                peakB.adducts,
                                            )
                                            removeAdductPair(
                                                pb,
                                                "[M+2H]++",
                                                peakB.adducts,
                                                pc,
                                                "[M+H+Na]++",
                                                peakC.adducts,
                                            )

                                        if (
                                            isAdductPairPresent(
                                                pa,
                                                "[M+2H]++",
                                                peakA.adducts,
                                                pb,
                                                "[M+H+K]++",
                                                peakB.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+H+K]++",
                                                peakA.adducts,
                                                pb,
                                                "[M+2K]++",
                                                peakB.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+H]++",
                                                peakA.adducts,
                                                pc,
                                                "[M+2K]++",
                                                peakC.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pb,
                                                "[M+2H]++",
                                                peakB.adducts,
                                                pc,
                                                "[M+H+K]++",
                                                peakC.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pb,
                                                "[M+H+K]++",
                                                peakB.adducts,
                                                pc,
                                                "[M+2K]++",
                                                peakC.adducts,
                                            )
                                            and peakA.mz < peakB.mz < peakC.mz
                                        ):
                                            removeAdductPair(
                                                pa,
                                                "[M+H+K]++",
                                                peakA.adducts,
                                                pb,
                                                "[M+2K]++",
                                                peakB.adducts,
                                            )
                                            removeAdductPair(
                                                pb,
                                                "[M+2H]++",
                                                peakB.adducts,
                                                pc,
                                                "[M+H+K]++",
                                                peakC.adducts,
                                            )

                                        if (
                                            isAdductPairPresent(
                                                pb,
                                                "[M+H]+",
                                                peakB.adducts,
                                                pc,
                                                "[M+Na]+",
                                                peakC.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+2Na-H]+",
                                                peakA.adducts,
                                                pc,
                                                "[M+CH3FeN]+",
                                                peakC.adducts,
                                            )
                                            and peakA.mz < peakB.mz < peakC.mz
                                        ):
                                            removeAdductPair(
                                                pa,
                                                "[M+2Na-H]+",
                                                peakA.adducts,
                                                pc,
                                                "[M+CH3FeN]+",
                                                peakC.adducts,
                                            )

                                        if (
                                            isAdductPairPresent(
                                                pa,
                                                "[M-H]-",
                                                peakA.adducts,
                                                pb,
                                                "[M+H]+",
                                                peakB.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+Na-2H]-",
                                                peakA.adducts,
                                                pc,
                                                "[M+Na]+",
                                                peakC.adducts,
                                            )
                                            and peakA.mz < peakB.mz < peakC.mz
                                        ):
                                            aFound = []
                                            for i, add in enumerate(peakA.adducts):
                                                if add.name == "[M+Na-2H]-":
                                                    aFound.append(i)
                                            aFound = reversed(sorted(list(set(aFound))))
                                            for i in aFound:
                                                del peakA.adducts[i]

                                        if (
                                            isAdductPairPresent(
                                                pa,
                                                "[M-H]-",
                                                peakA.adducts,
                                                pb,
                                                "[M+H]+",
                                                peakB.adducts,
                                            )
                                            and isAdductPairPresent(
                                                pa,
                                                "[M+Na]+",
                                                peakA.adducts,
                                                pc,
                                                "[M+Na]+",
                                                peakC.adducts,
                                            )
                                            and peakA.mz < peakB.mz < peakC.mz
                                        ):
                                            bFound = []
                                            for i, add in enumerate(peakB.adducts):
                                                if add.name == "[M+Na-2H]-":
                                                    bFound.append(i)
                                            bFound = reversed(sorted(list(set(bFound))))
                                            for i in bFound:
                                                del peakB.adducts[i]

            for pa in group:
                peakA = peaksInGroup[pa]
                if len(peakA.adducts) > 0:
                    for pb in group:
                        peakB = peaksInGroup[pb]
                        if len(peakB.adducts) > 0:
                            if peakA.mz < peakB.mz:
                                if (
                                    isAdductPairPresent(
                                        pa,
                                        "[M+H]+",
                                        peakA.adducts,
                                        pb,
                                        "[M+Na]+",
                                        peakB.adducts,
                                    )
                                    and isAdductPairPresent(
                                        pa,
                                        "[M+Na]+",
                                        peakA.adducts,
                                        pb,
                                        "[M+2Na-H]+",
                                        peakB.adducts,
                                    )
                                    and peakA.mz < peakB.mz
                                ):
                                    removeAdductPair(
                                        pa,
                                        "[M+Na]+",
                                        peakA.adducts,
                                        pb,
                                        "[M+2Na-H]+",
                                        peakB.adducts,
                                    )

                                if (
                                    isAdductPairPresent(
                                        pa,
                                        "[M+H]+",
                                        peakA.adducts,
                                        pb,
                                        "[M+K]+",
                                        peakB.adducts,
                                    )
                                    and isAdductPairPresent(
                                        pa,
                                        "[M+K]+",
                                        peakA.adducts,
                                        pb,
                                        "[M+2K-H]+",
                                        peakB.adducts,
                                    )
                                    and peakA.mz < peakB.mz
                                ):
                                    removeAdductPair(
                                        pa,
                                        "[M+K]+",
                                        peakA.adducts,
                                        pb,
                                        "[M+2K-H]+",
                                        peakB.adducts,
                                    )

                                if (
                                    isAdductPairPresent(
                                        pa,
                                        "[M-H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+Na]+",
                                        peakB.adducts,
                                    )
                                    and isAdductPairPresent(
                                        pa,
                                        "[M+Na-2H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+2Na-H]+",
                                        peakB.adducts,
                                    )
                                    and peakA.mz < peakB.mz
                                ):
                                    removeAdductPair(
                                        pa,
                                        "[M+Na-2H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+2Na-H]+",
                                        peakB.adducts,
                                    )

                                if (
                                    isAdductPairPresent(
                                        pa,
                                        "[M-H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+K]+",
                                        peakB.adducts,
                                    )
                                    and isAdductPairPresent(
                                        pa,
                                        "[M+K-2H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+2K-H]+",
                                        peakB.adducts,
                                    )
                                    and peakA.mz < peakB.mz
                                ):
                                    removeAdductPair(
                                        pa,
                                        "[M+K-2H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+2K-H]+",
                                        peakB.adducts,
                                    )

                                if (
                                    isAdductPairPresent(
                                        pa,
                                        "[M-H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+Na]+",
                                        peakB.adducts,
                                    )
                                    and isAdductPairPresent(
                                        pa,
                                        "[M+Na-2H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+H]+",
                                        peakB.adducts,
                                    )
                                    and peakA.mz < peakB.mz
                                ):
                                    removeAdductPair(
                                        pa,
                                        "[M+Na-2H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+Na]+",
                                        peakB.adducts,
                                    )

                                if (
                                    isAdductPairPresent(
                                        pa,
                                        "[M-H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+K]+",
                                        peakB.adducts,
                                    )
                                    and isAdductPairPresent(
                                        pa,
                                        "[M+K-2H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+H]+",
                                        peakB.adducts,
                                    )
                                    and peakA.mz < peakB.mz
                                ):
                                    removeAdductPair(
                                        pa,
                                        "[M+K-2H]-",
                                        peakA.adducts,
                                        pb,
                                        "[M+K]+",
                                        peakB.adducts,
                                    )

                                if (
                                    isAdductPairPresent(
                                        pa,
                                        "[M+2H]++",
                                        peakA.adducts,
                                        pb,
                                        "[M+H+Na]++",
                                        peakB.adducts,
                                    )
                                    and isAdductPairPresent(
                                        pa,
                                        "[M+H+Na]++",
                                        peakA.adducts,
                                        pb,
                                        "[M+2Na]++",
                                        peakB.adducts,
                                    )
                                    and peakA.mz < peakB.mz
                                ):
                                    removeAdductPair(
                                        pa,
                                        "[M+H+Na]++",
                                        peakA.adducts,
                                        pb,
                                        "[M+2Na]++",
                                        peakB.adducts,
                                    )

                                if (
                                    isAdductPairPresent(
                                        pa,
                                        "[M+2H]++",
                                        peakA.adducts,
                                        pb,
                                        "[M+H+K]++",
                                        peakB.adducts,
                                    )
                                    and isAdductPairPresent(
                                        pa,
                                        "[M+H+K]++",
                                        peakA.adducts,
                                        pb,
                                        "[M+2K]++",
                                        peakB.adducts,
                                    )
                                    and peakA.mz < peakB.mz
                                ):
                                    removeAdductPair(
                                        pa,
                                        "[M+H+Na]++",
                                        peakA.adducts,
                                        pb,
                                        "[M+2Na]++",
                                        peakB.adducts,
                                    )

            for pa in group:
                peakA = peaksInGroup[pa]
                if len(peakA.adducts) > 0:
                    if isAdductPairPresent(
                        pa,
                        "[M+H]+",
                        peakA.adducts,
                        pa,
                        "[2M+H]+",
                        peakA.adducts,
                        checkPartners=False,
                    ):
                        removeAdductFromFeaturePair(pa, "[M+H]+", peakA.adducts)
                    if isAdductPairPresent(
                        pa,
                        "[M+NH4]+",
                        peakA.adducts,
                        pa,
                        "[2M+NH4]+",
                        peakA.adducts,
                        checkPartners=False,
                    ):
                        removeAdductFromFeaturePair(pa, "[M+NH4]+", peakA.adducts)
                    if isAdductPairPresent(
                        pa,
                        "[M+Na]+",
                        peakA.adducts,
                        pa,
                        "[2M+Na]+",
                        peakA.adducts,
                        checkPartners=False,
                    ):
                        removeAdductFromFeaturePair(pa, "[M+Na]+", peakA.adducts)
                    if isAdductPairPresent(
                        pa,
                        "[M+K]+",
                        peakA.adducts,
                        pa,
                        "[2M+K]+",
                        peakA.adducts,
                        checkPartners=False,
                    ):
                        removeAdductFromFeaturePair(pa, "[M+K]+", peakA.adducts)

            inSourceFragments = {}

            # 3. iterate all peaks pairwise to find common in-source fragments
            for pa in group:
                peakA = peaksInGroup[pa]

                for pb in group:
                    peakB = peaksInGroup[pb]

                    # peakA shall always have the lower mass
                    if peakA.mz < peakB.mz and abs(peakA.mz - peakB.mz) <= 101:
                        mzDif = peakB.mz - peakA.mz
                        done = False

                        # generate putative in-source fragments (using the labelled carbon atoms)
                        if len(self.elements) > 0:
                            elemDictReq = {}
                            for elem in self.elements:
                                elemDictReq[elem] = [
                                    self.elements[elem].weight,
                                    self.elements[elem].numberValenzElectrons,
                                ]
                            t = SGRGenerator(atoms=elemDictReq)

                            useAtoms = []
                            if self.labellingElement in self.elements.keys():
                                useAtoms.append(self.labellingElement)
                                for k in self.elements.keys():
                                    if k != self.labellingElement:
                                        useAtoms.append(k)
                            else:
                                useAtoms = list(self.elements.keys())

                            if not (isinstance(peakA.xCount, str)):
                                atomsRange = []
                                if self.labellingElement in useAtoms:
                                    if self.metabolisationExperiment:
                                        elem = self.elements[self.labellingElement]
                                        atomsRange.append(
                                            [
                                                abs(peakA.xCount - peakB.xCount),
                                                abs(peakA.xCount - peakB.xCount + elem.maxCount - elem.minCount),
                                            ]
                                        )
                                    else:
                                        atomsRange.append(abs(peakA.xCount - peakB.xCount))

                                for elem in self.elements.keys():
                                    if elem != self.labellingElement:
                                        elem = self.elements[elem]
                                        atomsRange.append([elem.minCount, elem.maxCount])

                                corrFact = abs(peakB.mz - peakA.mz)
                                if corrFact <= 1:
                                    corrFact = 1.0
                                pFs = t.findFormulas(
                                    mzDif,
                                    ppm=self.ppm * 2.0 * peakA.mz / corrFact,
                                    useAtoms=useAtoms,
                                    atomsRange=atomsRange,
                                    fixed=self.labellingElement,
                                    useSevenGoldenRules=False,
                                )

                                for pF in pFs:
                                    if pa not in inSourceFragments.keys():
                                        inSourceFragments[pa] = {}
                                    if pb not in inSourceFragments[pa].keys():
                                        inSourceFragments[pa][pb] = []

                                    c = fT.parseFormula(pF)
                                    mw = fT.calcMolWeight(c)
                                    dif = abs(abs(peakB.mz - peakA.mz) - mw)
                                    sf = fT.flatToString(c, prettyPrintWithHTMLTags=False)
                                    inSourceFragments[pa][pb].append("%.4f-%s" % (peakB.mz, sf))

            if self.simplifyInSourceFragments:
                for pa in group:
                    peakA = peaksInGroup[pa]
                    peakA.fDesc = []

                    for pb in group:
                        peakB = peaksInGroup[pb]

                        for pc in group:
                            peakC = peaksInGroup[pc]

                            if peakA.mz < peakC.mz < peakB.mz:
                                if pa in inSourceFragments.keys() and pb in inSourceFragments[pa].keys() and pc in inSourceFragments[pa].keys() and pc in inSourceFragments.keys() and pb in inSourceFragments[pc].keys():
                                    del inSourceFragments[pa][pc]

            for pa in group:
                peakA = peaksInGroup[pa]

                if pa in inSourceFragments.keys():
                    for pb in inSourceFragments[pa].keys():
                        for inFrag in inSourceFragments[pa][pb]:
                            peakA.fDesc.append(inFrag)

            for pe in group:
                peak = peaksInGroup[pe]
                peak.fDesc = list(set(peak.fDesc))
                peak.adducts = list(set([a.name for a in peak.adducts]))

                if not hasattr(peak, "Ms"):
                    setattr(peak, "Ms", [])

                peak.Ms = []
                for assignedAdduct in peak.adducts:
                    if assignedAdduct in adductsDict.keys():
                        peak.Ms.append((peak.mz - adductsDict[assignedAdduct].mzoffset) / adductsDict[assignedAdduct].mCount / peak.loading)

    # data processing step 6: convolute different feature pairs into feature groups using the chromatographic
    # profiles of different metabolite ions
    def groupFeaturePairsUntargetedAndWriteToDB(self, chromPeaks, mzxml, tracer, tracerID, reportFunction=None):
        try:
            conn = connect(self.file + getDBSuffix(), isolation_level="DEFERRED")
            conn.execute("""PRAGMA synchronous = OFF""")
            conn.execute("""PRAGMA journal_mode = OFF""")
            curs = conn.cursor()

            nodes = {}
            correlations = {}

            allPeaks = {}
            for peak in chromPeaks:
                nodes[peak.id] = []
                allPeaks[peak.id] = peak
                peak.correlationsToOthers = []

            rengine = ro.r

            ## export model for R/3.5.3 from outside of MetExtract
            ## in R load the model and save it as the varible model
            ## execute
            ## load("C:/Users/cbueschl/Desktop/PeakCorrelation/model_main_2.Rdata")
            ## exp = list(predict = model$modelInfo$predict, finalModel = model$finalModel)
            ## save(exp, file = "C:/Users/cbueschl/Desktop/PeakCorrelation/model.Rdata")
            # rengine('load("C:/Users/cbueschl/Desktop/PeakCorrelation/model.Rdata")')
            # rengine('model = exp')
            # rengine('print("model is")')
            # rengine('print(model)')

            # compare all detected feature pairs at approximately the same retention time
            for piA in range(len(chromPeaks)):
                peakA = chromPeaks[piA]
                if reportFunction is not None:
                    reportFunction(
                        0.7 * piA / len(chromPeaks),
                        "Matching features (%d remaining)" % (len(chromPeaks) - piA),
                    )

                if peakA.id not in correlations.keys():
                    correlations[peakA.id] = {}

                for peakB in chromPeaks:
                    if peakB.id not in correlations.keys():
                        correlations[peakB.id] = {}

                    if peakA.mz < peakB.mz:
                        if abs(peakA.NPeakCenter - peakB.NPeakCenter) < self.peakCenterError:
                            bmin = int(
                                max(
                                    0,
                                    mean(
                                        [
                                            peakA.NPeakCenter - 1 * peakA.NBorderLeft,
                                            peakB.NPeakCenter - 1 * peakB.NBorderLeft,
                                            peakA.LPeakCenter - 1 * peakA.LBorderLeft,
                                            peakB.LPeakCenter - 1 * peakB.LBorderLeft,
                                        ]
                                    ),
                                )
                            )
                            bmax = int(
                                min(
                                    len(peakB.NXICSmoothed) - 1,
                                    mean(
                                        [
                                            peakB.NPeakCenter + 1 * peakB.NBorderRight,
                                            peakA.NPeakCenter + 1 * peakA.NBorderRight,
                                            peakB.LPeakCenter + 1 * peakB.LBorderRight,
                                            peakA.LPeakCenter + 1 * peakA.LBorderRight,
                                        ]
                                    ),
                                )
                            )

                            pb = corr(
                                peakA.NXICSmoothed[bmin:bmax],
                                peakB.NXICSmoothed[bmin:bmax],
                            )

                            if False:
                                ## development: test for shifted correlation, Kristina Missbach
                                testShifts = np.arange(-3, 3.00001, 0.1)
                                nonShiftCorr, shiftBestCorr, bestShift = getCorrelationShifted(
                                    np.transpose(np.array([peakA.times, peakA.NXICSmoothed]))[bmin:bmax],
                                    np.transpose(np.array([peakB.times, peakB.NXICSmoothed]))[bmin:bmax],
                                    peakA.times[bmin],
                                    peakA.times[bmax],
                                    rtShiftsMin=testShifts,
                                )
                                assignedClass = "NA"
                                try:
                                    ret = rengine(
                                        "as.character(model$predict(model$finalModel, data.frame(nonShiftedCorr = %f, bestShiftedCorr = %f, bestShift = %f)))"
                                        % (
                                            nonShiftCorr,
                                            shiftBestCorr,
                                            bestShift / 60.0,
                                        )
                                    )
                                    assignedClass = ret[0]
                                except Exception as ex:
                                    print("Exception: %s" % (ex.message))
                                print(
                                    "Correlation unshifted is %.3f/%.3f, best shifted correlation is %.3f by %.3f sec, verdict %s"
                                    % (
                                        pb,
                                        nonShiftCorr,
                                        shiftBestCorr,
                                        bestShift,
                                        assignedClass,
                                    )
                                )

                            if str(pb) == "nan":
                                pb = -1

                            correlations[peakA.id][peakB.id] = pb
                            correlations[peakB.id][peakA.id] = pb

                            silRatiosA = peakA.silRatios.silRatios
                            silRatiosB = peakB.silRatios.silRatios

                            meanSilRatioA = weightedMean(silRatiosA, peakA.silRatios.peakNInts)
                            meanSilRatioB = weightedMean(silRatiosB, peakB.silRatios.peakNInts)

                            silRatiosFold = 0
                            silRatiosSD = 1

                            try:
                                silRatiosFold = max(meanSilRatioA, meanSilRatioB) / min(meanSilRatioA, meanSilRatioB)
                                silRatiosSD = weightedSd(
                                    [abs(r - meanSilRatioA) / meanSilRatioA for r in silRatiosA if min(r, meanSilRatioA) > 0] + [abs(r - meanSilRatioB) / meanSilRatioB for r in silRatiosB if min(r, meanSilRatioB) > 0],
                                    peakA.silRatios.peakNInts + peakB.silRatios.peakNInts,
                                )

                                # check for similar chromatographic peak profile and similar native to labeled ratio
                                if pb >= self.minCorrelation and silRatiosFold <= 1 + max(0.25, 3 * silRatiosSD):
                                    nodes[peakA.id].append(peakB.id)
                                    nodes[peakB.id].append(peakA.id)

                                    peakB.correlationsToOthers.append(peakA.id)
                                    peakA.correlationsToOthers.append(peakB.id)

                            except Exception as e:
                                logging.error("Error while convoluting two feature pairs, skipping.. (%s)" % str(e))

                            try:
                                SQLInsert(
                                    curs,
                                    "featurefeatures",
                                    fID1=peakA.id,
                                    fID2=peakB.id,
                                    corr=pb,
                                    silRatioValue=silRatiosFold,
                                )
                            except Exception as e:
                                logging.error("Error while convoluting two feature pairs, skipping.. (%s)" % str(e))
                                SQLInsert(
                                    curs,
                                    "featurefeatures",
                                    fID1=peakA.id,
                                    fID2=peakB.id,
                                    corr=0,
                                    silRatioValue=0,
                                )

            self.postMessageToProgressWrapper("text", "%s: Convoluting feature groups" % tracer.name)

            for peak in chromPeaks:
                delattr(peak, "NXIC")
                delattr(peak, "LXIC")
                delattr(peak, "times")

            for k in nodes.keys():
                uniq = []
                for u in nodes[k]:
                    if u not in uniq:
                        uniq.append(u)
                nodes[k] = uniq

            # get subgraphs from the feature pair graph. Each subgraph represents one convoluted
            # feature group
            tGroups = getSubGraphs(nodes)

            def splitGroupWithHCA(tGroup, correlations):
                # if 1 or two feature pairs are in a group, automatically use them (no further splitting possible)
                if len(tGroup) <= 2:
                    return [tGroup]

                # construct 2-dimensional data matrix
                data = []
                for tG1 in tGroup:
                    t = []
                    for tG2 in tGroup:
                        if tG1 in correlations.keys() and tG2 in correlations[tG1].keys():
                            t.append(correlations[tG1][tG2])
                        elif tG1 == tG2:
                            t.append(1.0)
                        else:
                            t.append(0.0)
                    data.append(t)

                # calculate HCA tree using the correlations between all feature pairs
                hc = HCA_general.HCA_generic()
                tree = hc.generateTree(objs=data, ids=tGroup)
                # hc.plotTree(tree)
                # print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")

                # split the HCA tree in sub-clusters
                def checkSubCluster(tree, hca, corrThreshold, cutOffMinRatio):
                    if isinstance(tree, HCA_general.HCALeaf):
                        return False
                    elif isinstance(tree, HCA_general.HCAComposite):
                        corrsLKid = hca.link(tree.getLeftKid())
                        corrs = hca.link(tree)
                        corrsRKid = hca.link(tree.getRightKid())

                        aboveThresholdLKid = sum([corr > corrThreshold for corr in corrsLKid])
                        aboveThreshold = sum([corr > corrThreshold for corr in corrs])
                        aboveThresholdRKid = sum([corr > corrThreshold for corr in corrsRKid])

                        # print(aboveThresholdLKid, aboveThreshold, aboveThresholdRKid)

                    if (aboveThresholdLKid * 1.0 / len(corrs)) >= cutOffMinRatio and (aboveThreshold * 1.0 / len(corrs)) >= cutOffMinRatio and (aboveThresholdRKid * 1.0 / len(corrs)) >= cutOffMinRatio:
                        return False
                    else:
                        return True

                # subClusts=hc.splitTreeWithCallbackBottomUp(tree,
                #                                   CallBackMethod(_target=checkSubCluster, corrThreshold=self.minCorrelation, cutOffMinRatio=self.minCorrelationConnections).getRunMethod())
                subClusts = hc.splitTreeWithCallback(
                    tree,
                    CallBackMethod(
                        _target=checkSubCluster,
                        corrThreshold=self.minCorrelation,
                        cutOffMinRatio=self.minCorrelationConnections,
                    ).getRunMethod(),
                    recursive=False,
                )

                # convert the subclusters into arrays of feature pairs belonging to the same metabolite
                return [[leaf.getID() for leaf in subClust.getLeaves()] for subClust in subClusts]

            groups = []
            done = 0
            for tGroup in tGroups:
                # print("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
                cGroups = [tGroup]
                while len(cGroups) > 0:
                    # print("HCA with", len(cGroups))
                    gGroup = cGroups.pop(0)

                    ## TODO optimize this code, it recalculates the computationally expensive HCA too often for a high number of features
                    if False and len(gGroup) > 100:
                        groups.append(gGroup)
                        continue

                    sGroups = splitGroupWithHCA(gGroup, correlations)
                    # print("first group", len(gGroup), gGroup, "split into", sGroups)

                    if len(sGroups) == 1:
                        groups.append(sGroups[0])
                    else:
                        cGroups.extend(sGroups)

                done = done + 1
                self.postMessageToProgressWrapper(
                    "text",
                    "%s: Convoluting feature groups (%d/%d done)" % (tracer.name, done, len(tGroups)),
                )

            if True:
                done = 0
                for gi in range(len(groups)):
                    group = groups[gi]
                    if reportFunction is not None:
                        reportFunction(
                            0.7 + 0.3 * piA / len(chromPeaks),
                            "Searching for relationships (%d remaining)" % (len(groups) - gi),
                        )

                    # first, search for adduct relationships (same number of carbon atoms) in each convoluted
                    # feature group

                    peaksInGroup = {}
                    for a in group:
                        peaksInGroup[a] = allPeaks[a]

                        temp = []
                        for j in range(len(peaksInGroup[a].correlationsToOthers)):
                            if peaksInGroup[a].correlationsToOthers[j] in group:
                                temp.append(peaksInGroup[a].correlationsToOthers[j])

                        peaksInGroup[a].correlationsToOthers = temp

                    self.annotateChromPeaks(group, peaksInGroup)  # store feature pair annotation in the database

                    done = done + 1
                    self.postMessageToProgressWrapper(
                        "text",
                        "%s: Annotating feature groups (%d/%d done)" % (tracer.name, done, len(groups)),
                    )

            for peak in chromPeaks:
                adds = countEntries(peak.adducts)
                peak.adducts = list(adds.keys())
                curs.execute(
                    "UPDATE chromPeaks SET adducts=?, fDesc=?, correlationsToOthers=? WHERE id=?",
                    (
                        base64.b64encode(dumps(peak.adducts)).decode("utf-8"),
                        base64.b64encode(dumps(peak.fDesc)).decode("utf-8"),
                        base64.b64encode(dumps(peak.correlationsToOthers)).decode("utf-8"),
                        peak.id,
                    ),
                )
                curs.execute(
                    "UPDATE allChromPeaks SET adducts=?, fDesc=? WHERE id=?",
                    (
                        base64.b64encode(dumps(peak.adducts)).decode("utf-8"),
                        base64.b64encode(dumps(peak.fDesc)).decode("utf-8"),
                        peak.id,
                    ),
                )

                curs.execute(
                    "UPDATE chromPeaks SET heteroAtomsFeaturePairs=? WHERE id=?",
                    (
                        base64.b64encode(dumps(peak.heteroAtomsFeaturePairs)).decode("utf-8"),
                        peak.id,
                    ),
                )
                curs.execute(
                    "UPDATE allChromPeaks SET heteroAtomsFeaturePairs=? WHERE id=?",
                    (base64.b64encode(dumps(peak.heteroAtomsFeaturePairs)), peak.id),
                )

            # store feature group in the database
            for group in sorted(
                groups,
                key=lambda x: sum([allPeaks[p].NPeakCenterMin / 60.0 for p in x]) / len(x),
            ):
                SQLInsert(
                    curs,
                    "featureGroups",
                    id=self.curFeatureGroupID,
                    featureName="fg_%d" % self.curFeatureGroupID,
                    tracer=tracerID,
                )
                groupMeanElutionIndex = 0

                hasPos = False
                hasNeg = False

                for p in sorted(group, key=lambda x: allPeaks[x].mz):
                    hasPos = hasPos or allPeaks[p].ionMode == "+"
                    hasNeg = hasNeg or allPeaks[p].ionMode == "-"

                    groupMeanElutionIndex += allPeaks[p].NPeakCenter

                    allPeaks[p].fGroupID = self.curFeatureGroupID
                    SQLInsert(
                        curs,
                        "featureGroupFeatures",
                        fID=p,
                        fDesc="",
                        fGroupID=self.curFeatureGroupID,
                    )

                groupMeanElutionIndex = groupMeanElutionIndex / len(group)

                # store one positve and one negative ionisation mode MS scan in the database for
                # later visualisation (one for each convoluted feature group)
                if hasPos:
                    scan = mzxml.getIthMS1Scan(int(groupMeanElutionIndex), self.positiveScanEvent)

                    SQLInsert(
                        curs,
                        "massspectrum",
                        mID=self.curMassSpectrumID,
                        fgID=self.curFeatureGroupID,
                        time=scan.retention_time,
                        mzs=";".join([str(u) for u in scan.mz_list]),
                        intensities=";".join([str(u) for u in scan.intensity_list]),
                        ionMode="+",
                    )
                    self.curMassSpectrumID += 1
                if hasNeg:
                    scan = mzxml.getIthMS1Scan(int(groupMeanElutionIndex), self.negativeScanEvent)

                    SQLInsert(
                        curs,
                        "massspectrum",
                        mID=self.curMassSpectrumID,
                        fgID=self.curFeatureGroupID,
                        time=scan.retention_time,
                        mzs=";".join([str(u) for u in scan.mz_list]),
                        intensities=";".join([str(u) for u in scan.intensity_list]),
                        ionMode="-",
                    )
                    self.curMassSpectrumID += 1

                self.curFeatureGroupID += 1

            conn.commit()
            self.printMessage(
                "%s: Feature grouping done. " % tracer.name + str(len(groups)) + " feature groups",
                type="info",
            )

            conn.commit()
            curs.close()
            conn.close()

        except Exception as ex:
            import traceback

            traceback.print_exc()

            self.printMessage("Error in %s: %s" % (self.file, str(ex)), type="error")
            self.postMessageToProgressWrapper("failed", self.pID)

    # store one MS scan for each detected feature pair in the database
    def writeMassSpectraToDB(self, chromPeaks, mzxml, reportFunction=None):
        conn = connect(self.file + getDBSuffix(), isolation_level="DEFERRED")
        conn.execute("""PRAGMA synchronous = OFF""")
        conn.execute("""PRAGMA journal_mode = OFF""")
        curs = conn.cursor()
        massSpectraWrittenPos = {}
        massSpectraWrittenNeg = {}
        for pi in range(len(chromPeaks)):
            peak = chromPeaks[pi]
            if reportFunction is not None:
                reportFunction(
                    1.0 * pi / len(chromPeaks),
                    "%d feature pairs remaining" % (len(chromPeaks) - pi),
                )

            scanEvent = ""
            iMode = ""
            massSpectraWritten = None
            if peak.ionMode == "+":
                iMode = "+"
                scanEvent = self.positiveScanEvent
                massSpectraWritten = massSpectraWrittenPos
            elif peak.ionMode == "-":
                iMode = "-"
                scanEvent = self.negativeScanEvent
                massSpectraWritten = massSpectraWrittenNeg

            if peak.NPeakCenter not in massSpectraWritten.keys():
                scan = mzxml.getIthMS1Scan(peak.NPeakCenter, scanEvent)

                SQLInsert(
                    curs,
                    "massspectrum",
                    mID=self.curMassSpectrumID,
                    fgID=-1,
                    time=scan.retention_time,
                    mzs=";".join([str(u) for u in scan.mz_list]),
                    intensities=";".join([str(u) for u in scan.intensity_list]),
                    ionMode=iMode,
                )

                massSpectraWritten[peak.NPeakCenter] = self.curMassSpectrumID
                self.curMassSpectrumID = self.curMassSpectrumID + 1

            curs.execute("UPDATE chromPeaks SET massSpectrumID=%d WHERE id=%d" % (massSpectraWritten[peak.NPeakCenter], peak.id))

        conn.commit()
        curs.close()
        conn.close()

    ## write a new featureML file
    def writeResultsToFeatureML(self, forFile):
        conn = connect(forFile + getDBSuffix(), isolation_level="DEFERRED")
        conn.execute("""PRAGMA synchronous = OFF""")
        conn.execute("""PRAGMA journal_mode = OFF""")
        curs = conn.cursor()

        features = []

        for chromPeak in SQLSelectAsObject(
            curs,
            "SELECT c.id AS id, c.mz AS mz, c.lmz AS lmz, c.xcount AS xCount, c.Loading AS loading, c.NPeakCenterMin AS NPeakCenterMin, c.ionMode AS ionMode FROM chromPeaks c",
            newObject=ChromPeakPair,
        ):
            b = Bunch(
                id=chromPeak.id,
                ogroup="-1",
                mz=chromPeak.mz,
                rt=chromPeak.NPeakCenterMin,
                Xn=chromPeak.xCount,
                lmz=chromPeak.lmz,
                charge=chromPeak.loading,
                name=chromPeak.id,
                ionMode=chromPeak.ionMode,
            )
            features.append(b)

        exportAsFeatureML.writeFeatureListToFeatureML(features, forFile + ".featureML", ppmPM=self.ppm, rtPM=0.25 * 60)

    # write feature pairs detected in this LC-HRMS data file into a new TSV file.
    # Each row represent one feature pair
    def writeResultsToTSVFile(self, forFile):
        conn = connect(forFile + getDBSuffix(), isolation_level="DEFERRED")
        conn.execute("""PRAGMA synchronous = OFF""")
        conn.execute("""PRAGMA journal_mode = OFF""")
        curs = conn.cursor()

        chromPeaks = []
        configTracers = {}

        for chromPeak in SQLSelectAsObject(
            curs,
            "SELECT c.id AS id, g.fGroupID AS fGroupID, c.mz AS mz, c.lmz AS lmz, c.tmz AS tmz, c.xcount AS xCount, c.Loading AS loading, "
            "c.ionMode AS ionMode, c.NPeakCenter AS NPeakCenter, c.NPeakCenterMin AS NPeakCenterMin, "
            "c.NPeakScale AS NPeakScale, c.NPeakArea AS NPeakArea, c.NPeakAbundance AS NPeakAbundance, c.LPeakCenter AS LPeakCenter, "
            "c.LPeakCenterMin AS LPeakCenterMin, c.LPeakScale AS LPeakScale, c.LPeakArea AS LPeakArea, c.LPeakAbundance AS LPeakAbundance, "
            "c.NBorderLeft as NBorderLeft, c.NBorderRight as NBorderRight, c.LBorderLeft as LBorderLeft, c.LBorderRight as LBorderRight, "
            "c.peaksCorr AS peaksCorr, c.assignedMZs AS assignedMZs, c.heteroAtoms AS HAs, c.adducts AS ADs, "
            "c.fDesc AS DSc, t.id AS tracer, t.name AS tracerName, "
            "c.peaksRatio AS peaksRatio, c.peaksRatioMp1 AS peaksRatioMp1, c.peaksRatioMPm1 as peaksRatioMPm1, "
            "c.isotopesRatios AS isotopeRatios , c.mzDiffErrors AS mzDiffErrors , c.comments AS comments, c.artificialEICLShift AS artificialEICLShift FROM "
            "chromPeaks c INNER JOIN featureGroupFeatures g ON c.id=g.fID INNER JOIN tracerConfiguration t ON c.tracer=t.id",
            newObject=ChromPeakPair,
        ):
            chromPeak.NPeakCenterMin = chromPeak.NPeakCenterMin / 60.0
            chromPeak.LPeakCenterMin = chromPeak.LPeakCenterMin / 60.0
            setattr(chromPeak, "heteroIsotopologues", loads(base64.b64decode(chromPeak.HAs)))
            setattr(chromPeak, "adducts", loads(base64.b64decode(chromPeak.ADs)))
            setattr(chromPeak, "fDesc", loads(base64.b64decode(chromPeak.DSc)))
            setattr(
                chromPeak,
                "isotopeRatios",
                loads(base64.b64decode(chromPeak.isotopeRatios)),
            )
            setattr(
                chromPeak,
                "mzDiffErrors",
                loads(base64.b64decode(chromPeak.mzDiffErrors)),
            )
            setattr(
                chromPeak,
                "comments",
                "; ".join(loads(base64.b64decode(chromPeak.comments))),
            )
            chromPeaks.append(chromPeak)

        if self.metabolisationExperiment:
            for tracer in SQLSelectAsObject(
                curs,
                "SELECT t.id AS id, t.name AS name, t.elementCount AS elementCount, t.natural AS isotopeA, t.labelling AS isotopeB, t.deltaMZ AS mzDelta, t.purityN AS enrichmentA, t.purityL AS enrichmentB, t.amountN AS amountA, t.amountL AS amountB, t.monoisotopicRatio AS monoisotopicRatio, t.lowerError AS maxRelNegBias, t.higherError AS maxRelPosBias, t.tracerType AS tracerType FROM tracerConfiguration t",
                newObject=ConfiguredTracer,
            ):
                configTracers[tracer.id] = tracer

        curs.close()
        conn.close()

        csvFile = open(forFile + ".tsv", "w")
        csvFile.write(
            "\t".join(
                [
                    "Num",
                    "MZ",
                    "L_MZ",
                    "D_MZ_Error_ppm",
                    "FoundInScans",
                    "D_MZ_Peak_Error_mean_ppm",
                    "D_MZ_Peak_Error_sd_ppm",
                    "D_MZ_Min_ppm",
                    "D_MZ_Max_ppm",
                    "Quant20_MZ_Diff_ppm",
                    "Quant80_MZ_Diff_ppm",
                    "RT",
                    "Xn",
                    "Charge",
                    "ScanEvent",
                    "Ionisation_Mode",
                    "Tracer",
                    "Area_N",
                    "Area_L",
                    "Abundance_N",
                    "Abundance_L",
                    "Fold",
                    "PeakRatio",
                    "LeftBorder_N",
                    "RightBorder_N",
                    "LeftBorder_L",
                    "RightBorder_L",
                    "Group_ID",
                    "Corr",
                    "ArtificialEICLShift",
                    "Adducts",
                    "FDesc",
                    "Hetero_Elements",
                    "Comments",
                ]
            )
        )
        if len(chromPeaks) > 1:
            i = 1
            for isoRatio in chromPeaks[0].isotopeRatios:
                csvFile.write(
                    "\tObservedIsoRatioMean_%s_%d\tObservedIsoRatioSD_%s_%d\tTheoreticalIsoRatio_%s_%d"
                    % (
                        isoRatio.type,
                        abs(isoRatio.subs),
                        isoRatio.type,
                        abs(isoRatio.subs),
                        isoRatio.type,
                        abs(isoRatio.subs),
                    )
                )

        csvFile.write("\n")

        for chromPeak in sorted(chromPeaks, key=lambda x: x.LPeakCenter):
            hetIso = []
            for hetAtom in chromPeak.heteroIsotopologues:
                pIso = chromPeak.heteroIsotopologues[hetAtom]
                for hetAtomCount in pIso:
                    hetIso.append(
                        "%s%d (scans: %d, obs: %.1f%%, exp: %.1f%%)"
                        % (
                            hetAtom,
                            hetAtomCount,
                            len(pIso[hetAtomCount]),
                            100.0 * mean([d[1] for d in pIso[hetAtomCount]]),
                            100.0 * mean([d[2] for d in pIso[hetAtomCount]]),
                        )
                    )
            if len(hetIso) == 0:
                hetIso = ""
            else:
                hetIso = ", ".join(hetIso)

            addsF = ",".join(chromPeak.adducts)
            if len(addsF) == 0:
                addsF = "-"

            fDesc = ",".join(chromPeak.fDesc)

            scanEvent = ""
            if chromPeak.ionMode == "+":
                scanEvent = self.positiveScanEvent
            elif chromPeak.ionMode == "-":
                scanEvent = self.negativeScanEvent

            mzDelta = 0.0
            if self.metabolisationExperiment:
                mzDelta = configTracers[chromPeak.tracer].mzDelta
            else:
                mzDelta = self.xOffset

            minMZDiffError = -1
            maxMZDiffError = -1
            quantLow20 = -1
            quantHig20 = -1
            try:
                minMZDiffError = min(chromPeak.mzDiffErrors.vals)
                maxMZDiffError = max(chromPeak.mzDiffErrors.vals)
                quantLow20 = np.percentile(chromPeak.mzDiffErrors.vals, 20)
                quantHig20 = np.percentile(chromPeak.mzDiffErrors.vals, 80)
            except:
                pass

            csvFile.write(
                "\t".join(
                    [
                        str(x)
                        for x in [
                            chromPeak.id,
                            chromPeak.mz,
                            chromPeak.lmz,
                            (chromPeak.lmz - chromPeak.mz - chromPeak.tmz) * 1000000.0 / chromPeak.mz,
                            chromPeak.assignedMZs,
                            chromPeak.mzDiffErrors.mean,
                            chromPeak.mzDiffErrors.sd,
                            minMZDiffError,
                            maxMZDiffError,
                            quantLow20,
                            quantHig20,
                            chromPeak.NPeakCenterMin,
                            chromPeak.xCount,
                            chromPeak.loading,
                            scanEvent,
                            chromPeak.ionMode,
                            chromPeak.tracerName,
                            chromPeak.NPeakArea,
                            chromPeak.LPeakArea,
                            chromPeak.NPeakAbundance,
                            chromPeak.LPeakAbundance,
                            chromPeak.NPeakArea / chromPeak.LPeakArea,
                            chromPeak.peaksRatio,
                            chromPeak.NBorderLeft,
                            chromPeak.NBorderRight,
                            chromPeak.LBorderLeft,
                            chromPeak.LBorderRight,
                            chromPeak.fGroupID,
                            chromPeak.peaksCorr,
                            chromPeak.artificialEICLShift,
                            addsF,
                            fDesc,
                            hetIso,
                            chromPeak.comments,
                        ]
                    ]
                )
            )
            for isoRatio in chromPeak.isotopeRatios:
                observedMean = isoRatio.observedRatioMean
                if observedMean is None:
                    observedMean = -1.0
                observedSD = isoRatio.observedRatioSD
                if observedSD is None:
                    observedSD = -1.0
                theoreticalRatio = isoRatio.theoreticalRatio
                if theoreticalRatio is None:
                    theoreticalRatio = -1.0

                csvFile.write("\t%.5f\t%.5f\t%.5f" % (observedMean, observedSD, theoreticalRatio))
            csvFile.write("\n")

        csvFile.write(
            "## MetExtract II %s\n"
            % (
                Bunch(
                    MetExtractVersion=self.meVersion,
                    RVersion=self.rVersion,
                    UUID_ext=self.processingUUID,
                )
                .dumpAsJSon()
                .replace('"', "'")
            )
        )

        processingParams = Bunch()
        processingParams.experimentOperator = self.experimentOperator
        processingParams.experimentID = self.experimentID
        processingParams.experimentComments = self.experimentComments
        processingParams.experimentName = self.experimentName
        csvFile.write("## Experiment parameters %s\n" % (processingParams.dumpAsJSon().replace('"', "'")))

        processingParams = Bunch()
        processingParams.startTime = self.startTime
        processingParams.stopTime = self.stopTime
        processingParams.positiveScanEvent = self.positiveScanEvent
        processingParams.negativeScanEvent = self.negativeScanEvent
        processingParams.metabolisationExperiment = self.metabolisationExperiment
        processingParams.intensityThreshold = self.intensityThreshold
        processingParams.intensityCutoff = self.intensityCutoff
        processingParams.maxLoading = self.maxLoading
        processingParams.xCounts = self.xCountsString
        processingParams.xoffset = self.xOffset
        processingParams.ppm = self.ppm
        processingParams.isotopicPatternCountLeft = self.isotopicPatternCountLeft
        processingParams.isotopicPatternCountRight = self.isotopicPatternCountRight
        processingParams.lowAbundanceIsotopeCutoff = self.lowAbundanceIsotopeCutoff
        processingParams.intensityErrorN = self.intensityErrorN
        processingParams.intensityErrorL = self.intensityErrorL
        processingParams.purityN = self.purityN
        processingParams.purityL = self.purityL

        # 2. Results clustering
        processingParams.minSpectraCount = self.minSpectraCount
        processingParams.clustPPM = self.clustPPM

        # 3. Peak detection
        processingParams.chromPeakPPM = self.chromPeakPPM

        processingParams.eicSmoothingWindow = self.eicSmoothingWindow
        processingParams.eicSmoothingWindowSize = self.eicSmoothingWindowSize
        processingParams.eicSmoothingPolynom = self.eicSmoothingPolynom
        processingParams.scales = self.scales
        processingParams.minCorr = self.minPeakCorr
        processingParams.minCorrelationConvolution = self.minCorrelation
        processingParams.minCorrelationConnections = self.minCorrelationConnections
        csvFile.write("## Data processing parameters %s\n" % (processingParams.dumpAsJSon().replace('"', "'")))

        csvFile.close()

    # plot all detected feature pairs as the first results page in the PDF
    def generateFeaturePairOverviewMap(self, chromPeaks, pdf):
        drawing = Drawing(900, 450)
        sc = ScatterPlot()
        sc.x = 50
        sc.x = 50
        sc.height = 750
        sc.width = 510
        dd = []
        for gi in set([peak.fGroupID for peak in chromPeaks]):
            ddd = []
            for peak in sorted(
                [peak for peak in chromPeaks if peak.fGroupID == gi],
                key=functools.cmp_to_key(lambda x, y: int(100 * (x.NPeakCenter - y.NPeakCenter)) if x.fGroupID == y.fGroupID else int(x.fGroupID - y.fGroupID)),
            ):
                ddd.append((peak.NPeakCenterMin, peak.mz))
            dd.append(ddd)
        sc.data = dd
        sc.yLabel = "M/Z"
        sc.xLabel = "Retention Time [min]"
        sc.lineLabelFormat = noLabel
        drawing.add(sc)
        renderPDF.draw(drawing, pdf, 5, 5)
        pdf.drawString(20, 810, "Plots with the same symbol/color belong to the same group.")
        pdf.drawString(
            20,
            790,
            "Every third group the symbols/colors are repeated but denote a different group.",
        )
        pdf.showPage()

    # include the grouped feature pairs as a) an overlaid EIC plot and b) a list
    # illustrating the relationships between the different feature pairs in a feature group
    def createPDFGroupPages(self, curs, pdf, gPeaks, lGroupID):
        drawing = Drawing(500, 350)
        lp = LinePlot()
        lp.x = 50
        lp.y = 50
        lp.height = 330
        lp.width = 500
        dd = []
        oMinDb = 1000000
        oMaxDb = 0
        for gPeakA in gPeaks:
            gPeak = gPeakA[0]
            eic = gPeakA[1]
            eicN = gPeakA[2]
            times = gPeakA[3]

            minD = int(max(0, gPeak.NPeakCenter - 1 * gPeak.NBorderLeft))
            maxD = int(min(len(eic) - 1, gPeak.NPeakCenter + 1 * gPeak.NBorderRight))
            oMinDb = int(min(oMinDb, int(max(0, gPeak.NPeakCenter - 3 * gPeak.NBorderLeft))))
            oMaxDb = int(
                max(
                    oMaxDb,
                    int(min(len(eic) - 1, gPeak.NPeakCenter + 3 * gPeak.NBorderRight)),
                )
            )

            maxPeakVal = max(eic[int(max(0, minD)) : int(min(len(eic) - 1, maxD))])

            for i in range(0, len(eic)):
                if i < minD or i > maxD:
                    eicN[i] = 0

            if maxPeakVal > 0:
                eic = [e / maxPeakVal for e in eic]
                eicN = [e / maxPeakVal for e in eicN]

            dd.append([(times[i] / 60.0, eic[i]) for i in range(0, len(eic))])
            dd.append([(times[i] / 60.0, eicN[i]) for i in range(0, len(eicN))])
        ddd = []
        for d in dd:
            ddd.append(d[oMinDb:oMaxDb])
        lp.data = ddd
        for i in range(len(ddd)):
            if (i % 2) > 0:
                lp.lines[i].strokeColor = Color(178 / 255.0, 34 / 255.0, 34 / 255.0)
                lp.lines[i].strokeWidth = 0.2
            else:
                lp.lines[i].strokeColor = Color(47 / 255.0, 79 / 255.0, 79 / 255.0)
                lp.lines[i].strokeWidth = 0.1
        lp.joinedLines = 1
        drawing.add(lp)
        renderPDF.draw(drawing, pdf, 5, 5)

        frow = [""]
        data = [frow]
        idCols = {}
        for i in range(len(gPeaks)):
            te = "%.4f/%d" % (gPeaks[i][0].mz, gPeaks[i][0].xCount)
            frow.append(TTR(te))
            arow = [te]
            for j in range(len(gPeaks)):
                arow.append("")
            data.append(arow)
            idCols[gPeaks[i][0].id] = i + 1

        style = [
            ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.white),
            ("BOX", (0, 0), (-1, -1), 0.25, colors.lightslategray),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("LEFTPADDING", (0, 0), (-1, -1), 0),
            ("RIGHTPADDING", (0, 0), (-1, -1), 0),
            ("TOPPADDING", (0, 0), (-1, -1), 0),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ]

        cIds = ",".join(["%d" % f[0].id for f in gPeaks])
        for row in curs.execute("SELECT ff.fID1, ff.fID2, ff.corr FROM featurefeatures ff, chrompeaks c, chrompeaks f WHERE ff.fID1==c.id AND ff.fID2==f.id AND ff.fID1 IN (%s) AND ff.fID2 IN (%s) ORDER BY c.mz" % (cIds, cIds)):
            fI1 = idCols[row[0]]
            fI2 = idCols[row[1]]
            correlation = row[2]
            if fI1 > fI2:
                a = fI2
                fI2 = fI1
                fI1 = a

            bgColor = colors.orange
            if correlation >= self.minCorrelation:
                bgColor = colors.Color(
                    154 / 255.0,
                    205 / 255.0,
                    50 / 255.0,
                    0.5 + 0.4 * (correlation - self.minCorrelation) / (1.0 - self.minCorrelation),
                )  # olivedrab
            elif correlation >= 0:
                bgColor = colors.Color(
                    178 / 255.0,
                    34 / 255.0,
                    34 / 255.0,
                    0.05 + 0.15 * correlation / self.minCorrelation,
                )  # slategrey
            else:
                bgColor = colors.Color(47 / 255.0, 79 / 255.0, 79 / 255.0, 0.2 - 0.5 * correlation)  # firebrick

            style.append(("BACKGROUND", (fI1, fI2), (fI1, fI2), bgColor))
            style.append(("BACKGROUND", (fI2, fI1), (fI2, fI1), bgColor))

        for i in range(len(idCols)):
            style.append(
                (
                    "BACKGROUND",
                    (i + 1, i + 1),
                    (i + 1, i + 1),
                    colors.Color(154 / 255.0, 205 / 255.0, 50 / 255.0, 1),
                )
            )

        style.append(("VALIGN", (0, 0), (-1, 0), "BOTTOM"))
        style.append(("ALIGN", (0, 0), (-1, 0), "CENTER"))
        style.append(("BOTTOMPADDING", (0, 0), (-1, 0), 4))
        style.append(("LEFTPADDING", (0, 0), (-1, 0), 3))

        awidth = 10
        if len(data) > 36:
            style.append(("FONTSIZE", (0, 0), (-1, -1), 5))
            awidth = 5

        widths = [50]
        widths.extend([awidth] * (len(data) - 1))
        table = Table(data, colWidths=widths, rowHeights=widths, style=style)

        table.wrapOn(pdf, 2 * len(data), 2 * len(data))
        table.drawOn(pdf, 60, 400)

        pdf.drawString(
            20,
            810,
            "Feature group %d. %d participating mzs (next page)" % (lGroupID, len(gPeaks)),
        )
        pdf.showPage()

        ptext = ""
        j = 0
        for f in sorted(gPeaks, key=lambda x: x[0].mz):
            gPeak = f[0]
            fDesc = ", ".join(gPeak.fDesc)
            fAdducts = ""
            if len(gPeak.adducts) == 1:
                fAdducts = gPeak.adducts[0]
            elif len(gPeak.adducts) > 1:
                fAdducts = ", ".join(gPeak.adducts)
            ptext = ptext + "\u2022 %.6f (%s<sub>%d</sub>): <u>%s</u>  %s<br/>" % (
                gPeak.mz,
                self.labellingElement,
                gPeak.xCount,
                fAdducts,
                fDesc,
            )
            j = j + 1
        p = Paragraph(ptext, style=getSampleStyleSheet()["Normal"])
        w, h = p.wrap(pagesizes.A4[0] * 0.9, pagesizes.A4[1] * 0.9)
        p.wrapOn(pdf, pagesizes.A4[0] * 0.9, pagesizes.A4[1] * 0.9)
        p.drawOn(pdf, *coord(10 * mm, pagesizes.A4[1] - h - 10 * mm))
        pdf.showPage()

    # plot the MS scan stored in the database for each detected feature pair
    def plotMSScanForFeaturePair(self, peak, tracer, scanEvent, mzxml, pdf):
        try:
            mInt = 0
            drawings = {}
            lps = {}
            iscans = [0]
            for iscan in iscans:
                drawing = Drawing(500, 350)
                drawings[iscan] = drawing

                lp = LinePlot()
                lps[iscan] = lp

                lp.x = 40
                lp.y = 30
                lp.height = 150
                lp.width = 520

                peak.mz + peak.xCount * tracer.mzDelta / peak.loading
                scan = mzxml.getIthMS1Scan(peak.NPeakCenter + iscan, scanEvent)

                dd = []
                colors = []
                strokes = []

                for i in range(len(scan.mz_list)):
                    cmz = scan.mz_list[i]
                    cint = scan.intensity_list[i]
                    dup = peak.mz + peak.xCount * tracer.mzDelta / peak.loading

                    for o in [0, 1, peak.xCount - 1, peak.xCount, peak.xCount + 1]:
                        dup = peak.mz + o * tracer.mzDelta / peak.loading
                        if (dup * (1.0 - self.ppm / 1000000.0)) <= cmz <= (dup * (1.0 + self.ppm / 1000000.0)):
                            dd.append([(cmz, 0), (cmz, cint)])
                            colors.append(Color(178 / 255.0, 34 / 255.0, 34 / 255.0))  # firebrick
                            strokes.append(1.2)

                    if len(peak.heteroIsotopologues) > 0:
                        for pIso in peak.heteroIsotopologues:
                            himz = 0
                            if self.heteroAtoms[pIso].mzOffset < 0:
                                himz = peak.mz + self.heteroAtoms[pIso].mzOffset / peak.loading  # delta m/z is negative
                            else:
                                himz = peak.mz + peak.xCount * tracer.mzDelta / peak.loading + self.heteroAtoms[pIso].mzOffset / peak.loading  # delta m/z is positive

                            if (himz * (1.0 - self.ppm / 1000000.0)) <= cmz <= (himz * (1.0 + self.ppm / 1000000.0)):
                                dd.append([(cmz, 0), (cmz, cint)])
                                colors.append(Color(154 / 255.0, 205 / 255.0, 50 / 255.0))  # olivedrab
                                strokes.append(1.2)

                dd.append([(peak.mz - 5, 0), (peak.mz - 5, 0)])
                colors.append(Color(47 / 255.0, 79 / 255.0, 79 / 255.0))  # slategrey
                strokes.append(0.1)
                for i in range(len(scan.mz_list)):
                    cmz = scan.mz_list[i]
                    cint = scan.intensity_list[i]
                    dup = peak.mz + peak.xCount * tracer.mzDelta / peak.loading

                    if (peak.mz - 5) < cmz < (dup + 5):
                        dd.append([(cmz, 0), (cmz, cint)])
                        colors.append(Color(47 / 255.0, 79 / 255.0, 79 / 255.0))  # slategrey
                        strokes.append(0.1)
                        mInt = max(mInt, cint)

                dd.append(
                    [
                        (peak.mz + peak.xCount * tracer.mzDelta / peak.loading + 5, 0),
                        (peak.mz + peak.xCount * tracer.mzDelta / peak.loading + 5, 0),
                    ]
                )
                colors.append(Color(47 / 255.0, 79 / 255.0, 79 / 255.0))  # slategrey
                strokes.append(0.1)

                assert len(dd) == len(colors) == len(strokes)

                lp.data = dd
                for i in range(len(dd)):
                    lp.lines[i].strokeColor = colors[i]
                    lp.lines[i].strokeWidth = strokes[i]

                if iscan == 0:
                    lp.xValueAxis.labelTextFormat = "%.4f"
                    lp.xValueAxis.valueSteps = [
                        peak.mz,
                        peak.mz + peak.xCount * tracer.mzDelta / peak.loading,
                    ]
                else:
                    lp.yValueAxis.valueSteps = [-1000000]
                    lp.xValueAxis.valueSteps = [-100]
                lp.joinedLines = 1

                drawing.add(lp)
            for iscan in iscans:
                drawing = drawings[iscan]
                lp = lps[iscan]

                lp.yValueAxis.valueMin = 0
                lp.yValueAxis.valueMax = mInt
                lp.xValueAxis.valueMin = peak.mz - 5
                lp.xValueAxis.valueMax = peak.mz + peak.xCount * tracer.mzDelta / peak.loading + 5

                if iscan != 0:
                    lp.xValueAxis.visible = False
                    lp.yValueAxis.visible = False
                renderPDF.draw(drawing, pdf, 15 + iscan * 2, 475 + iscan * 2)

                if iscan == 0:
                    p = Paragraph(
                        "Mass spectrum at %.2f min" % (scan.retention_time / 60.0),
                        style=getSampleStyleSheet()["Normal"],
                    )
                    w, h = p.wrap(pagesizes.A4[0] * 0.9, pagesizes.A4[1] * 0.9)
                    p.wrapOn(pdf, pagesizes.A4[0] * 0.9, pagesizes.A4[1] * 0.9)
                    p.drawOn(pdf, *coord(30, 660))
        except Exception as ex:
            import traceback

            traceback.print_exc()

            self.printMessage("Error in %s: %s" % (self.file, str(ex)), type="error")
            self.postMessageToProgressWrapper("failed", self.pID)

    # plot the EIC stored in the database for each detected feature pair
    def plotXICForFeature(
        self,
        center,
        centerD,
        eic,
        eicIso,
        leftFlank,
        rightFlank,
        leftFlankD,
        rightFlankD,
        times,
        yPosPlot,
        pdf,
    ):
        assert len(eic) == len(eicIso) == len(times)

        drawing = Drawing(500, 350)
        lp = LinePlot()
        lp.x = 40
        lp.y = 30
        lp.height = 200
        lp.width = 300
        minD = int(center - 1 * leftFlank)
        maxD = int(center + 1 * rightFlank)

        dd = [
            [(times[i] / 60.0, eicIso[i]) for i in range(0, len(eic))],
            [(times[i] / 60.0, eic[i]) for i in range(0, len(eic))],
            [(times[i] / 60.0, eic[i]) for i in range(minD, maxD + 1)],
        ]

        lp.data = dd
        lp.lines[1].strokeColor = Color(47 / 255.0, 79 / 255.0, 79 / 255.0)
        lp.lines[2].strokeColor = Color(178 / 255.0, 34 / 255.0, 34 / 255.0)
        lp.lines[0].strokeColor = Color(193 / 255.0, 205 / 255.0, 193 / 255.0)
        lp.lines[0].strokeWidth = 0.1
        lp.lines[1].strokeWidth = 0.1
        lp.lines[2].strokeWidth = 0.1
        lp.joinedLines = 1
        drawing.add(lp)
        renderPDF.draw(drawing, pdf, 15, yPosPlot)
        drawing = Drawing(500, 350)
        lp = LinePlot()
        lp.x = 40
        lp.y = 30
        lp.height = 200
        lp.width = 200
        minDsmall = int(max(0, min(int(centerD - 3 * leftFlankD), int(centerD - 3 * leftFlankD))))
        maxDsmall = int(
            min(
                max(int(centerD + 3 * rightFlankD), int(centerD + 3 * rightFlankD)),
                len(eic) - 1,
            )
        )
        dd = [
            [(times[i] / 60.0, eicIso[i]) for i in range(minDsmall, maxDsmall + 1)],
            [(times[i] / 60.0, eic[i]) for i in range(minDsmall, maxDsmall + 1)],
            [(times[i] / 60.0, eic[i]) for i in range(minD, maxD + 1)],
        ]
        lp.data = dd
        lp.lines[1].strokeColor = Color(47 / 255.0, 79 / 255.0, 79 / 255.0)
        lp.lines[2].strokeColor = Color(178 / 255.0, 34 / 255.0, 34 / 255.0)
        lp.lines[0].strokeColor = Color(193 / 255.0, 205 / 255.0, 193 / 255.0)
        lp.lines[0].strokeWidth = 0.1
        lp.lines[1].strokeWidth = 0.1
        lp.lines[2].strokeWidth = 0.2
        lp.joinedLines = 1
        drawing.add(lp)
        renderPDF.draw(drawing, pdf, 330, yPosPlot)

    # create a PDF file illustrating the detected feature pairs and convoluted feature groups
    def writeResultsToPDF(self, mzxml, reportFunction=None):
        conn = connect(self.file + getDBSuffix(), isolation_level="DEFERRED")
        conn.execute("""PRAGMA synchronous = OFF""")
        conn.execute("""PRAGMA journal_mode = OFF""")
        curs = conn.cursor()

        allChromPeaks = []
        configTracers = {}
        try:
            for chromPeak in SQLSelectAsObject(
                curs,
                "SELECT c.id AS id, g.fGroupID AS fGroupID, c.mz AS mz, c.lmz AS lmz, c.xcount AS xCount, c.Loading AS loading, "
                "c.ionMode AS ionMode, c.NPeakCenter AS NPeakCenter, c.NPeakCenterMin AS NPeakCenterMin, "
                "c.NPeakScale AS NPeakScale, c.NPeakArea AS NPeakArea, c.LPeakCenter AS LPeakCenter, "
                "c.LPeakCenterMin AS LPeakCenterMin, c.LPeakScale AS LPeakScale, c.LPeakArea AS LPeakArea, "
                "c.NBorderLeft as NBorderLeft, c.NBorderRight as NBorderRight, c.LBorderLeft as LBorderLeft, c.LBorderRight as LBorderRight,"
                "c.peaksCorr AS peaksCorr, c.assignedMZs AS assignedMZs, c.heteroAtoms AS HAs, c.adducts AS ADs, "
                "c.fDesc AS DSc, t.id AS tracer, t.name AS tracerName, "
                "c.peaksRatio AS peaksRatio, c.peaksRatioMp1 AS peaksRatioMp1, c.peaksRatioMPm1 AS peaksRatioMPm1, c.comments AS comments FROM "
                "chromPeaks c INNER JOIN featureGroupFeatures g ON c.id=g.fID INNER JOIN tracerConfiguration t ON c.tracer=t.id",
                newObject=ChromPeakPair,
            ):
                chromPeak.NPeakCenterMin = chromPeak.NPeakCenterMin / 60.0
                chromPeak.LPeakCenterMin = chromPeak.LPeakCenterMin / 60.0
                setattr(
                    chromPeak,
                    "heteroIsotopologues",
                    loads(base64.b64decode(chromPeak.HAs)),
                )
                setattr(chromPeak, "adducts", loads(base64.b64decode(chromPeak.ADs)))
                setattr(chromPeak, "fDesc", loads(base64.b64decode(chromPeak.DSc)))
                setattr(chromPeak, "comments", loads(base64.b64decode(chromPeak.comments)))
                chromPeak.assignedMZs = loads(base64.b64decode(chromPeak.assignedMZs))

                allChromPeaks.append(chromPeak)

            if self.metabolisationExperiment:
                for tracer in SQLSelectAsObject(
                    curs,
                    "SELECT t.id AS id, t.name AS name, t.elementCount AS elementCount, t.natural AS isotopeA, t.labelling AS isotopeB, t.deltaMZ AS mzDelta, t.purityN AS enrichmentA, t.purityL AS enrichmentB, t.amountN AS amountA, t.amountL AS amountB, t.monoisotopicRatio AS monoisotopicRatio, t.lowerError AS maxRelNegBias, t.higherError AS maxRelPosBias, t.tracerType AS tracerType FROM tracerConfiguration t",
                    newObject=ConfiguredTracer,
                ):
                    configTracers[tracer.id] = tracer
            else:
                for tracer in SQLSelectAsObject(
                    curs,
                    "SELECT t.id AS id, t.name AS name FROM tracerConfiguration t",
                    newObject=ConfiguredTracer,
                ):
                    configTracers[tracer.id] = tracer

            self.postMessageToProgressWrapper("text", "Creating PDF")
            pdf = canvas.Canvas(self.file + ".pdf", pagesize=pagesizes.A4)
            self.writeSettingsToPDF(pdf)

            for tracerID, tracer in configTracers.items():
                chromPeaks = [c for c in allChromPeaks if c.tracer == tracerID]

                self.writeCurrentTracerToPDF(pdf, tracer)
                self.generateFeaturePairOverviewMap(chromPeaks, pdf)

                lGroupID = -1
                gPeaks = []

                cpeak = 0
                sortedChromPeaks = sorted(
                    chromPeaks,
                    key=functools.cmp_to_key(lambda x, y: int(100 * (x.NPeakCenter - y.NPeakCenter)) if x.fGroupID == y.fGroupID else int(x.fGroupID - y.fGroupID)),
                )

                for si in range(len(sortedChromPeaks)):
                    peak = sortedChromPeaks[si]
                    if reportFunction is not None:
                        reportFunction(
                            1.0 * si / len(sortedChromPeaks),
                            "%d feature pairs remaining" % (len(sortedChromPeaks) - si),
                        )

                    scanEvent = ""
                    if peak.ionMode == "+":
                        scanEvent = self.positiveScanEvent
                    elif peak.ionMode == "-":
                        scanEvent = self.negativeScanEvent

                    self.postMessageToProgressWrapper(
                        "text",
                        "%s: Creating PDF (%d features remaining)" % (tracer.name, len(chromPeaks) - cpeak),
                    )
                    cpeak += 1

                    if lGroupID != -1 and lGroupID != peak.fGroupID and len(gPeaks) > 1 and platform.system() != "Darwin":
                        self.createPDFGroupPages(curs, pdf, gPeaks, lGroupID)

                    if lGroupID != peak.fGroupID:
                        lGroupID = peak.fGroupID
                        gPeaks = []

                    eicN, times, scanIds, mzsN = mzxml.getEIC(peak.mz, self.chromPeakPPM, filterLine=scanEvent)
                    eicN = smoothDataSeries(
                        times,
                        eicN,
                        windowLen=self.eicSmoothingWindowSize,
                        window=self.eicSmoothingWindow,
                        polynom=self.eicSmoothingPolynom,
                    )
                    eicN_cropped = copy(eicN)

                    eicL, times, scanIds, mzsL = mzxml.getEIC(
                        peak.mz + peak.xCount * tracer.mzDelta / peak.loading,
                        self.chromPeakPPM,
                        filterLine=scanEvent,
                    )
                    eicL = smoothDataSeries(
                        times,
                        eicL,
                        windowLen=self.eicSmoothingWindowSize,
                        window=self.eicSmoothingWindow,
                        polynom=self.eicSmoothingPolynom,
                    )

                    eicNP1, times, scanIds, mzsNP1 = mzxml.getEIC(
                        peak.mz + 1 * tracer.mzDelta / peak.loading,
                        self.chromPeakPPM,
                        filterLine=scanEvent,
                    )
                    eicNP1 = smoothDataSeries(
                        times,
                        eicNP1,
                        windowLen=self.eicSmoothingWindowSize,
                        window=self.eicSmoothingWindow,
                        polynom=self.eicSmoothingPolynom,
                    )

                    eicLm1, times, scanIds, mzsLm1 = mzxml.getEIC(
                        peak.mz + (peak.xCount - 1) * tracer.mzDelta / peak.loading,
                        self.chromPeakPPM,
                        filterLine=scanEvent,
                    )
                    eicLm1 = smoothDataSeries(
                        times,
                        eicLm1,
                        windowLen=self.eicSmoothingWindowSize,
                        window=self.eicSmoothingWindow,
                        polynom=self.eicSmoothingPolynom,
                    )

                    gPeaks.append([peak, eicN, eicN_cropped, times, scanIds])

                    if True:  # Plot of Mass spectra for the feature pair
                        self.plotMSScanForFeaturePair(peak, tracer, scanEvent, mzxml, pdf)

                    if True:  # Plots for the non-labelled XIC
                        self.plotXICForFeature(
                            peak.NPeakCenter,
                            peak.NPeakCenter,
                            eicN,
                            eicNP1,
                            peak.NBorderLeft,
                            peak.NBorderRight,
                            peak.NBorderLeft,
                            peak.NBorderRight,
                            times,
                            240,
                            pdf,
                        )

                        p = Paragraph(
                            "XIC of non-labelled isotopologue",
                            style=getSampleStyleSheet()["Normal"],
                        )
                        w, h = p.wrap(pagesizes.A4[0] * 0.9, pagesizes.A4[1] * 0.9)
                        p.wrapOn(pdf, pagesizes.A4[0] * 0.9, pagesizes.A4[1] * 0.9)
                        p.drawOn(pdf, *coord(30, 475))

                    if True:  # Plots for the labelled XIC
                        self.plotXICForFeature(
                            peak.LPeakCenter,
                            peak.NPeakCenter,
                            eicL,
                            eicLm1,
                            peak.LBorderLeft,
                            peak.LBorderRight,
                            peak.LBorderLeft,
                            peak.LBorderRight,
                            times,
                            5,
                            pdf,
                        )

                        p = Paragraph(
                            "XIC of labelled isotopologue",
                            style=getSampleStyleSheet()["Normal"],
                        )
                        w, h = p.wrap(pagesizes.A4[0] * 0.9, pagesizes.A4[1] * 0.9)
                        p.wrapOn(pdf, pagesizes.A4[0] * 0.9, pagesizes.A4[1] * 0.9)
                        p.drawOn(pdf, *coord(30, 240))

                    l = int(
                        max(
                            peak.NPeakCenter - peak.NBorderLeft,
                            peak.LPeakCenter - peak.LBorderLeft,
                        )
                    )
                    h = int(
                        min(
                            peak.NPeakCenter + peak.NBorderRight,
                            peak.LPeakCenter + peak.LBorderRight,
                        )
                    )

                    lines = []

                    lines.append("<u>%s</u> Ion Mode: <u>%s</u>" % (scanEvent, peak.ionMode))
                    lines.append(
                        "Native: m/z: <u>%.6f</u> RT: <u>%.2f</u> Area: <u>%.2f</u> LeftBorder: <u>%.1f</u> RightBorder: <u>%.1f</u>"
                        % (
                            peak.mz,
                            peak.NPeakCenterMin,
                            peak.NPeakArea,
                            peak.NBorderLeft,
                            peak.NBorderRight,
                        )
                    )
                    lines.append(
                        "Labelled: m/z: <u>%.6f</u> (Error [ppm]: <u>%.2f</u>) RT: <u>%.2f</u> Area: <u>%.2f</u> LeftBorder: <u>%.1f</u> RightBorder: <u>%.1f</u>"
                        % (
                            peak.lmz,
                            (peak.lmz - peak.mz - peak.xCount * tracer.mzDelta / peak.loading) * 1000000.0 / peak.mz,
                            peak.LPeakCenterMin,
                            peak.LPeakArea,
                            peak.LBorderLeft,
                            peak.LBorderRight,
                        )
                    )
                    lines.append(
                        "Scanratios M/M': <u>%.4f</u> (Area Ratio <u>%.4f</u>) M+1/M: <u>%.4f</u> M'-1/M': <u>%.4f</u>"
                        % (
                            peak.peaksRatio,
                            peak.NPeakArea / peak.LPeakArea,
                            peak.peaksRatioMp1,
                            peak.peaksRatioMPm1,
                        )
                    )
                    lines.append(
                        "ID: <u>%d</u> %s<sub>n</sub>: <u>%d</u> Z: <u>%d</u> Correlation: <u>%.4f</u> (<u>%d</u> scans) valid SIL patterns: <u>%d</u> scans Feature group: <u>%d</u>"
                        % (
                            peak.id,
                            self.labellingElement,
                            peak.xCount,
                            peak.loading,
                            peak.peaksCorr,
                            h - l + 1,
                            len(peak.assignedMZs),
                            peak.fGroupID,
                        )
                    )
                    if peak.artificialEICLShift != 0:
                        lines.append("Artificial EIC (M') shift: <u>%d</u> scans" % peak.artificialEICLShift)
                    if len(peak.comments) > 0:
                        lines.extend(peak.comments)

                    lines.append("<br/>Putative hetero atoms: ")
                    if len(peak.heteroIsotopologues) > 0:
                        for hetAtom in peak.heteroIsotopologues:
                            pIso = peak.heteroIsotopologues[hetAtom]
                            for hetAtomCount in pIso:
                                isotopeNumber, element = re.match(r"([0-9]+)([a-zA-Z]+)", hetAtom, re.I).groups()
                                isotopeNumber = int(isotopeNumber)
                                lines.append(
                                    "<super>%d</super>%s<sub>%d</sub> (%s m/z: %.4f, scans: %d, obs: %.1f%%, exp: %.1f%%)"
                                    % (
                                        isotopeNumber,
                                        element,
                                        hetAtomCount,
                                        "\u0394",
                                        self.heteroAtoms[hetAtom].mzOffset,
                                        len(pIso[hetAtomCount]),
                                        100.0 * mean([d[1] for d in pIso[hetAtomCount]]),
                                        100.0 * mean([d[2] for d in pIso[hetAtomCount]]),
                                    )
                                )
                    else:
                        lines.append("- (none detected)")

                    p = Paragraph("<br/>".join(lines), style=getSampleStyleSheet()["Normal"])
                    pw, ph = p.wrap(pagesizes.A4[0] * 0.9, pagesizes.A4[1] * 0.9)
                    p.wrapOn(pdf, pagesizes.A4[0] * 0.9, -pagesizes.A4[1] * 0.15)
                    p.drawOn(pdf, *coord(20, 820 - ph))

                    pdf.showPage()

                if lGroupID != -1 and len(gPeaks) > 1 and platform.system() != "Darwin":
                    self.createPDFGroupPages(curs, pdf, gPeaks, lGroupID)

            pdf.save()

            curs.close()
            conn.close()

        except Exception as ex:
            import traceback

            traceback.print_exc()

            self.printMessage("Error in %s: %s" % (self.file, str(ex)), type="error")
            self.postMessageToProgressWrapper("failed")

    # store the detected feature pairs in a new mzXML file. Only those MS peaks will be included, which contribute to
    # the chromatographic peaks of a valid feature pair
    def writeResultsToNewMZXMLIntermediateObject(self, mzxml, newMZXMLData, chromPeaks):
        for peak in chromPeaks:
            # writeMZXML: 0001/1: 12C  0010/2: 12C-Iso  0100/4: 13C-Iso  1000/8: 13C

            if peak.ionMode == "+":
                scanEv = self.positiveScanEvent
            elif peak.ionMode == "-":
                scanEv = self.negativeScanEvent

            if self.writeMZXML & 1:
                scans, times, scanIds = mzxml.getArea(
                    peak.NPeakCenter - peak.NBorderLeft,
                    peak.NPeakCenter + peak.NBorderRight,
                    peak.mz,
                    self.chromPeakPPM,
                    scanEv,
                )

                for w in range(len(scans)):
                    curScan = scans[w]
                    scanid = scanIds[w]
                    if not (scanid in newMZXMLData):
                        newMZXMLData[scanid] = Bunch(mzs=[], ints=[], rt=mzxml.getScanByID(scanid).retention_time)
                    for mz, inte in curScan:
                        newMZXMLData[scanid].mzs.append(mz)
                        newMZXMLData[scanid].ints.append(inte)
            if self.writeMZXML & 2:
                scans, times, scanIds = mzxml.getArea(
                    peak.NPeakCenter - peak.NBorderLeft,
                    peak.NPeakCenter + peak.NBorderRight,
                    peak.mz + self.xOffset / peak.loading,
                    self.chromPeakPPM,
                    scanEv,
                )

                for w in range(len(scans)):
                    curScan = scans[w]
                    scanid = scanIds[w]
                    if not (scanid in newMZXMLData):
                        newMZXMLData[scanid] = Bunch(mzs=[], ints=[], rt=mzxml.getScanByID(scanid).retention_time)
                    for mz, inte in curScan:
                        newMZXMLData[scanid].mzs.append(mz)
                        newMZXMLData[scanid].ints.append(inte)
            if self.writeMZXML & 4:
                scans, times, scanIds = mzxml.getArea(
                    peak.NPeakCenter - peak.NBorderLeft,
                    peak.NPeakCenter + peak.NBorderRight,
                    peak.mz + (peak.xCount - 1) * self.xOffset / peak.loading,
                    self.chromPeakPPM,
                    scanEv,
                )

                for w in range(len(scans)):
                    curScan = scans[w]
                    scanid = scanIds[w]
                    if not (scanid in newMZXMLData):
                        newMZXMLData[scanid] = Bunch(mzs=[], ints=[], rt=mzxml.getScanByID(scanid).retention_time)
                    for mz, inte in curScan:
                        newMZXMLData[scanid].mzs.append(mz)
                        newMZXMLData[scanid].ints.append(inte)
            if self.writeMZXML & 8:
                scans, times, scanIds = mzxml.getArea(
                    peak.NPeakCenter - peak.NBorderLeft,
                    peak.NPeakCenter + peak.NBorderRight,
                    peak.mz + peak.xCount * self.xOffset / peak.loading,
                    self.chromPeakPPM,
                    scanEv,
                )

                for w in range(len(scans)):
                    curScan = scans[w]
                    scanid = scanIds[w]
                    if not (scanid in newMZXMLData):
                        newMZXMLData[scanid] = Bunch(mzs=[], ints=[], rt=mzxml.getScanByID(scanid).retention_time)
                    for mz, inte in curScan:
                        newMZXMLData[scanid].mzs.append(mz)
                        newMZXMLData[scanid].ints.append(inte)

    # helper method if more than one tracer substance is used
    def writeIntermediateMZXMLDataToNewMZXMLFile(self, mzxml, newMZXMLData):
        try:
            self.printMessage("Writing MzXML file..", type="info")
            if ".mzxml" in self.file.lower():
                toFile = self.file[0 : self.file.lower().rfind(".mzxml")] + ".proc.mzXML"
                mzxml.resetMZData(self.file, toFile, newMZXMLData)
            elif ".mzml" in self.file.lower() and False:  ## mzml export currently not supported
                toFile = self.file[0 : self.file.lower().rfind(".mzml")] + "proc.mzML"
                self.printMessage("Only mzXML files can be written", type="error")
            else:
                return RuntimeError("Invalid file. Cannot write new data")

        except Exception as ex:
            self.printMessage("Cannot write MzXML file.. (%s)" % ex, type="error")

    # stores the TICs of the LC-HRMS data in the database
    def writeTICsToDB(self, mzxml, scanEvents):
        conn = connect(self.file + getDBSuffix(), isolation_level="DEFERRED")
        conn.execute("""PRAGMA synchronous = OFF""")
        conn.execute("""PRAGMA journal_mode = OFF""")
        curs = conn.cursor()

        ## save TICs
        ## save mean, sd scan time
        i = 1
        for pol, scanEvent in scanEvents.items():
            if scanEvent != "None":
                TIC, times, scanIds = mzxml.getTIC(filterLine=scanEvent)
                curs.execute(
                    "INSERT INTO tics(id, loading, scanevent, times, intensities) VALUES(%d, '%s', '%s', '%s', '%s')"
                    % (
                        i,
                        pol,
                        scanEvent,
                        ";".join([str(t) for t in times]),
                        ";".join(["%.0f" % t for t in TIC]),
                    )
                )
                i = i + 1

                scanTimes = [times[i + 1] - times[i] for i in range(len(times) - 1)]
                curs.execute("INSERT INTO stats(key, value) VALUES('%s', '%s')" % ("MeanScanTime_%s" % pol, str(mean(scanTimes))))
                curs.execute("INSERT INTO stats(key, value) VALUES('%s', '%s')" % ("SDScanTime_%s" % pol, str(sd(scanTimes))))
                curs.execute("INSERT INTO stats(key, value) VALUES('%s', '%s')" % ("NumberOfScans_%s" % pol, len(times)))
                curs.execute(
                    "INSERT INTO stats(key, value) VALUES('%s', '%s')"
                    % (
                        "TotalNumberOfSignals_%s" % pol,
                        mzxml.getSignalCount(filterLine=scanEvent),
                    )
                )

                minInt, maxInt, avgInt = mzxml.getMinMaxAvgSignalIntensities(filterLine=scanEvent)
                curs.execute("INSERT INTO stats(key, value) VALUES('%s', '%s')" % ("MinSignalInt_%s" % pol, minInt))
                curs.execute("INSERT INTO stats(key, value) VALUES('%s', '%s')" % ("MaxSignalInt_%s" % pol, maxInt))
                curs.execute("INSERT INTO stats(key, value) VALUES('%s', '%s')" % ("AvgSignalInt_%s" % pol, avgInt))

        conn.commit()
        curs.close()
        conn.close()

    # method, which is called by the multiprocessing module to actually process the LC-HRMS data
    def identify(self):
        try:
            start = time.time()

            # region Initialize data processing pipeline
            ######################################################################################

            self.printMessage("File: %s" % self.file, type="info")

            self.postMessageToProgressWrapper("start")
            self.postMessageToProgressWrapper("max", 100.0)
            self.postMessageToProgressWrapper("value", 0.0)
            self.postMessageToProgressWrapper("text", "Initialising")

            if not USEGRADIENTDESCENDPEAKPICKING:
                self.CP = MassSpecWavelet(
                    self.chromPeakFile,
                    scales=self.scales,
                    snrTh=self.snrTh,
                    minScans=self.hAMinScans,
                )
            else:
                from .chromPeakPicking.GradientPeaks import GradientPeaks

                self.CP = GradientPeaks()  ## generic gradient descend peak picking - do not use. Parameters need to be optimized
                self.CP = GradientPeaks(minInt=1000, minIntFlanks=1, minIncreaseRatio=0.01)  ## LTQ Orbitrap XL
                self.CP = GradientPeaks(
                    minInt=10000,
                    minIntFlanks=10,
                    minIncreaseRatio=0.15,
                    expTime=[10, 250],
                )  ## Swiss Orbitrap HF data
                self.CP = GradientPeaks(
                    minInt=1000,
                    minIntFlanks=10,
                    minIncreaseRatio=0.05,
                    minDelta=10000,
                    expTime=[5, 150],
                )  ## Bernhard HSS
                self.CP = GradientPeaks(
                    minInt=1000,
                    minIntFlanks=100,
                    minIncreaseRatio=0.5,
                    minDelta=100,
                    expTime=[5, 150],
                )  ##Lin
                # self.CP=GradientPeaks(minInt=5, minIntFlanks=2, minIncreaseRatio=.05, expTime=[15, 150], minDelta=1, minInflectionDelta=2) ## Roitinger
                # self.CP=GradientPeaks(minInt=10000, minIntFlanks=10, minIncreaseRatio=.05, expTime=[5, 45])       ## RNA

            self.BL = Baseline.Baseline()

            self.curPeakId = 1
            self.curEICId = 1
            self.curMZId = 1
            self.curMZBinId = 1
            self.curFeatureGroupID = 1
            self.curMassSpectrumID = 1
            # endregion

            # region Create results database
            ######################################################################################
            self.postMessageToProgressWrapper("text", "Creating results DB")

            if os.path.exists(self.file + getDBSuffix()) and os.path.isfile(self.file + getDBSuffix()):
                os.remove(self.file + getDBSuffix())
            conn = connect(self.file + getDBSuffix(), isolation_level="DEFERRED")
            conn.execute("""PRAGMA synchronous = OFF""")
            conn.execute("""PRAGMA journal_mode = OFF""")
            curs = conn.cursor()

            self.writeConfigurationToDB(conn, curs)

            curs.close()
            conn.close()
            # endregion

            # region Parse mzXML file
            ######################################################################################

            self.postMessageToProgressWrapper("text", "Parsing chromatogram file")

            mzxml = self.parseMzXMLFile()
            newMZXMLData = {}

            self.writeTICsToDB(mzxml, {"+": self.positiveScanEvent, "-": self.negativeScanEvent})

            # endregion

            # Start calculation
            ######################################################################################

            self.postMessageToProgressWrapper("text", "Starting data processing")
            tracerProgressWidth = 100.0 / 1

            tracerNum = 1

            # region Process configured tracer
            ######################################################################################

            tracer = self.configuredTracer

            curTracerProgress = tracerNum / 1

            tracerID = 0
            if self.metabolisationExperiment:
                if tracer is None:
                    self.printMessage("Metabolisation experiment requires a configured tracer, but none was available in the identify method.", type="error")
                    raise Exception("Metabolisation experiment requires a configured tracer, but none was available in the identify method.")

                tracerID = tracer.id
                self.printMessage("Tracer: %s" % tracer.name, type="info")

                ##################################################################################################
                # Attention: delta mz for one labelling atom is always saved in the member variable self.xOffset #
                ##################################################################################################

                self.xOffset = getIsotopeMass(tracer.isotopeB)[0] - getIsotopeMass(tracer.isotopeA)[0]

            else:
                # Full metabolome labelling experiment
                pass
            # endregion

            # region 1. Find 12C 13C partners in the mz dimension (0-25%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", curTracerProgress + 0 * tracerProgressWidth)
            self.postMessageToProgressWrapper("text", "%s: Extracting signal pairs" % tracer.name)

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper(
                    "value",
                    curTracerProgress + 0 * tracerProgressWidth + 0.25 * curVal * tracerProgressWidth,
                )
                self.postMessageToProgressWrapper("text", "%s: Extracting signal pairs (%s)" % (tracer.name, text))

            mzs, negFound, posFound = self.findSignalPairs(curTracerProgress, mzxml, tracer, reportFunction)
            self.writeSignalPairsToDB(mzs, mzxml, tracerID)

            self.printMessage(
                "%s: Extracting signal pairs done. pos: %d neg: %d mzs (including mismatches)" % (tracer.name, posFound, negFound),
                type="info",
            )
            # endregion

            # region 2. Cluster found mz values according to mz value and number of x atoms (25-35%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", curTracerProgress + 0.25 * tracerProgressWidth)
            self.postMessageToProgressWrapper("text", "%s: Clustering found signal pairs" % tracer.name)

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper(
                    "value",
                    curTracerProgress + 0.25 * tracerProgressWidth + 0.1 * curVal * tracerProgressWidth,
                )
                self.postMessageToProgressWrapper(
                    "text",
                    "%s: Clustering found signal pairs (%s)" % (tracer.name, text),
                )

            mzbins = self.clusterFeaturePairs(mzs, reportFunction)
            self.writeFeaturePairClustersToDB(mzbins)
            mzbins = self.removeImpossibleFeaturePairClusters(mzbins)

            self.printMessage(
                "%s: Clustering found signal pairs done. pos: %d neg: %d mz bins (including mismatches)" % (tracer.name, len(mzbins["+"]), len(mzbins["-"])),
                type="info",
            )
            # endregion

            # region 3. Extract chromatographic peaks (35-65%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", curTracerProgress + 0.35 * tracerProgressWidth)
            self.postMessageToProgressWrapper("text", "%s: Separating feature pairs" % tracer.name)

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper(
                    "value",
                    curTracerProgress + 0.35 * tracerProgressWidth + 0.3 * curVal * tracerProgressWidth,
                )
                self.postMessageToProgressWrapper("text", "%s: Separating feature pairs (%s)" % (tracer.name, text))

            chromPeaks = self.findChromatographicPeaksAndWriteToDB(mzbins, mzxml, tracerID, reportFunction)

            self.printMessage(
                "%s: Separating feature pairs done. pos: %d neg: %d chromatographic peaks (including mismatches)"
                % (
                    tracer.name,
                    len([c for c in chromPeaks if c.ionMode == "+"]),
                    len([c for c in chromPeaks if c.ionMode == "-"]),
                ),
                type="info",
            )
            # endregion

            # region 4. Remove isotopolog feature pairs and other false positive findings (65-70%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", curTracerProgress + 0.65 * tracerProgressWidth)
            self.postMessageToProgressWrapper("text", "%s: Removing false positive feature pairs" % tracer.name)

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper(
                    "value",
                    curTracerProgress + 0.65 * tracerProgressWidth + 0.05 * curVal * tracerProgressWidth,
                )
                self.postMessageToProgressWrapper(
                    "text",
                    "%s: Removing false positive feature pairs (%s)" % (tracer.name, text),
                )

            self.removeFalsePositiveFeaturePairsAndUpdateDB(chromPeaks, reportFunction)

            self.printMessage(
                "%s: Removing false positive feature pairs done. %d chromatographic peaks" % (tracer.name, len(chromPeaks)),
                type="info",
            )
            # endregion

            # region 5. Search for hetero atoms (70-75%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", curTracerProgress + 0.7 * tracerProgressWidth)
            self.postMessageToProgressWrapper("text", "%s: Searching for hetero atoms" % tracer.name)

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper(
                    "value",
                    curTracerProgress + 0.7 * tracerProgressWidth + 0.05 * curVal * tracerProgressWidth,
                )
                self.postMessageToProgressWrapper("text", "%s: Annotating feature pairs (%s)" % (tracer.name, text))

            self.annotateFeaturePairs(chromPeaks, mzxml, tracer, reportFunction)
            # endregion

            # region 6. Group feature pairs untargeted using chromatographic peak shape into feature groups (75-80%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", curTracerProgress + 0.75 * tracerProgressWidth)
            self.postMessageToProgressWrapper("text", "%s: Grouping feature pairs" % tracer.name)

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper(
                    "value",
                    curTracerProgress + 0.75 * tracerProgressWidth + 0.05 * curVal * tracerProgressWidth,
                )
                self.postMessageToProgressWrapper("text", "%s: Grouping feature pairs (%s)" % (tracer.name, text))

            self.groupFeaturePairsUntargetedAndWriteToDB(chromPeaks, mzxml, tracer, tracerID, reportFunction)
            # endregion

            # region 7. Extract mass spectra for feature pairs and feature groups (80-95%)
            ######################################################################################

            self.postMessageToProgressWrapper("value", curTracerProgress + 0.8 * tracerProgressWidth)
            self.postMessageToProgressWrapper("text", "%s: Extracting mass spectra to DB" % tracer.name)

            def reportFunction(curVal, text):
                self.postMessageToProgressWrapper(
                    "value",
                    curTracerProgress + 0.8 * tracerProgressWidth + 0.15 * curVal * tracerProgressWidth,
                )
                self.postMessageToProgressWrapper(
                    "text",
                    "%s: Extracting mass spectra to DB (%s)" % (tracer.name, text),
                )

            self.writeMassSpectraToDB(chromPeaks, mzxml, reportFunction)

            # Log time used for processing of individual files
            elapsed = (time.time() - start) / 60.0
            hours = ""
            if elapsed >= 60.0:
                if elapsed < 120.0:
                    hours = "1 hour "
                else:
                    hours = "%d hours " % (elapsed // 60)
            mins = "%.2f min(s)" % (elapsed % 60.0)

            self.printMessage(
                "%s: Calculations finished (%s%s).." % (tracer.name, hours, mins),
                type="info",
            )
            # endregion

            # region 8. Write results to files (95-100%, without progress indicator)
            ######################################################################################

            # W.1 Save results to new MzXML file (intermediate step) (95-100%)

            if self.writeMZXML:
                self.postMessageToProgressWrapper("value", curTracerProgress + 0.95 * tracerProgressWidth)
                self.postMessageToProgressWrapper("text", "%s: Writing results to mzXML.." % tracer.name)

                self.writeResultsToNewMZXMLIntermediateObject(mzxml, newMZXMLData, chromPeaks)
            # endregion

            self.postMessageToProgressWrapper("text", "Writing results to mzXML..")

            # region W.2 Write results to TSV File
            ##########################################################################################
            if self.writeTSV:
                self.postMessageToProgressWrapper("text", "Writing results to TSV..")

                self.writeResultsToTSVFile(self.file)
            # endregion

            # region W.2 Write results to TSV File
            ##########################################################################################
            if self.writeFeatureML:
                self.postMessageToProgressWrapper("text", "Writing results to featureML..")

                self.writeResultsToFeatureML(self.file)
            # endregion

            # region W.3 Write resutls to PDF
            ##########################################################################################
            if self.writePDF:
                self.postMessageToProgressWrapper("text", "Writing results to PDF")

                def reportFunction(curVal, text):
                    self.postMessageToProgressWrapper("value", 90.0 + 5.0 * curVal)
                    self.postMessageToProgressWrapper("text", "Writing results to PDF (%s)" % text)

                self.writeResultsToPDF(mzxml, reportFunction=reportFunction)
            # endregion

            # region W.4 Save all results to new MzXML file
            ##########################################################################################
            if self.writeMZXML:
                self.postMessageToProgressWrapper("text", "Writing results to new mzXML file")

                self.writeIntermediateMZXMLDataToNewMZXMLFile(mzxml, newMZXMLData)
            # endregion

            mzxml.freeMe()
            self.printMessage("%s done.." % self.file, type="info")
            self.postMessageToProgressWrapper("end")

        except Exception as ex:
            import traceback

            traceback.print_exc()

            self.printMessage("Error in %s: %s" % (self.file, str(ex)), type="error")
            self.postMessageToProgressWrapper("failed")
