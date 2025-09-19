from __future__ import print_function, division, absolute_import
import logging

from sqlite3 import *
from operator import itemgetter
import os
from math import floor
import ast
from collections import defaultdict

from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF

from XICAlignment import XICAlignment

from utils import (
    ChromPeakPair,
    getSubGraphs,
    Bunch,
    SQLInsert,
    natSort,
    get_main_dir,
    getSubGraphsFromDictDict,
    CallBackMethod,
)
from runIdentification import ChromPeakPair, getDBSuffix
from TableUtils import TableUtils
from MZHCA import HierarchicalClustering, cutTreeSized

import base64
from pickle import dumps, loads

from math import isnan


import time

import HCA_general
from Chromatogram import Chromatogram
from utils import mean, sd, corr


import exportAsFeatureML


# HELPER METHOD for writing first page of PDF (unused)
def _writeFirstPage(pdf, groupSizePPM, maxTimeDeviation, align, nPolynom):
    curHeight = 780

    pdf.drawString(60, curHeight, "Max. Group Size")
    pdf.drawString(200, curHeight, "%.1f" % groupSizePPM)
    curHeight -= 20

    pdf.drawString(60, curHeight, "Max. Time Deviation")
    pdf.drawString(200, curHeight, "%.2f seconds" % maxTimeDeviation)
    curHeight -= 20

    pdf.drawString(60, curHeight, "Align")
    pdf.drawString(200, curHeight, "%s" % ("yes" if align else "no"))
    curHeight -= 20

    if align:
        pdf.drawString(60, curHeight, "Polynom Order")
        pdf.drawString(200, curHeight, "%d" % nPolynom)
        curHeight -= 20

    pdf.showPage()


# store used configuration to DB file
def writeConfigToDB(
    curs,
    align,
    file,
    groupSizePPM,
    maxLoading,
    maxTimeDeviation,
    xCounts,
    nPolynom,
    negativeScanEvent,
    positiveScanEvent,
    rVersion,
    meVersion,
):
    SQLInsert(curs, "config", key="MEVersion", value=str(meVersion))
    SQLInsert(curs, "config", key="RVersion", value=str(rVersion))

    SQLInsert(curs, "config", key="FPBRACK_xCounts", value=str(xCounts))
    SQLInsert(curs, "config", key="FPBRACK_groupSizePPM", value=str(groupSizePPM))
    SQLInsert(
        curs, "config", key="FPBRACK_positiveScanEvent", value=str(positiveScanEvent)
    )
    SQLInsert(
        curs, "config", key="FPBRACK_negativeScanEvent", value=str(negativeScanEvent)
    )
    SQLInsert(
        curs, "config", key="FPBRACK_maxTimeDeviation", value=str(maxTimeDeviation)
    )
    SQLInsert(curs, "config", key="FPBRACK_maxLoading", value=str(maxLoading))
    SQLInsert(curs, "config", key="FPBRACK_file", value=str(file))
    SQLInsert(curs, "config", key="FPBRACK_align", value=str(align))
    SQLInsert(curs, "config", key="FPBRACK_nPolynom", value=str(nPolynom))


# bracket results
def bracketResults(
    indGroups,
    xCounts,
    groupSizePPM,
    positiveScanEvent=None,
    negativeScanEvent=None,
    maxTimeDeviation=0.36 * 60,
    maxLoading=1,
    file="./results.tsv",
    align=True,
    nPolynom=1,
    expPeakArea=True,
    expApexIntensity=True,
    expPeakSNR=False,
    pwMaxSet=None,
    pwValSet=None,
    pwTextSet=None,
    rVersion="",
    meVersion="",
    generalProcessingParams=Bunch(),
    start=0,
):
    # create results SQLite tables
    resDB = Bunch(conn=None, curs=None)
    if os.path.exists(file + getDBSuffix()) and os.path.isfile(file + getDBSuffix()):
        # Try to remove the file, with retry in case it's still locked
        import time

        for attempt in range(3):
            try:
                os.remove(file + getDBSuffix())
                break
            except PermissionError:
                if attempt < 2:
                    time.sleep(0.1)  # Wait 100ms before retry
                    continue
                else:
                    print(
                        f"Warning: Could not remove existing database file {file + getDBSuffix()}, proceeding anyway"
                    )
                    break
    resDB.conn = connect(file + getDBSuffix())
    # conn.execute('''PRAGMA synchronous = OFF''')
    # conn.execute('''PRAGMA journal_mode = OFF''')
    resDB.curs = resDB.conn.cursor()

    try:
        writePDF = False
        # used for debug purposes
        colos = [colors.red]

        cpf = (
            get_main_dir() + "/XICAlignment.r"
        )  # initialise chromatographic alignment script
        xicAlign = XICAlignment(cpf)

        resDB.curs.execute("DROP TABLE IF EXISTS GroupResults")
        resDB.curs.execute(
            "CREATE TABLE GroupResults(id INTEGER PRIMARY KEY, mz FLOAT, lmz FLOAT, dmz FLOAT, rt FLOAT, xn INTEGER, charge INTEGER, scanEvent TEXT, ionisationMode TEXT, tracer TEXT, OGroup INTEGER, Ion TEXT, LOSS TEXT, M TEXT)"
        )

        resDB.curs.execute("DROP TABLE IF EXISTS FoundFeaturePairs")
        resDB.curs.execute(
            "CREATE TABLE FoundFeaturePairs(resID INTEGER, file TEXT, featurePairID INTEGER, featureGroupID INTEGER, NPeakScale FLOAT, LPeakScale FLOAT, NBorderLeft FLOAT, NBorderRight FLOAT, LBorderLeft FLOAT, LBorderRight FLOAT, areaN FLOAT, areaL FLOAT, featureType TEXT)"
        )

        resDB.curs.execute("DROP TABLE IF EXISTS FileMapping")
        resDB.curs.execute(
            "CREATE TABLE FileMapping(fileName TEXT, filePath TEXT, groupID INTEGER)"
        )

        resDB.curs.execute("DROP TABLE IF EXISTS FileGroups")
        resDB.curs.execute(
            "CREATE TABLE FileGroups(groupName TEXT, id INTEGER PRIMARY KEY)"
        )

        resDB.curs.execute("DROP TABLE IF EXISTS config")
        resDB.curs.execute(
            "CREATE TABLE config(id INTEGER PRIMARY KEY AUTOINCREMENT, key TEXT, value TEXT)"
        )

        writeConfigToDB(
            resDB.curs,
            align,
            file,
            groupSizePPM,
            maxLoading,
            maxTimeDeviation,
            xCounts,
            nPolynom,
            negativeScanEvent,
            positiveScanEvent,
            rVersion,
            meVersion,
        )

        results = []

        # generate results for each input file
        i = 1
        for key in natSort(indGroups.keys()):
            files = indGroups[key]
            SQLInsert(resDB.curs, "FileGroups", groupName=key, id=i)

            for ident in files:
                conn = None
                curs = None
                if os.path.isfile(ident + getDBSuffix()):
                    conn = connect(ident + getDBSuffix(), isolation_level="DEFERRED")
                    conn.execute("""PRAGMA synchronous = OFF""")
                    conn.execute("""PRAGMA journal_mode = OFF""")
                    curs = conn.cursor()
                fname = ident
                if ".mzxml" in ident.lower():
                    fname = ident[
                        max(
                            ident.rfind("/") + 1, ident.rfind("\\") + 1
                        ) : ident.lower().rfind(".mzxml")
                    ]
                if ".mzml" in ident.lower():
                    fname = ident[
                        max(
                            ident.rfind("/") + 1, ident.rfind("\\") + 1
                        ) : ident.lower().rfind(".mzml")
                    ]
                b = Bunch(
                    filePath=ident,
                    fileName=fname,
                    featurePairs=[],
                    conn=conn,
                    curs=curs,
                )

                results.append(b)

                SQLInsert(
                    resDB.curs,
                    "FileMapping",
                    fileName=b.fileName,
                    filePath=b.filePath,
                    groupID=i,
                )

            i += 1

        with open(file, "w") as f:
            # initialise TSV results file

            f.write(
                "Num\tComment\tMZ\tL_MZ\tD_MZ\tMZ_Range\tRT\tRT_Range\tPeakScalesNL\tXn\tCharge\tScanEvent\tIonisation_Mode\tTracer\tOGroup\tIon\tLoss\tM"
            )
            for res in results:
                fname = res.fileName

                if ".mzxml" in fname.lower():
                    fname = fname[
                        max(
                            fname.rfind("/") + 1, fname.rfind("\\") + 1
                        ) : fname.lower().rfind(".mzxml")
                    ]
                if ".mzml" in fname.lower():
                    fname = fname[
                        max(
                            fname.rfind("/") + 1, fname.rfind("\\") + 1
                        ) : fname.lower().rfind(".mzml")
                    ]

                if expPeakArea:
                    f.write("\t" + res.fileName + "_Area_N")
                    f.write("\t" + res.fileName + "_Area_L")
                if expApexIntensity:
                    f.write("\t" + res.fileName + "_Abundance_N")
                    f.write("\t" + res.fileName + "_Abundance_L")
                if expPeakSNR:
                    f.write("\t" + res.fileName + "_SNR_N")
                    f.write("\t" + res.fileName + "_SNR_L")
                if True:
                    f.write("\t" + res.fileName + "_FID")
                    f.write("\t" + res.fileName + "_GroupID")

            f.write("\tdoublePeak")
            f.write("\n")

            # get ionisation modes scan events
            ionModes = {}
            if positiveScanEvent != "None":
                ionModes["+"] = positiveScanEvent
            if negativeScanEvent != "None":
                ionModes["-"] = negativeScanEvent
            assert len(ionModes) > 0

            totalChromPeaks = 0
            totalChromPeaksProc = 0
            tracersDeltaMZ = {}

            ## TODO this used the last open file. May be invalid. improve
            isMetabolisationExperiment = False
            for res in results:
                if res.curs is not None:
                    for row in res.curs.execute(
                        "SELECT key, value FROM config WHERE key='metabolisationExperiment'"
                    ):
                        isMetabolisationExperiment = str(row[1]).lower() == "true"
                        break

            ## TODO this used the last open file. May be invalid. improve
            if isMetabolisationExperiment:
                for res in results:
                    if res.curs is not None:
                        for row in res.curs.execute(
                            "SELECT id, name, deltaMZ FROM tracerConfiguration"
                        ):
                            tracerName = str(row[1])
                            dmz = float(row[2])

                            if tracerName not in tracersDeltaMZ:
                                tracersDeltaMZ[tracerName] = dmz

                            if abs(tracersDeltaMZ[tracerName] - dmz) < 0.00001:
                                logging.warning(
                                    "Warning: Tracers have not been configured identical in all measurement files"
                                )
            else:
                ## TODO this used the last open file. May be invalid. improve
                for res in results:
                    if res.curs is not None:
                        for row in res.curs.execute(
                            "SELECT key, value FROM config WHERE key='xOffset'"
                        ):
                            dmz = float(row[1])

                            if "FLE" not in tracersDeltaMZ:
                                tracersDeltaMZ["FLE"] = dmz

                            if abs(tracersDeltaMZ["FLE"] - dmz) < 0.00001:
                                logging.warning(
                                    "Warning: Tracers have not been configured identical in all measurement files"
                                )

            grp = 1
            if writePDF:
                pdf = canvas.Canvas(file + "_0.pdf")
            if writePDF:
                _writeFirstPage(pdf, groupSizePPM, maxTimeDeviation, align, nPolynom)
            xy = 1

            curNum = 1

            xCounts = set()
            for res in results:
                if res.curs is not None:
                    for row in res.curs.execute(
                        "SELECT DISTINCT xcount FROM chromPeaks"
                    ):
                        xCounts.add(str(row[0]))
            xCounts = sorted(list(xCounts))

            # xCountshelp = []
            # a=xCounts.replace(" ", "").replace(";", ",").split(",")
            # for j in a:
            #     if "-" in j:
            #         xCountshelp.extend(range(int(j[0:j.find("-")]), int(j[j.find("-")+1:len(j)])+1))
            #     else:
            #         xCountshelp.append(int(j))
            # xCounts=sorted(list(set(xCountshelp)))

            # prefetch number of bracketing steps for the progress bar
            totalThingsToDo = 0
            for ionMode in ionModes:
                for tracer in tracersDeltaMZ:
                    for xCount in xCounts:
                        for cLoading in range(maxLoading, 0, -1):
                            totalThingsToDo += 1

            doneSoFar = 0
            doneSoFarPercent = 0
            if pwMaxSet is not None:
                pwMaxSet.put(Bunch(mes="max", val=totalThingsToDo))

            # bracket data
            # process ionModes, tracers, xCount and loading separately
            for ionMode in ionModes:
                scanEvent = ionModes[ionMode]
                for xCount in xCounts:
                    for cLoading in range(maxLoading, 0, -1):
                        # check each processed LC-HRMS file if the used data processing parameters match
                        for res in results:
                            res.featurePairs = []

                            if pwTextSet is not None:
                                # Log time used for bracketing
                                elapsed = (time.time() - start) / 60.0
                                hours = ""
                                if elapsed >= 60.0:
                                    hours = "%d hours " % (elapsed // 60)
                                mins = "%.2f mins" % (elapsed % 60.0)
                                pwTextSet.put(
                                    Bunch(
                                        mes="text",
                                        val="<p align='right' >%s%s elapsed</p>\n\n\n\nProcessing \n  (Ionmode: %s, XCount: %s, Charge: %d) \n  File: %s"
                                        % (
                                            hours,
                                            mins,
                                            ionMode,
                                            xCount,
                                            cLoading,
                                            res.fileName,
                                        ),
                                    )
                                )
                            if res.curs is not None:
                                for row in res.curs.execute(
                                    "SELECT c.id, c.tracer, c.NPeakCenter, c.NPeakCenterMin, c.NPeakScale, c.NSNR, c.NPeakArea, c.mz, c.lmz, c.tmz, c.xcount, c.LPeakCenter, c.LPeakCenterMin, c.LPeakScale, c.LSNR, c.LPeakArea, "
                                    "c.Loading, (SELECT f.fGroupId FROM featureGroupFeatures f WHERE f.fID=c.id) AS fGroupId, (SELECT t.name FROM tracerConfiguration t WHERE t.id=c.tracer) AS tracerName, c.ionMode, "
                                    "c.NPeakAbundance, c.LPeakAbundance, c.NBorderLeft, c.NBorderRight, c.LBorderLeft, c.LBorderRight "
                                    "FROM chromPeaks c WHERE c.ionMode='%s' AND c.xcount='%s' and c.Loading=%d"
                                    % (ionMode, xCount, cLoading)
                                ):
                                    try:
                                        cp = ChromPeakPair(
                                            id=row[0],
                                            tracer=row[1],
                                            NPeakCenter=row[2],
                                            NPeakCenterMin=row[3],
                                            NPeakScale=row[4],
                                            NSNR=row[5],
                                            NPeakArea=row[6],
                                            mz=row[7],
                                            lmz=row[8],
                                            tmz=row[9],
                                            xCount=str(row[10]),
                                            LPeakCenter=row[11],
                                            LPeakCenterMin=row[12],
                                            LPeakScale=row[13],
                                            LSNR=row[14],
                                            LPeakArea=row[15],
                                            loading=row[16],
                                            fGroupID=row[17],
                                            tracerName=row[18],
                                            ionMode=str(row[19]),
                                            NPeakAbundance=float(row[20]),
                                            LPeakAbundance=float(row[21]),
                                            NBorderLeft=float(row[22]),
                                            NBorderRight=float(row[23]),
                                            LBorderLeft=float(row[24]),
                                            LBorderRight=float(row[25]),
                                        )

                                        assert cp.ionMode in ionModes.keys()
                                        res.featurePairs.append(cp)
                                    except TypeError as err:
                                        print(
                                            "  TypeError in file %s, id %s, (Ionmode: %s, XCount: %s, Charge: %d)"
                                            % (
                                                str(res.fileName),
                                                str(row[0]),
                                                ionMode,
                                                xCount,
                                                cLoading,
                                            ),
                                            err.message,
                                        )
                                    except:
                                        print(
                                            "  some general error in file %s, id %s, (Ionmode: %s, XCount: %s, Charge: %d)"
                                            % (
                                                str(res.fileName),
                                                str(row[0]),
                                                ionMode,
                                                xCount,
                                                cLoading,
                                            )
                                        )

                            totalChromPeaks = totalChromPeaks + len(res.featurePairs)

                        if pwTextSet is not None:
                            # Log time used for bracketing
                            elapsed = (time.time() - start) / 60.0
                            hours = ""
                            if elapsed >= 60.0:
                                hours = "%d hours " % (elapsed // 60)
                            mins = "%.2f mins" % (elapsed % 60.0)
                            pwTextSet.put(
                                Bunch(
                                    mes="text",
                                    val="<p align='right' >%s%s elapsed</p>\n\n\n\nClustering \n  (Ionmode: %s, XCount: %s, Charge: %d)"
                                    % (hours, mins, ionMode, xCount, cLoading),
                                )
                            )

                            totalChromPeaks = totalChromPeaks + len(res.featurePairs)

                        for tracer in tracersDeltaMZ:
                            ## TODO this does not work correctly. improve
                            if pwValSet is not None:
                                pwValSet.put(Bunch(mes="value", val=doneSoFar))
                            doneSoFar += 1

                            # get all results that match current bracketing criteria
                            chromPeaksAct = 0
                            allmz = []
                            for res in results:
                                chromPeaksAct = chromPeaksAct + len(
                                    [
                                        chromPeak.mz
                                        for chromPeak in res.featurePairs
                                        if chromPeak.xCount == xCount
                                        and chromPeak.loading == cLoading
                                        and chromPeak.tracerName == tracer
                                        and chromPeak.ionMode == ionMode
                                    ]
                                )

                            allMZs = []
                            for res in results:
                                allMZs.extend(
                                    [
                                        chromPeak.mz
                                        for chromPeak in res.featurePairs
                                        if chromPeak.xCount == xCount
                                        and chromPeak.loading == cLoading
                                        and chromPeak.tracerName == tracer
                                        and chromPeak.ionMode == ionMode
                                    ]
                                )

                            # cluster all current results with HCA
                            allMZs = sorted(allMZs)
                            lastMZ = None
                            u = 0
                            curAllMZs = []

                            ## pre-separate detected feature pairs to improve speed of HCA
                            while len(allMZs) > 0:
                                procCurSet = False
                                if u < len(allMZs) and (
                                    lastMZ is None
                                    or (allMZs[u] - lastMZ)
                                    < 3 * (groupSizePPM * allMZs[u] / 1e6)
                                ):
                                    curAllMZs.append(allMZs[u])
                                    lastMZ = allMZs[u]
                                    u += 1
                                else:
                                    lastMZ = None
                                    allMZs = allMZs[u:]
                                    u = 0
                                    procCurSet = True

                                if procCurSet or len(allMZs) == u:
                                    if len(curAllMZs) > 0:
                                        hc = HierarchicalClustering(
                                            curAllMZs,
                                            dist=lambda x, y: x.getValue()
                                            - y.getValue(),
                                            val=lambda x: x,
                                            mean=lambda x, y: x / y,
                                            add=lambda x, y: x + y,
                                        )

                                        # cut HCA tree in subclusters
                                        for n in cutTreeSized(
                                            hc.getTree(), groupSizePPM
                                        ):
                                            try:
                                                lowMZ = min(
                                                    [
                                                        kid.getValue()
                                                        for kid in n.getKids()
                                                    ]
                                                )
                                                minMZ = max(
                                                    [
                                                        kid.getValue()
                                                        for kid in n.getKids()
                                                    ]
                                                )

                                                maxMZInGroup = lowMZ

                                                partChromPeaks = {}
                                                partXICs = {}
                                                toDel = {}
                                                allChromPeaks = []
                                                for res in results:
                                                    chromPeaks = []
                                                    toDel[res.filePath] = set()
                                                    for i in range(
                                                        len(res.featurePairs)
                                                    ):
                                                        chromPeak = res.featurePairs[i]
                                                        if (
                                                            chromPeak.xCount == xCount
                                                            and chromPeak.loading
                                                            == cLoading
                                                            and chromPeak.tracerName
                                                            == tracer
                                                            and chromPeak.ionMode
                                                            == ionMode
                                                            and chromPeak.mz <= minMZ
                                                        ):
                                                            toDel[res.filePath].add(i)
                                                            chromPeaks.append(chromPeak)
                                                            allChromPeaks.append(
                                                                chromPeak
                                                            )
                                                            maxMZInGroup = max(
                                                                maxMZInGroup,
                                                                chromPeak.mz,
                                                            )

                                                    if align:
                                                        xics = []
                                                        # get EICs of current subCluster
                                                        for chromPeak in chromPeaks:
                                                            res.curs.execute(
                                                                "SELECT x.xicL, x.times, x.scanCount, eicid FROM XICs x, chromPeaks c WHERE c.id=%d AND c.eicID=x.id"
                                                                % chromPeak.id
                                                            )
                                                            i = 0
                                                            for row in res.curs:
                                                                scanCount = float(
                                                                    row[2]
                                                                )
                                                                times = [
                                                                    float(r)
                                                                    for r in row[
                                                                        1
                                                                    ].split(";")
                                                                ]
                                                                xicL = [
                                                                    float(r)
                                                                    for r in row[
                                                                        0
                                                                    ].split(";")
                                                                ]
                                                                eicID = int(row[3])
                                                                i = i + 1
                                                            assert i == 1
                                                            xics.append(
                                                                [
                                                                    scanCount,
                                                                    times,
                                                                    xicL,
                                                                    eicID,
                                                                ]
                                                            )
                                                    partChromPeaks[res.filePath] = (
                                                        chromPeaks
                                                    )
                                                    if align:
                                                        partXICs[res.filePath] = xics

                                                xy = xy + 1
                                                if (xy % 50) == 0:
                                                    if writePDF:
                                                        pdf.save()
                                                        pdf = canvas.Canvas(
                                                            file + "_%d.pdf" % xy
                                                        )
                                                        _writeFirstPage(
                                                            pdf,
                                                            groupSizePPM,
                                                            maxTimeDeviation,
                                                            align,
                                                            nPolynom,
                                                        )
                                                        print(
                                                            "\r   new pdf %d metabolic features.."
                                                            % xy,
                                                        )

                                                # debug purposes: create PDF for current subcluster
                                                if writePDF:
                                                    drawing = Drawing(500, 350)

                                                    lp = LinePlot()
                                                    lp.x = 50
                                                    lp.y = 30
                                                    lp.height = 350
                                                    lp.width = 500

                                                    dd = []

                                                    tr = 0
                                                    pdf.drawString(
                                                        40,
                                                        810,
                                                        "Tracer: %s IonMode: %s"
                                                        % (tracer, ionMode),
                                                    )
                                                    pdf.drawString(
                                                        40,
                                                        795,
                                                        "MZ: (min: %f - max: %f (tolerated: %f ppm: %.1f)) "
                                                        % (
                                                            lowMZ,
                                                            maxMZInGroup,
                                                            minMZ,
                                                            (maxMZInGroup - lowMZ)
                                                            * 1000000
                                                            / lowMZ,
                                                        ),
                                                    )
                                                    pdf.drawString(
                                                        40,
                                                        780,
                                                        "xCount: %s Z: %d Feature pairs: %d"
                                                        % (
                                                            xCount,
                                                            cLoading,
                                                            len(partXICs),
                                                        ),
                                                    )
                                                    pdj = []
                                                    pdjm = []
                                                    for r in partXICs.keys():
                                                        xics = partXICs[r]
                                                        for xic in xics:
                                                            pdj.append(xic[2])
                                                            pdjm.append(xic[1])
                                                            dd.append(
                                                                [
                                                                    (
                                                                        xic[1][i]
                                                                        / 60.0,
                                                                        xic[2][i],
                                                                    )
                                                                    for i in range(
                                                                        0, len(xic[1])
                                                                    )
                                                                ]
                                                            )
                                                            tr = tr + 1
                                                    pdf.drawString(
                                                        40,
                                                        765,
                                                        "Aligning %d EICs" % (tr),
                                                    )

                                                    lp.data = dd

                                                    i = 0
                                                    for k in partXICs.keys():
                                                        lp.lines[i].strokeColor = colos[
                                                            i % len(colos)
                                                        ]
                                                        lp.lines[i].strokeWidth = 0.01
                                                        i = i + 1

                                                    lp.joinedLines = 1
                                                    drawing.add(lp)
                                                    renderPDF.draw(
                                                        drawing, pdf, 10, 380
                                                    )

                                                    alignedEICs = (
                                                        xicAlign.getAligendXICs(
                                                            pdj,
                                                            pdjm,
                                                            align=align,
                                                            nPolynom=nPolynom,
                                                        )[0]
                                                    )
                                                    drawing = Drawing(500, 350)
                                                    lp = LinePlot()
                                                    lp.x = 50
                                                    lp.y = 30
                                                    lp.height = 350
                                                    lp.width = 500

                                                    dd = []

                                                    tr = 0
                                                    for xic in alignedEICs:
                                                        dd.append(
                                                            [
                                                                (i / 60.0, xic[i])
                                                                for i in range(
                                                                    0, len(xic)
                                                                )
                                                            ]
                                                        )
                                                        tr = tr + 1

                                                    lp.data = dd

                                                    i = 0
                                                    for k in range(len(alignedEICs)):
                                                        lp.lines[i].strokeColor = colos[
                                                            i % len(colos)
                                                        ]
                                                        lp.lines[i].strokeWidth = 0.01
                                                        i = i + 1

                                                    lp.joinedLines = 1
                                                    drawing.add(lp)
                                                    renderPDF.draw(drawing, pdf, 10, 0)

                                                    pdf.showPage()

                                                eics = []
                                                peaks = []
                                                scantimes = []
                                                dd = []
                                                do = []
                                                for k in partChromPeaks.keys():
                                                    for i in range(
                                                        len(partChromPeaks[k])
                                                    ):
                                                        do.append(k)
                                                        dd.append(
                                                            partChromPeaks[k][i].mz
                                                        )
                                                        if align:
                                                            eics.append(
                                                                partXICs[k][i][2]
                                                            )
                                                            scantimes.append(
                                                                partXICs[k][i][1]
                                                            )
                                                        peaks.append(
                                                            [partChromPeaks[k][i]]
                                                        )

                                                # optional: align EICs; get bracketed chromatographic peaks
                                                aligned = xicAlign.alignXIC(
                                                    eics,
                                                    peaks,
                                                    scantimes,
                                                    align=align,
                                                    maxTimeDiff=maxTimeDeviation,
                                                    nPolynom=nPolynom,
                                                )
                                                aligned = [
                                                    (x[0][0], int(x[0][1]))
                                                    for x in aligned
                                                ]

                                                if writePDF:
                                                    uz = 800
                                                    o = 0
                                                    for peak in peaks:
                                                        for p in peak:
                                                            pdf.drawString(
                                                                40,
                                                                uz,
                                                                "%35s: %d %.2f -> %d %d"
                                                                % (
                                                                    do[o][
                                                                        (
                                                                            1
                                                                            + do[
                                                                                o
                                                                            ].rfind("/")
                                                                        ) :
                                                                    ],
                                                                    int(p.NPeakCenter),
                                                                    p.NPeakCenterMin
                                                                    / 60.0,
                                                                    int(aligned[o][0]),
                                                                    aligned[o][1],
                                                                ),
                                                            )
                                                        uz -= 20
                                                        o += 1

                                                maxGroup = max(
                                                    aligned, key=itemgetter(1)
                                                )[1]
                                                minGroup = min(
                                                    aligned, key=itemgetter(1)
                                                )[1]

                                                groupedChromPeaks = []
                                                groupedChromPeaksAVGMz = []
                                                groupedChromPeaksAVGLMz = []
                                                groupedChromPeaksAVGTMz = []
                                                groupedChromPeaksAVGTimes = []
                                                groupedChromPeaksAVGNPeakScale = []
                                                groupedChromPeaksAVGLPeakScale = []

                                                for i in range(maxGroup + 1):
                                                    groupedChromPeaks.append({})
                                                    groupedChromPeaksAVGMz.append([])
                                                    groupedChromPeaksAVGLMz.append([])
                                                    groupedChromPeaksAVGTMz.append([])
                                                    groupedChromPeaksAVGTimes.append([])
                                                    groupedChromPeaksAVGNPeakScale.append(
                                                        []
                                                    )
                                                    groupedChromPeaksAVGLPeakScale.append(
                                                        []
                                                    )

                                                # calculate average values (RT and mz) for feature pairs in the sub-subclusters
                                                j = 0
                                                for k in partChromPeaks.keys():
                                                    for i in range(
                                                        len(partChromPeaks[k])
                                                    ):
                                                        if not (
                                                            k
                                                            in groupedChromPeaks[
                                                                aligned[j][1]
                                                            ]
                                                        ):
                                                            groupedChromPeaks[
                                                                aligned[j][1]
                                                            ][k] = []
                                                        groupedChromPeaks[
                                                            aligned[j][1]
                                                        ][k].append(
                                                            (
                                                                aligned[j],
                                                                partChromPeaks[k][i],
                                                            )
                                                        )
                                                        groupedChromPeaksAVGMz[
                                                            aligned[j][1]
                                                        ].append(
                                                            partChromPeaks[k][i].mz
                                                        )
                                                        groupedChromPeaksAVGLMz[
                                                            aligned[j][1]
                                                        ].append(
                                                            partChromPeaks[k][i].lmz
                                                        )
                                                        groupedChromPeaksAVGTMz[
                                                            aligned[j][1]
                                                        ].append(
                                                            partChromPeaks[k][i].tmz
                                                        )
                                                        groupedChromPeaksAVGTimes[
                                                            aligned[j][1]
                                                        ].append(
                                                            partChromPeaks[k][
                                                                i
                                                            ].NPeakCenterMin
                                                        )
                                                        groupedChromPeaksAVGNPeakScale[
                                                            aligned[j][1]
                                                        ].append(
                                                            partChromPeaks[k][
                                                                i
                                                            ].NPeakScale
                                                        )
                                                        groupedChromPeaksAVGLPeakScale[
                                                            aligned[j][1]
                                                        ].append(
                                                            partChromPeaks[k][
                                                                i
                                                            ].LPeakScale
                                                        )

                                                        j = j + 1

                                                assert j == len(aligned)
                                                assert (
                                                    len(groupedChromPeaks)
                                                    == len(groupedChromPeaksAVGMz)
                                                    == len(groupedChromPeaksAVGTimes)
                                                )

                                                # write results to data matrix and SQLite DB
                                                for i in range(minGroup, maxGroup + 1):
                                                    if len(groupedChromPeaks[i]) > 0:
                                                        avgmz = sum(
                                                            groupedChromPeaksAVGMz[i]
                                                        ) / len(
                                                            groupedChromPeaksAVGMz[i]
                                                        )
                                                        avglmz = sum(
                                                            groupedChromPeaksAVGLMz[i]
                                                        ) / len(
                                                            groupedChromPeaksAVGLMz[i]
                                                        )
                                                        avgtmz = sum(
                                                            groupedChromPeaksAVGTMz[i]
                                                        ) / len(
                                                            groupedChromPeaksAVGTMz[i]
                                                        )
                                                        avgtime = sum(
                                                            groupedChromPeaksAVGTimes[i]
                                                        ) / len(
                                                            groupedChromPeaksAVGTimes[i]
                                                        )
                                                        avgNPeakScale = sum(
                                                            groupedChromPeaksAVGNPeakScale[
                                                                i
                                                            ]
                                                        ) / len(
                                                            groupedChromPeaksAVGNPeakScale[
                                                                i
                                                            ]
                                                        )
                                                        avgLPeakScale = sum(
                                                            groupedChromPeaksAVGLPeakScale[
                                                                i
                                                            ]
                                                        ) / len(
                                                            groupedChromPeaksAVGLPeakScale[
                                                                i
                                                            ]
                                                        )

                                                        f.write(str(curNum))
                                                        f.write("\t")
                                                        f.write("")
                                                        f.write("\t")
                                                        f.write(str(avgmz))
                                                        f.write("\t")
                                                        f.write(str(avglmz))
                                                        f.write("\t")
                                                        f.write(str(avgtmz))
                                                        f.write("\t")
                                                        f.write(
                                                            "%.6f - %.6f"
                                                            % (
                                                                min(
                                                                    groupedChromPeaksAVGMz[
                                                                        i
                                                                    ]
                                                                ),
                                                                max(
                                                                    groupedChromPeaksAVGMz[
                                                                        i
                                                                    ]
                                                                ),
                                                            )
                                                        )
                                                        f.write("\t")
                                                        f.write(
                                                            "%.2f" % (avgtime / 60.0)
                                                        )
                                                        f.write("\t")
                                                        f.write(
                                                            "%.3f - %.3f"
                                                            % (
                                                                min(
                                                                    groupedChromPeaksAVGTimes[
                                                                        i
                                                                    ]
                                                                )
                                                                / 60.0,
                                                                max(
                                                                    groupedChromPeaksAVGTimes[
                                                                        i
                                                                    ]
                                                                )
                                                                / 60.0,
                                                            )
                                                        )
                                                        f.write("\t")
                                                        f.write(
                                                            "%.1f:%.1f"
                                                            % (
                                                                avgNPeakScale,
                                                                avgLPeakScale,
                                                            )
                                                        )
                                                        f.write("\t")
                                                        f.write("%s" % xCount)
                                                        f.write("\t")
                                                        f.write(str(cLoading))
                                                        f.write("\t")
                                                        f.write(scanEvent)
                                                        f.write("\t")
                                                        f.write(ionMode)
                                                        f.write("\t")
                                                        f.write(str(tracer))

                                                        f.write("\t")
                                                        f.write("")  # OGroup
                                                        f.write("\t")
                                                        f.write("")  # Ion
                                                        f.write("\t")
                                                        f.write("")  # Loss
                                                        f.write("\t")
                                                        f.write("")  # M

                                                        SQLInsert(
                                                            resDB.curs,
                                                            "GroupResults",
                                                            id=curNum,
                                                            mz=avgmz,
                                                            lmz=avglmz,
                                                            dmz=avgtmz,
                                                            rt=avgtime,
                                                            xn=xCount,
                                                            charge=cLoading,
                                                            ionisationMode=ionMode,
                                                            scanEvent=scanEvent,
                                                            tracer=str(tracer),
                                                        )

                                                        doublePeak = 0
                                                        for j in range(len(results)):
                                                            res = results[j].filePath
                                                            if (
                                                                res
                                                                in groupedChromPeaks[i]
                                                                and len(
                                                                    groupedChromPeaks[
                                                                        i
                                                                    ][res]
                                                                )
                                                                > 0
                                                            ):
                                                                for (
                                                                    peak
                                                                ) in groupedChromPeaks[
                                                                    i
                                                                ][res]:
                                                                    SQLInsert(
                                                                        resDB.curs,
                                                                        "FoundFeaturePairs",
                                                                        file=results[
                                                                            j
                                                                        ].fileName,
                                                                        featurePairID=peak[
                                                                            1
                                                                        ].id,
                                                                        featureGroupID=peak[
                                                                            1
                                                                        ].fGroupID,
                                                                        resID=curNum,
                                                                        NPeakScale=peak[
                                                                            1
                                                                        ].NPeakScale,
                                                                        LPeakScale=peak[
                                                                            1
                                                                        ].LPeakScale,
                                                                        NBorderLeft=peak[
                                                                            1
                                                                        ].NBorderLeft,
                                                                        NBorderRight=peak[
                                                                            1
                                                                        ].NBorderRight,
                                                                        LBorderLeft=peak[
                                                                            1
                                                                        ].LBorderLeft,
                                                                        LBorderRight=peak[
                                                                            1
                                                                        ].LBorderRight,
                                                                        areaN=peak[
                                                                            1
                                                                        ].NPeakArea,
                                                                        areaL=peak[
                                                                            1
                                                                        ].LPeakArea,
                                                                        featureType="foundPattern",
                                                                    )

                                                                doublePeak += (
                                                                    1
                                                                    if len(
                                                                        groupedChromPeaks[
                                                                            i
                                                                        ][res]
                                                                    )
                                                                    > 1
                                                                    else 0
                                                                )

                                                                if expPeakArea:
                                                                    f.write("\t")
                                                                    # Area of native, monoisotopic isotopolog
                                                                    f.write(
                                                                        ";".join(
                                                                            [
                                                                                str(
                                                                                    peak[
                                                                                        1
                                                                                    ].NPeakArea
                                                                                )
                                                                                for peak in groupedChromPeaks[
                                                                                    i
                                                                                ][res]
                                                                            ]
                                                                        )
                                                                    )

                                                                    f.write("\t")
                                                                    # Area of labeled isotopolog
                                                                    f.write(
                                                                        ";".join(
                                                                            [
                                                                                str(
                                                                                    peak[
                                                                                        1
                                                                                    ].LPeakArea
                                                                                )
                                                                                for peak in groupedChromPeaks[
                                                                                    i
                                                                                ][res]
                                                                            ]
                                                                        )
                                                                    )

                                                                if expApexIntensity:
                                                                    f.write("\t")
                                                                    # Abundance of native, monoisotopic isotopolog
                                                                    f.write(
                                                                        ";".join(
                                                                            [
                                                                                str(
                                                                                    peak[
                                                                                        1
                                                                                    ].NPeakAbundance
                                                                                )
                                                                                for peak in groupedChromPeaks[
                                                                                    i
                                                                                ][res]
                                                                            ]
                                                                        )
                                                                    )

                                                                    f.write("\t")
                                                                    # Abundance of labeled isotopolog
                                                                    f.write(
                                                                        ";".join(
                                                                            [
                                                                                str(
                                                                                    peak[
                                                                                        1
                                                                                    ].LPeakAbundance
                                                                                )
                                                                                for peak in groupedChromPeaks[
                                                                                    i
                                                                                ][res]
                                                                            ]
                                                                        )
                                                                    )

                                                                if expPeakSNR:
                                                                    f.write("\t")
                                                                    # SNR of native, monoisotopic isotopolog
                                                                    f.write(
                                                                        ";".join(
                                                                            [
                                                                                str(
                                                                                    peak[
                                                                                        1
                                                                                    ].NSNR
                                                                                )
                                                                                for peak in groupedChromPeaks[
                                                                                    i
                                                                                ][res]
                                                                            ]
                                                                        )
                                                                    )

                                                                    f.write("\t")
                                                                    # SNR of labeled isotopolog
                                                                    f.write(
                                                                        ";".join(
                                                                            [
                                                                                str(
                                                                                    peak[
                                                                                        1
                                                                                    ].LSNR
                                                                                )
                                                                                for peak in groupedChromPeaks[
                                                                                    i
                                                                                ][res]
                                                                            ]
                                                                        )
                                                                    )

                                                                if True:
                                                                    f.write("\t")
                                                                    # Feature ID
                                                                    f.write(
                                                                        ";".join(
                                                                            [
                                                                                str(
                                                                                    peak[
                                                                                        1
                                                                                    ].id
                                                                                )
                                                                                for peak in groupedChromPeaks[
                                                                                    i
                                                                                ][res]
                                                                            ]
                                                                        )
                                                                    )

                                                                    f.write("\t")
                                                                    # Group ID
                                                                    f.write(
                                                                        ";".join(
                                                                            [
                                                                                str(
                                                                                    peak[
                                                                                        1
                                                                                    ].fGroupID
                                                                                )
                                                                                for peak in groupedChromPeaks[
                                                                                    i
                                                                                ][res]
                                                                            ]
                                                                        )
                                                                    )

                                                            else:
                                                                if expPeakArea:
                                                                    f.write("\t\t")
                                                                if expApexIntensity:
                                                                    f.write("\t\t")
                                                                if expPeakSNR:
                                                                    f.write("\t\t")
                                                                if True:
                                                                    f.write("\t\t")
                                                        f.write("\t%d" % doublePeak)
                                                        f.write("\n")

                                                        curNum = curNum + 1

                                                if writePDF:
                                                    pdf.showPage()

                                            except Exception as e:
                                                import traceback

                                                traceback.print_exc()
                                                logging.error(str(traceback))

                                                logging.error("Error: %s" % (str(e)))
                                            grp = grp + 1

                                            for res in results:
                                                if res.filePath in toDel:
                                                    toDel[res.filePath] = [
                                                        a for a in toDel[res.filePath]
                                                    ]
                                                    toDel[res.filePath].sort()
                                                    toDel[res.filePath].reverse()
                                                    for i in toDel[res.filePath]:
                                                        res.featurePairs.pop(i)
                                                        chromPeaksAct = (
                                                            chromPeaksAct - 1
                                                        )
                                                        totalChromPeaksProc = (
                                                            totalChromPeaksProc + 1
                                                        )

                                            if pwTextSet is not None:
                                                # Log time used for bracketing
                                                elapsed = (time.time() - start) / 60.0
                                                hours = ""
                                                if elapsed >= 60.0:
                                                    hours = "%d hours " % (
                                                        elapsed // 60
                                                    )
                                                mins = "%.2f mins" % (elapsed % 60.0)
                                                pwTextSet.put(
                                                    Bunch(
                                                        mes="text",
                                                        val="<p align='right' >%s%s elapsed</p>\n\n\n\nBracketing results\n%d feature pairs (%d individual pairs processed).. (Ionmode: %s, XCount: %s, Charge: %d)"
                                                        % (
                                                            hours,
                                                            mins,
                                                            curNum,
                                                            totalChromPeaksProc,
                                                            ionMode,
                                                            xCount,
                                                            cLoading,
                                                        ),
                                                    )
                                                )
                                    curAllMZs = []
            if writePDF:
                pdf.save()

            import uuid
            import platform
            import datetime

            identifier = "%s_%s_%s" % (
                str(uuid.uuid1()),
                str(platform.node()),
                str(datetime.datetime.now()),
            )

            f.write(
                "## MetExtract II %s\n"
                % (
                    Bunch(
                        MetExtractVersion=meVersion,
                        RVersion=rVersion,
                        UUID_ext=identifier,
                    )
                    .dumpAsJSon()
                    .replace('"', "'")
                )
            )
            f.write(
                "## Individual files processing parameters %s\n"
                % (generalProcessingParams.dumpAsJSon().replace('"', "'"))
            )
            processingParams = Bunch()

            xCountsFmt = []
            try:
                for xn in sorted(xCounts):
                    if xn - 1 in xCounts and xn + 1 in xCounts:
                        xCountsFmt.append("-")
                    else:
                        if xn - 1 not in xCounts:
                            xCountsFmt.append(", ")
                        xCountsFmt.append(xn)
                xCountsFmt.pop(0)
                seen = set()
                seen_add = seen.add
                xCountsFmt = [
                    x for x in xCountsFmt if not (x in seen or seen_add(x)) or x == ", "
                ]
            except Exception:
                xCountsFmt = xCounts

            processingParams.FPBracketing_xCounts = "".join(
                [str(t) for t in xCountsFmt]
            )
            processingParams.FPBracketing_groupSizePPM = groupSizePPM
            processingParams.FPBracketing_positiveScanEvent = positiveScanEvent
            processingParams.FPBracketing_negativeScanEvent = negativeScanEvent
            processingParams.FPBracketing_maxTimeDeviation = maxTimeDeviation
            processingParams.FPBracketing_maxLoading = maxLoading
            processingParams.FPBracketing_resultsFile = file
            processingParams.FPBracketing_align = align
            if align:
                processingParams.FPBracketing_nPolynom = nPolynom
            f.write(
                "## Bracketing files processing parameters %s\n"
                % (processingParams.dumpAsJSon().replace('"', "'"))
            )

        exportAsFeatureML.convertMEMatrixToFeatureMLSepPolarities(file, postfix="_ab")

    except Exception as ex:
        import traceback

        traceback.print_exc()
        print(str(traceback))

        logging.error("Error during bracketing of files")

    finally:
        resDB.conn.commit()
        resDB.curs.close()
        resDB.conn.close()

        logging.info("Bracketing done...")


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


# store used configuration to DB file
def writeMetaboliteGroupingConfigToDB(
    curs, minConnectionsInFiles, minConnectionRate, groups
):
    SQLInsert(curs, "config", key="MEGROUP_groups", value=str(groups))
    SQLInsert(
        curs,
        "config",
        key="MEGROUP_minConnectionsInFiles",
        value=str(minConnectionsInFiles),
    )
    SQLInsert(
        curs, "config", key="MEGROUP_minConnectionRate", value=str(minConnectionRate)
    )


def getBordersFor(curs, fID, file):
    for row in curs.execute(
        "SELECT NBorderLeft, NBorderRight, LBorderLeft, LBorderRight FROM foundfeaturepairs WHERE file=? and resID=?",
        (file, fID),
    ):
        return (float(row[0]), float(row[1]), float(row[2]), float(row[3]))
    return None


def getPeak(times, rt, borderleft, borderright):
    ## find best matching rt
    ind, tim = min(enumerate(times), key=lambda x: abs(x[1] - rt))
    minind = int(max(0, ind - int(borderleft)))
    maxind = int(min(len(times) - 1, ind + borderright))

    return (minind, maxind)


def calculateMetaboliteGroups(
    file="./results.tsv",
    groups=[],
    eicPPM=10.0,
    maxAnnotationTimeWindow=0.05,
    minConnectionsInFiles=1,
    minConnectionRate=0.4,
    minPeakCorrelation=0.85,
    useRatio=False,
    runIdentificationInstance=None,
    pwMaxSet=None,
    pwValSet=None,
    pwTextSet=None,
    cpus=1,
    toFile=None,
):
    if toFile is None:
        toFile = file

    logging.info("Starting convoluting feature pairs..")

    resDB = Bunch(conn=None, curs=None)
    resDB.conn = connect(file + getDBSuffix())
    # conn.execute('''PRAGMA synchronous = OFF''')
    # conn.execute('''PRAGMA journal_mode = OFF''')
    resDB.curs = resDB.conn.cursor()

    try:
        # Select only those groups that shall be used for metabolite group grouping
        ## load mzxml files for improved metabolite convolution
        useGroups = []
        useGroupsForConfig = []
        numFilesForConvolution = 0
        for group in groups:
            if group.useForMetaboliteGrouping:
                useGroups.append(group)
                useGroupsForConfig.append(str(group.name) + ":" + str(group.files))

                for fi in group.files:
                    numFilesForConvolution += 1

        if pwMaxSet is not None:
            pwMaxSet.put(Bunch(mes="max", val=numFilesForConvolution))
        if pwValSet is not None:
            pwValSet.put(Bunch(mes="value", val=0))

        # connect to results db
        writeMetaboliteGroupingConfigToDB(
            resDB.curs,
            minConnectionsInFiles,
            minConnectionRate,
            str(useGroupsForConfig).replace("'", "").replace('"', ""),
        )

        # read results table
        table = TableUtils.readFile(file, delim="\t")
        cols = [col.getName() for col in table.getColumns()]

        # check if necessary columns are present
        ## overall columns
        assert "Num" in cols
        assert "OGroup" in cols
        assert "Ion" in cols
        assert "Loss" in cols
        assert "M" in cols
        assert "doublePeak" in cols

        doublePeaks = 0
        for row in table.getData(cols=["COUNT(*)"], where="doublePeak>0"):
            doublePeaks = int(row)

        if doublePeaks > 0:
            logging.info("  found double peaks: %d" % doublePeaks)
        writeCols = []
        if doublePeaks > 0:
            TableUtils.saveFile(
                table,
                toFile.replace(".tsv", ".doublePeaks.tsv"),
                cols=writeCols,
                order="OGroup, MZ, Xn",
                where="doublePeak>0",
            )

        table.executeSQL("DELETE FROM :table: WHERE doublePeak>0")
        ## file specific columns

        for group in useGroups:
            for f in group.files:
                fShort = f[f.rfind("/") + 1 : f.rfind(".")]
                if fShort[0].isdigit():
                    fShort = "_" + fShort

                assert "%s_FID" % fShort in cols
                assert "%s_GroupID" % fShort in cols

        # fetch all feature pairs from the results file
        nodes = {}
        for fpNum, mz, lmz, rt, xn, charge, scanEvent, ionMode in table.getData(
            cols=[
                "Num",
                "MZ",
                "L_MZ",
                "RT",
                "Xn",
                "Charge",
                "ScanEvent",
                "Ionisation_Mode",
            ]
        ):
            nodes[fpNum] = Bunch(
                fpNum=fpNum,
                mz=mz,
                lmz=lmz,
                rt=rt,
                xn=xn,
                charge=charge,
                scanEvent=scanEvent,
                ionMode=ionMode,
                oGroup=None,
                correlationToOtherFPs={},
            )

        ## generate all correlations and sil ratios
        fileCorrelations = {}
        fileSILRatios = {}

        borders = {}
        for row in resDB.curs.execute(
            "SELECT NBorderLeft, NBorderRight, LBorderLeft, LBorderRight, file, resID FROM foundfeaturepairs"
        ):
            fil = str(row[4])
            resID = int(row[5])
            if fil not in borders.keys():
                borders[fil] = {}
            borders[fil][resID] = (
                float(row[0]),
                float(row[1]),
                float(row[2]),
                float(row[3]),
            )

        ## Iterate all files
        processedFiles = 0
        procObjects = []
        for group in useGroups:
            for fi in group.files:
                if fi not in fileCorrelations.keys():
                    # if pwTextSet is not None: pwTextSet.put(Bunch(mes="text", val="Convoluting feature pairs in file %s" % (fi)))
                    # if pwValSet is not None: pwValSet.put(Bunch(mes="value", val=processedFiles))

                    s = ConvoluteFPsInFile(
                        eicPPM, fi, maxAnnotationTimeWindow, nodes, borders
                    )
                    procObjects.append(s)

                    # processedFiles += 1

        from multiprocessing import Pool, Manager

        p = Pool(
            processes=cpus, maxtasksperchild=1
        )  # only in python >=2.7; experimental
        manager = Manager()
        queue = manager.Queue()

        start = time.time()
        # start the multiprocessing pool
        res = p.imap_unordered(processConvoluteFPsInFile, procObjects)

        # wait until all subprocesses have finished re-integrating their respective LC-HRMS data file
        loop = True
        freeSlots = range(min(len(procObjects), cpus))
        assignedThreads = {}
        while loop:
            completed = res._index
            if completed == len(procObjects):
                loop = False
            else:
                if pwValSet is not None:
                    pwValSet.put(Bunch(mes="value", val=completed))

                elapsed = (time.time() - start) / 60.0
                hours = ""
                if elapsed >= 60.0:
                    hours = "%d hours " % (elapsed // 60)
                mins = "%.2f mins" % (elapsed % 60.0)

                if pwTextSet is not None:
                    pwTextSet.put(
                        Bunch(
                            mes="text",
                            val="<p align='right' >%s%s elapsed</p>\n\n%d / %d files done (%d parallel)"
                            % (hours, mins, completed, len(procObjects), cpus),
                        )
                    )

                time.sleep(0.5)

        p.close()

        ## fetch multiprocessing results
        for rei, re in enumerate(res):
            fileCorrelations[procObjects[rei].fi] = re[0]
            fileSILRatios[procObjects[rei].fi] = re[1]

        logging.info("Testing SIL ratios")
        ## test all feature pair pairs for co-elution and similar SIL ratio
        connections = {}
        for fpNumA in nodes.keys():
            connections[fpNumA] = {}
        for fpNumA in nodes.keys():
            nodeA = nodes[fpNumA]

            for fpNumB in nodes.keys():
                nodeB = nodes[fpNumB]

                if (
                    nodeA.fpNum != nodeB.fpNum
                    and nodeA.fpNum < nodeB.fpNum
                    and abs(nodeB.rt - nodeA.rt) <= maxAnnotationTimeWindow
                ):
                    ## test feature pair pair convolution in all samples
                    allCorrelations = []
                    allSILRatios = []
                    for group in useGroups:
                        for fi in group.files:
                            if (
                                fi in fileCorrelations.keys()
                                and fpNumA in fileCorrelations[fi].keys()
                                and fpNumB in fileCorrelations[fi][fpNumA].keys()
                            ):
                                co = fileCorrelations[fi][fpNumA][fpNumB]
                                allCorrelations.append(co)
                            if (
                                fi in fileSILRatios.keys()
                                and fpNumA in fileSILRatios[fi].keys()
                                and fpNumB in fileSILRatios[fi][fpNumA].keys()
                            ):
                                silRatio = fileSILRatios[fi][fpNumA][fpNumB]
                                allSILRatios.append(silRatio)

                    if len(allCorrelations) > 0 and len(allSILRatios) > 0:
                        if sum(
                            [co >= minPeakCorrelation for co in allCorrelations]
                        ) >= minConnectionRate * len(allCorrelations) and (
                            not useRatio
                            or sum([rat for rat in allSILRatios])
                            >= minConnectionRate * len(allSILRatios)
                        ):
                            nodes[nodeA.fpNum].correlationToOtherFPs[nodeB.fpNum] = True
                            nodes[nodeB.fpNum].correlationToOtherFPs[nodeA.fpNum] = True

                        m = mean(allCorrelations)
                        connections[fpNumA][fpNumB] = m
                        connections[fpNumB][fpNumA] = m
                    else:
                        connections[fpNumA][fpNumB] = 0
                        connections[fpNumB][fpNumA] = 0

        logging.info(" .. done")

        nodes2 = {}
        for fpNumA in nodes.keys():
            if fpNumA not in nodes2.keys():
                nodes2[fpNumA] = []
            for fpNumB in nodes[fpNumA].correlationToOtherFPs.keys():
                nodes2[fpNumA].append(fpNumB)

        for k in nodes2.keys():
            uniq = []
            for u in nodes2[k]:
                if u not in uniq:
                    uniq.append(u)
            nodes2[k] = uniq

        # get subgraphs from the feature pair graph. Each subgraph represents one convoluted
        # feature group
        tGroups = getSubGraphs(nodes2)

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

            # split the HCA tree in sub-clusters
            def checkSubCluster(tree, hca, corrThreshold, cutOffMinRatio):
                if isinstance(tree, HCA_general.HCALeaf):
                    return False
                elif isinstance(tree, HCA_general.HCAComposite):
                    corrs = hca.getLinkFor(tree)
                    inds = hca.getIndsFor(tree)

                    aboveThreshold = sum(
                        [
                            corr > corrThreshold
                            for i, corr in enumerate(corrs)
                            if i in inds
                        ]
                    )

                    # hc.plotTree(tree)
                    # print(corrs)
                    # print([corrs[i] for i in inds])
                    # print(aboveThreshold)
                    # print(".....")
                    # print("")
                    return not (aboveThreshold * 1.0 / len(inds)) >= cutOffMinRatio

                else:
                    raise Error("Unknown if-branch")

            subClusts = hc.splitTreeWithCallback(
                tree,
                CallBackMethod(
                    _target=checkSubCluster,
                    corrThreshold=minPeakCorrelation,
                    cutOffMinRatio=minConnectionRate,
                ).getRunMethod(),
                recursive=True,
            )

            # convert the subclusters into arrays of feature pairs belonging to the same metabolite
            return [
                [leaf.getID() for leaf in subClust.getLeaves()]
                for subClust in subClusts
            ]

        groups = []
        done = 0
        for tGroup in tGroups:
            sGroups = splitGroupWithHCA(tGroup, connections)
            groups.extend(sGroups)

            done = done + 1

        # Separate feature pair clusters; softer
        curGroup = 1
        for tGroup in groups:
            for i in range(len(tGroup)):
                table.setData(
                    cols=["OGroup"], vals=[curGroup], where="Num = %s" % (tGroup[i])
                )
                resDB.curs.execute(
                    "UPDATE GroupResults SET OGroup=%d WHERE id = %s"
                    % (curGroup, tGroup[i])
                )

            curGroup += 1

        logging.info("Annotating ions in feature groups")
        # Annotate groups with common adducts and in-source fragments
        if runIdentificationInstance is not None:
            groups = defaultdict(list)
            for row in table.getData(
                cols=[
                    "Num",
                    "OGroup",
                    "MZ",
                    "Ionisation_Mode",
                    "Charge",
                    "Ion",
                    "Xn",
                    "Loss",
                    "M",
                ]
            ):
                num, ogrp, mz, ionMode, loading, adducts, xCount, fDesc, ms = row
                groups[ogrp].append(
                    ChromPeakPair(
                        id=num,
                        fGroupID=ogrp,
                        mz=mz,
                        ionMode=ionMode,
                        loading=loading,
                        adducts=[],
                        heteroAtomsFeaturePairs=[],
                        xCount=xCount,
                        fDesc=[],
                    )
                )

            for ogrp in groups.keys():
                chromPeaks = {}
                for fp in groups[ogrp]:
                    chromPeaks[fp.id] = fp

                runIdentificationInstance.annotateChromPeaks(
                    chromPeaks.keys(), chromPeaks
                )

            for ogrp in groups.keys():
                for fp in groups[ogrp]:
                    table.setData(
                        cols=["Ion", "Loss", "M"],
                        vals=[
                            ",".join([str(a) for a in fp.adducts]),
                            ",".join([str(a) for a in fp.fDesc]),
                            ",".join([str(a) for a in fp.Ms]),
                        ],
                        where="Num = %d" % (fp.id),
                    )
                    resDB.curs.execute(
                        "UPDATE GroupResults SET Ion='%s', Loss='%s', M='%s' WHERE id = %d"
                        % (
                            ",".join([str(a) for a in fp.adducts]),
                            ",".join([str(a) for a in fp.fDesc]),
                            ",".join([str(a) for a in fp.Ms]),
                            fp.id,
                        )
                    )

        resDB.curs.execute(
            "DELETE FROM GroupResults WHERE id NOT IN (%s)"
            % (",".join([str(num) for num in table.getData(cols=["Num"])]))
        )

        ## reassign feature group ids
        resDB.curs.execute("UPDATE GroupResults SET OGroup='X'||OGroup")
        table.executeSQL("UPDATE :table: SET OGroup='X'||OGroup")
        oGrps = []
        curGrp = 1
        curFP = 1
        for row in resDB.curs.execute(
            "SELECT OGroup, AVG(rt) FROM GroupResults GROUP BY OGroup ORDER BY AVG(rt)"
        ):
            oGrps.append(str(row[0]))
        for tgrp in oGrps:
            resDB.curs.execute(
                "UPDATE GroupResults SET OGroup='%d' WHERE OGroup='%s'" % (curGrp, tgrp)
            )
            table.setData(cols=["OGroup"], vals=[curGrp], where="OGroup='%s'" % (tgrp))
            curGrp = curGrp + 1

        processingParams = Bunch()
        processingParams.MEConvoluting_groups = (
            str(useGroupsForConfig).replace("'", "").replace('"', "")
        )
        processingParams.MEConvoluting_connThreshold = minConnectionRate
        processingParams.MEConvoluting_minPeakCorrelation = minPeakCorrelation
        table.addComment(
            "## Convolution FPs processing parameters %s"
            % (processingParams.dumpAsJSon().replace('"', "'"))
        )

        writeCols = []
        for column in table.getColumns():
            if not column.name.endswith("_CorrelatesTo"):
                writeCols.append(column.name)

        TableUtils.saveFile(table, toFile, cols=writeCols, order="OGroup, MZ, Xn")

        exportAsFeatureML.convertMEMatrixToFeatureMLSepPolarities(file, postfix="_ac")

    except Exception as ex:
        import traceback

        traceback.print_exc()

        raise ex

    finally:
        resDB.conn.commit()
        resDB.curs.close()
        resDB.conn.close()


from utils import mapArrayToRefTimes


class ConvoluteFPsInFile:
    def __init__(self, eicPPM, fi, maxAnnotationTimeWindow, nodes, borders):
        self.eicPPM = eicPPM
        self.fi = fi
        self.maxAnnotationTimeWindow = maxAnnotationTimeWindow
        self.nodes = nodes
        self.borders = borders

        self.fileCorrelations = None
        self.fileSILRatios = None

    def getConvolutionInFile(self):
        startProc = time.time()

        eicPPM = self.eicPPM
        fi = self.fi
        maxAnnotationTimeWindow = self.maxAnnotationTimeWindow
        nodes = self.nodes
        borders = self.borders

        fileCorrelations = {}
        fileSILRatios = {}

        logging.info("  Convoluting feature pairs in file %s" % (fi))
        b = fi.replace("\\", "/")
        fiName = ""
        if ".mzXML" in b:
            fiName = b[(b.rfind("/") + 1) : b.rfind(".mzXML")]
        if ".mzML" in b:
            fiName = b[(b.rfind("/") + 1) : b.rfind(".mzML")]
        mzXML = Chromatogram()
        mzXML.parse_file(fi)

        filterLines = mzXML.getFilterLines(
            includeMS1=True,
            includeMS2=False,
            includePosPolarity=True,
            includeNegPolarity=True,
        )

        for fpNumA in nodes.keys():
            fileCorrelations[fpNumA] = {}
            fileSILRatios[fpNumA] = {}

        for fpNumA in nodes.keys():
            nodeA = nodes[fpNumA]

            for fpNumB in nodes.keys():
                nodeB = nodes[fpNumB]

                if (
                    nodeA.fpNum != nodeB.fpNum
                    and nodeA.fpNum < nodeB.fpNum
                    and abs(nodeB.rt - nodeA.rt) <= maxAnnotationTimeWindow
                ):
                    ra = None
                    if (
                        fiName in borders.keys()
                        and nodeA.fpNum in borders[fiName].keys()
                    ):
                        ra = borders[fiName][nodeA.fpNum]
                    rb = None
                    if (
                        fiName in borders.keys()
                        and nodeB.fpNum in borders[fiName].keys()
                    ):
                        rb = borders[fiName][nodeB.fpNum]

                    if ra != None and rb != None:
                        nanbl, nanbr, nalbl, nalbr = ra
                        nbnbl, nbnbr, nblbl, nblbr = rb

                        if (
                            nodeA.scanEvent in filterLines
                            and nodeB.scanEvent in filterLines
                        ):
                            meanRT = mean([nodeA.rt, nodeB.rt])

                            eicAL, timesA, scanIdsA, mzs = mzXML.getEIC(
                                nodeA.lmz,
                                ppm=eicPPM,
                                filterLine=nodeA.scanEvent,
                                removeSingles=True,
                                intThreshold=0,
                                useMS1=True,
                                useMS2=False,
                                startTime=meanRT * 60 - 120,
                                endTime=meanRT * 60 + 120,
                            )
                            eicBL, timesB, scanIdsB, mzs = mzXML.getEIC(
                                nodeB.lmz,
                                ppm=eicPPM,
                                filterLine=nodeB.scanEvent,
                                removeSingles=True,
                                intThreshold=0,
                                useMS1=True,
                                useMS2=False,
                                startTime=meanRT * 60 - 120,
                                endTime=meanRT * 60 + 120,
                            )

                            eicA, timesA, scanIdsA, mzs = mzXML.getEIC(
                                nodeA.mz,
                                ppm=eicPPM,
                                filterLine=nodeA.scanEvent,
                                removeSingles=True,
                                intThreshold=0,
                                useMS1=True,
                                useMS2=False,
                                startTime=meanRT * 60 - 120,
                                endTime=meanRT * 60 + 120,
                            )
                            eicB, timesB, scanIdsB, mzs = mzXML.getEIC(
                                nodeB.mz,
                                ppm=eicPPM,
                                filterLine=nodeB.scanEvent,
                                removeSingles=True,
                                intThreshold=0,
                                useMS1=True,
                                useMS2=False,
                                startTime=meanRT * 60 - 120,
                                endTime=meanRT * 60 + 120,
                            )

                            timesMin = [t for t in timesA]
                            timesMin.extend([t for t in timesB])
                            timesMin = sorted(timesMin)

                            eicAL = mapArrayToRefTimes(eicAL, timesA, timesMin)
                            eicBL = mapArrayToRefTimes(eicBL, timesB, timesMin)
                            eicA = mapArrayToRefTimes(eicA, timesA, timesMin)
                            eicB = mapArrayToRefTimes(eicB, timesB, timesMin)

                            timesMin = [t / 60.0 for t in timesMin]

                            ## A) Test correlation of different feature pairs
                            try:
                                lI, rI = getPeak(
                                    timesMin,
                                    meanRT,
                                    min(nanbl, nbnbl) * 2,
                                    min(nanbr, nbnbr) * 2,
                                )

                                eicAC = eicAL[lI:rI]
                                eicBC = eicBL[lI:rI]

                                co = corr(eicAC, eicBC)

                                if not (isnan(co)) and co != None:
                                    fileCorrelations[fpNumA][fpNumB] = co
                            except Exception as err:
                                logging.error(
                                    "  Error during convolution of feature pairs (Peak-correlation, Nums: %s and %s, message: %s).."
                                    % (fpNumA, fpNumB, err.message)
                                )

                            ## B) Test similarity of native:labeled ratio
                            try:
                                lISIL, rISIL = getPeak(
                                    timesMin,
                                    nodeA.rt,
                                    min(nanbl, nanbr) * 0.8,
                                    min(nalbl, nalbr) * 0.8,
                                )

                                eicANCSIL = eicA[lISIL:rISIL]
                                eicALCSIL = eicAL[lISIL:rISIL]

                                folds = [
                                    eicANCSIL[i] / eicALCSIL[i]
                                    for i in range(len(eicANCSIL))
                                    if eicALCSIL[i] > 0
                                ]
                                ma = mean(folds)
                                sa = sd(folds)

                                lISIL, rISIL = getPeak(
                                    timesMin,
                                    nodeB.rt,
                                    min(nanbl, nanbr) * 0.8,
                                    min(nalbl, nalbr) * 0.8,
                                )

                                eicBNCSIL = eicB[lISIL:rISIL]
                                eicBLCSIL = eicBL[lISIL:rISIL]

                                folds = [
                                    eicBNCSIL[i] / eicBLCSIL[i]
                                    for i in range(len(eicBNCSIL))
                                    if eicBLCSIL[i] > 0
                                ]
                                mb = mean(folds)
                                sb = sd(folds)

                                if (
                                    ma != None
                                    and mb != None
                                    and sa != None
                                    and sb != None
                                ):
                                    silRatioFold = 1
                                    try:
                                        silRatioFold = max([ma, mb]) / min([ma, mb])
                                    except RuntimeWarning as ex:
                                        pass
                                    # print(ma, mb, sa, sb, max([ma, mb]), min([ma, mb]))
                                    fileSILRatios[fpNumA][fpNumB] = (
                                        silRatioFold <= 1 + max(0.5, 3 * mean([sa, sb]))
                                    )

                                    # print(fiName, nodeA.mz, nodeB.mz, meanRT, silRatioFold, co)
                            except Exception as err:
                                logging.error(
                                    "  Error during convolution of feature pairs (SIL-ratio, Nums: %s and %s, message: %s).."
                                    % (fpNumA, fpNumB, err.message)
                                )

        logging.info(
            "  Finished convoluting feature pairs in file %s (%.1f minutes)"
            % (fi, (time.time() - startProc) / 60.0)
        )
        return (fileCorrelations, fileSILRatios)


def processConvoluteFPsInFile(cif):
    se = cif.getConvolutionInFile()
    return se
