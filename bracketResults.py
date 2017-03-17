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

from utils import ChromPeakPair, Bunch, SQLInsert, natSort, get_main_dir, getSubGraphsFromDictDict, CallBackMethod
from runIdentification import ChromPeakPair, getDBSuffix
from TableUtils import TableUtils
from MZHCA import HierarchicalClustering, cutTreeSized

import base64
from pickle import dumps, loads

import time

import HCA_general

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
def writeConfigToDB(curs, align, file, groupSizePPM, maxLoading, maxTimeDeviation, xCounts, nPolynom,
                    negativeScanEvent, positiveScanEvent, rVersion, meVersion):
    SQLInsert(curs, "config", key="MEVersion", value=str(meVersion))
    SQLInsert(curs, "config", key="RVersion", value=str(rVersion))

    SQLInsert(curs, "config", key="FPBRACK_xCounts", value=str(xCounts))
    SQLInsert(curs, "config", key="FPBRACK_groupSizePPM", value=str(groupSizePPM))
    SQLInsert(curs, "config", key="FPBRACK_positiveScanEvent", value=str(positiveScanEvent))
    SQLInsert(curs, "config", key="FPBRACK_negativeScanEvent", value=str(negativeScanEvent))
    SQLInsert(curs, "config", key="FPBRACK_maxTimeDeviation", value=str(maxTimeDeviation))
    SQLInsert(curs, "config", key="FPBRACK_maxLoading", value=str(maxLoading))
    SQLInsert(curs, "config", key="FPBRACK_file", value=str(file))
    SQLInsert(curs, "config", key="FPBRACK_align", value=str(align))
    SQLInsert(curs, "config", key="FPBRACK_nPolynom", value=str(nPolynom))

# bracket results
def bracketResults(indGroups, xCounts, groupSizePPM, positiveScanEvent=None, negativeScanEvent=None,
                 maxTimeDeviation=0.36 * 60, maxLoading=1, file="./results.tsv", align=True, nPolynom=1,
                 pwMaxSet=None, pwValSet=None, pwTextSet=None, rVersion="", meVersion="", generalProcessingParams=Bunch(), start=0):


    # create results SQLite tables
    resDB=Bunch(conn=None, curs=None)
    if os.path.exists(file+getDBSuffix()) and os.path.isfile(file+getDBSuffix()):
        os.remove(file+getDBSuffix())
    resDB.conn=connect(file+getDBSuffix())
    resDB.curs=resDB.conn.cursor()

    try:
        writePDF = False
        # used for debug purposes
        colos = [colors.red]

        cpf = get_main_dir() + "./XICAlignment.r"       # initialise chromatographic alignment script
        xicAlign = XICAlignment(cpf)


        resDB.curs.execute("DROP TABLE IF EXISTS GroupResults")
        resDB.curs.execute("CREATE TABLE GroupResults(id INTEGER PRIMARY KEY, mz FLOAT, lmz FLOAT, dmz FLOAT, rt FLOAT, xn INTEGER, charge INTEGER, scanEvent TEXT, ionisationMode TEXT, tracer TEXT, OGroup INTEGER, Ion TEXT, LOSS TEXT, M TEXT)")

        resDB.curs.execute("DROP TABLE IF EXISTS FoundFeaturePairs")
        resDB.curs.execute("CREATE TABLE FoundFeaturePairs(resID INTEGER, file TEXT, featurePairID INTEGER, featureGroupID INTEGER, areaN FLOAT, areaL FLOAT, featureType TEXT)")

        resDB.curs.execute("DROP TABLE IF EXISTS FileMapping")
        resDB.curs.execute("CREATE TABLE FileMapping(fileName TEXT, filePath TEXT, groupID INTEGER)")

        resDB.curs.execute("DROP TABLE IF EXISTS FileGroups")
        resDB.curs.execute("CREATE TABLE FileGroups(groupName TEXT, id INTEGER PRIMARY KEY)")

        resDB.curs.execute("DROP TABLE IF EXISTS config")
        resDB.curs.execute("CREATE TABLE config(id INTEGER PRIMARY KEY AUTOINCREMENT, key TEXT, value TEXT)")

        writeConfigToDB(resDB.curs, align, file, groupSizePPM, maxLoading, maxTimeDeviation, xCounts, nPolynom,
                        negativeScanEvent, positiveScanEvent, rVersion, meVersion)


        results = []

        # generate results for each input file
        i=1
        for key in natSort(indGroups.keys()):
            files=indGroups[key]
            SQLInsert(resDB.curs, "FileGroups", groupName=key, id=i)

            for ident in files:
                conn = connect(ident + getDBSuffix())
                curs = conn.cursor()
                fname=ident
                if ".mzxml" in ident.lower():
                    fname=ident[max(ident.rfind("/") + 1, ident.rfind("\\") + 1):ident.lower().rfind(".mzxml")]
                if ".mzml" in ident.lower():
                    fname=ident[max(ident.rfind("/") + 1, ident.rfind("\\") + 1):ident.lower().rfind(".mzml")]
                b=Bunch(filePath=ident, fileName=fname,featurePairs=[], conn=conn, curs=curs)

                results.append(b)

                SQLInsert(resDB.curs, "FileMapping", fileName=b.fileName, filePath=b.filePath, groupID=i)

            i+=1

        with open(file, "wb") as f:
            # initialise TSV results file

            f.write("Num\tComment\tMZ\tL_MZ\tD_MZ\tMZ_Range\tRT\tRT_Range\tXn\tCharge\tScanEvent\tIonisation_Mode\tTracer\tOGroup\tIon\tLoss\tM")
            for res in results:
                f.write("\t" + res.fileName + "_Area_N")
                f.write("\t" + res.fileName + "_Area_L")
                f.write("\t" + res.fileName + "_Abundance_N")
                f.write("\t" + res.fileName + "_Abundance_L")
                f.write("\t" + res.fileName + "_FID")
                f.write("\t" + res.fileName + "_GroupID")
                f.write("\t" + res.fileName + "_CorrelatesTo")
            for res in results:
                f.write("\t" + res.fileName + "_fold")

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
            rows = []
            for row in res.curs.execute("SELECT key, value FROM config WHERE key='metabolisationExperiment'"):
                rows.append(str(row[1]))
            assert len(rows) == 1
            isMetabolisationExperiment = rows[0].lower() == "true"

            ## TODO this used the last open file. May be invalid. improve
            if isMetabolisationExperiment:
                for row in res.curs.execute("SELECT id, name, deltaMZ FROM tracerConfiguration"):
                    tracerName = str(row[1])
                    dmz = float(row[2])
                    if tracerName not in tracersDeltaMZ:
                        tracersDeltaMZ[tracerName] = dmz

                    if tracersDeltaMZ[tracerName] != dmz:
                        logging.warning("Warning: Tracers have not been configured identical in all measurement files")
            else:
                rows = []
                ## TODO this used the last open file. May be invalid. improve
                for row in res.curs.execute("SELECT key, value FROM config WHERE key='xOffset'"):
                    rows.append(str(row[1]))
                assert len(rows) == 1

                dmz = float(rows[0])
                if "FLE" not in tracersDeltaMZ:
                    tracersDeltaMZ["FLE"] = dmz

                if tracersDeltaMZ["FLE"] != dmz:
                    logging.warning("Warning: Tracers have not been configured identical in all measurement files")

            if pwMaxSet is not None: pwMaxSet.put(Bunch(mes="max", val=totalChromPeaks))

            grp = 1
            if writePDF: pdf = canvas.Canvas(file + "_0.pdf")
            if writePDF: _writeFirstPage(pdf, groupSizePPM, maxTimeDeviation, align, nPolynom)
            xy = 1

            curNum = 1

            xCountshelp = []
            a=xCounts.replace(" ", "").replace(";", ",").split(",")
            for j in a:
                if "-" in j:
                    xCountshelp.extend(range(int(j[0:j.find("-")]), int(j[j.find("-")+1:len(j)])+1))
                else:
                    xCountshelp.append(int(j))
            xCounts=sorted(list(set(xCountshelp)))

            # prefetch number of bracketing steps for the progress bar
            totalThingsToDo=0
            for ionMode in ionModes:
                for tracer in tracersDeltaMZ:
                    for xCount in xCounts:
                        for cLoading in range(maxLoading, 0, -1):
                            totalThingsToDo+=1

            doneSoFar=0
            doneSoFarPercent=0


            # bracket data
            # process ionModes, tracers, xCount and loading separately
            for ionMode in ionModes:
                scanEvent = ionModes[ionMode]
                for xCount in xCounts:
                    for cLoading in range(maxLoading, 0, -1):

                        # check each processed LC-HRMS file if the used data processing parameters match
                        for res in results:
                            res.featurePairs=[]

                            if pwTextSet is not None:
                                # Log time used for bracketing
                                elapsed = (time.time() - start) / 60.
                                hours = ""
                                if elapsed >= 60.:
                                    hours = "%d hours " % (elapsed // 60)
                                mins = "%.2f mins" % (elapsed % 60.)
                                pwTextSet.put(Bunch(mes="text", val="<p align='right' >%s%s elapsed</p>\n\n\n\nProcessing \n  (Ionmode: %s, XCount: %d, Charge: %d) \n  File: %s" % (
                                                    hours, mins, ionMode, xCount, cLoading, res.fileName)))

                            for row in res.curs.execute(
                                    "SELECT c.id, c.tracer, c.NPeakCenter, c.NPeakCenterMin, c.NPeakScale, c.NSNR, c.NPeakArea, c.mz, c.xcount, c.LPeakCenter, c.LPeakCenterMin, c.LPeakScale, c.LSNR, c.LPeakArea, c.Loading, (SELECT f.fGroupId FROM featureGroupFeatures f WHERE f.fID=c.id) AS fGroupId, (SELECT t.name FROM tracerConfiguration t WHERE t.id=c.tracer) AS tracerName, c.ionMode, c.correlationsToOthers, c.NPeakAbundance, c.LPeakAbundance "
                                    "FROM chromPeaks c WHERE c.ionMode='%s' AND c.xcount=%d and c.Loading=%d"%(ionMode, xCount, cLoading)):
                                try:
                                    cp = ChromPeakPair(id=row[0], tracer=row[1], NPeakCenter=row[2], NPeakCenterMin=row[3],
                                                   NPeakScale=row[4], NSNR=row[5], NPeakArea=row[6], mz=row[7], xCount=row[8],
                                                   LPeakCenter=row[9], LPeakCenterMin=row[10], LPeakScale=row[11], LSNR=row[12],
                                                   LPeakArea=row[13], loading=row[14], fGroupID=row[15], tracerName=row[16],
                                                   ionMode=str(row[17]), correlationsToOthers=loads(base64.b64decode(str(row[18]))),
                                                   NPeakAbundance=float(row[19]), LPeakAbundance=float(row[20]))

                                    assert cp.ionMode in ionModes.keys()
                                    res.featurePairs.append(cp)
                                except TypeError as err:
                                    print "  TypeError in file %s, id %s, (Ionmode: %s, XCount: %d, Charge: %d)"%(str(res.fileName), str(row[0]), ionMode, xCount, cLoading), err.message
                                except:
                                    print "  some general error in file %s, id %s, (Ionmode: %s, XCount: %d, Charge: %d)"%(str(res.fileName), str(row[0]), ionMode, xCount, cLoading)

                            totalChromPeaks = totalChromPeaks + len(res.featurePairs)



                        if pwTextSet is not None:
                            # Log time used for bracketing
                            elapsed = (time.time() - start) / 60.
                            hours = ""
                            if elapsed >= 60.:
                                hours = "%d hours " % (elapsed // 60)
                            mins = "%.2f mins" % (elapsed % 60.)
                            pwTextSet.put(Bunch(mes="text", val="<p align='right' >%s%s elapsed</p>\n\n\n\nClustering \n  (Ionmode: %s, XCount: %d, Charge: %d)" % (
                                                hours, mins, ionMode, xCount, cLoading)))


                            totalChromPeaks = totalChromPeaks + len(res.featurePairs)


                        for tracer in tracersDeltaMZ:
                            ## TODO this does not work correctly. improve
                            doneSoFar+=1
                            if floor(doneSoFar/totalThingsToDo*100)>doneSoFarPercent:
                                doneSoFarPercent=floor(doneSoFar/totalThingsToDo*100)
                                #print "\r   %d%% performed"%doneSoFarPercent,

                            # get all results that match current bracketing criteria
                            chromPeaksAct = 0
                            allmz = []
                            for res in results:
                                chromPeaksAct = chromPeaksAct + len([chromPeak.mz for chromPeak in res.featurePairs if chromPeak.xCount == xCount and chromPeak.loading == cLoading and chromPeak.tracerName == tracer and chromPeak.ionMode == ionMode])

                            allMZs=[]
                            for res in results:
                                allMZs.extend([chromPeak.mz for chromPeak in res.featurePairs if
                                                               chromPeak.xCount == xCount and chromPeak.loading == cLoading and chromPeak.tracerName == tracer and chromPeak.ionMode == ionMode])

                            # cluster all current results with HCA
                            if len(allMZs)>0:
                                hc = HierarchicalClustering(allMZs,
                                                        dist=lambda x, y: x.getValue() - y.getValue(),
                                                        val=lambda x: x, mean=lambda x, y: x / y,
                                                        add=lambda x, y: x + y)

                                # cut HCA tree in subclusters
                                for n in cutTreeSized(hc.getTree(), groupSizePPM):
                                    try:
                                        lowMZ=min([kid.getValue() for kid in n.getKids()])
                                        minMZ=max([kid.getValue() for kid in n.getKids()])

                                        maxMZInGroup = lowMZ

                                        partChromPeaks = {}
                                        partXICs = {}
                                        toDel = {}
                                        allChromPeaks=[]
                                        for res in results:
                                            chromPeaks = []
                                            toDel[res.filePath] = set()
                                            for i in range(len(res.featurePairs)):
                                                chromPeak = res.featurePairs[i]
                                                if chromPeak.xCount == xCount and chromPeak.loading == cLoading and chromPeak.tracerName == tracer and chromPeak.ionMode == ionMode and chromPeak.mz <= minMZ:
                                                    toDel[res.filePath].add(i)
                                                    chromPeaks.append(chromPeak)
                                                    allChromPeaks.append(chromPeak)
                                                    maxMZInGroup = max(maxMZInGroup, chromPeak.mz)

                                            xics = []
                                            # get EICs of current subCluster
                                            for chromPeak in chromPeaks:

                                                res.curs.execute("SELECT x.xicL, x.times, x.scanCount, eicid FROM XICs x, chromPeaks c WHERE c.id=%d AND c.eicID=x.id" % chromPeak.id)
                                                i = 0
                                                for row in res.curs:
                                                    scanCount = float(row[2])
                                                    times = [float(r) for r in row[1].split(";")]
                                                    xicL = [float(r) for r in row[0].split(";")]
                                                    eicID = int(row[3])
                                                    i = i + 1
                                                assert (i == 1)
                                                xics.append([scanCount, times, xicL, eicID])
                                            partChromPeaks[res.filePath] = chromPeaks
                                            partXICs[res.filePath] = xics


                                        xy = xy + 1
                                        if (xy % 50) == 0:
                                            if writePDF:
                                                pdf.save()
                                                pdf = canvas.Canvas(file + "_%d.pdf" % xy)
                                                _writeFirstPage(pdf, groupSizePPM, maxTimeDeviation, align, nPolynom)
                                                print "\r   new pdf %d metabolic features.." % xy,

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
                                            pdf.drawString(40, 810, "Tracer: %s IonMode: %s" % (tracer, ionMode))
                                            pdf.drawString(40, 795, "MZ: (min: %f - max: %f (tolerated: %f ppm: %.1f)) " % (lowMZ, maxMZInGroup, minMZ, (maxMZInGroup - lowMZ) * 1000000 / lowMZ))
                                            pdf.drawString(40, 780, "xCount: %d Z: %d Feature pairs: %d" % (xCount, cLoading, len(partXICs)))
                                            pdj = []
                                            pdjm = []
                                            for r in partXICs.keys():
                                                xics = partXICs[r]
                                                for xic in xics:
                                                    pdj.append(xic[2])
                                                    pdjm.append(xic[1])
                                                    dd.append([(xic[1][i] / 60., xic[2][i]) for i in range(0, len(xic[1]))])
                                                    tr = tr + 1
                                            pdf.drawString(40, 765, "Aligning %d EICs" % (tr))

                                            lp.data = dd

                                            i = 0
                                            for k in partXICs.keys():
                                                lp.lines[i].strokeColor = colos[i % len(colos)]
                                                lp.lines[i].strokeWidth = 0.01
                                                i = i + 1

                                            lp.joinedLines = 1
                                            drawing.add(lp)
                                            renderPDF.draw(drawing, pdf, 10, 380)

                                            alignedEICs = xicAlign.getAligendXICs(pdj, pdjm, align=align, nPolynom=nPolynom)[0]
                                            drawing = Drawing(500, 350)
                                            lp = LinePlot()
                                            lp.x = 50
                                            lp.y = 30
                                            lp.height = 350
                                            lp.width = 500

                                            dd = []

                                            tr = 0
                                            for xic in alignedEICs:
                                                dd.append([(i / 60., xic[i]) for i in range(0, len(xic))])
                                                tr = tr + 1

                                            lp.data = dd

                                            i = 0
                                            for k in range(len(alignedEICs)):
                                                lp.lines[i].strokeColor = colos[i % len(colos)]
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
                                            for i in range(len(partChromPeaks[k])):
                                                do.append(k)
                                                dd.append(partChromPeaks[k][i].mz)
                                                eics.append(partXICs[k][i][2])
                                                peaks.append([partChromPeaks[k][i]])
                                                scantimes.append(partXICs[k][i][1])

                                        # optional: align EICs; get bracketed chromatographic peaks
                                        aligned = xicAlign.alignXIC(eics, peaks, scantimes, align=align, maxTimeDiff=maxTimeDeviation, nPolynom=nPolynom)
                                        aligned = [(x[0][0], int(x[0][1])) for x in aligned]

                                        if writePDF:
                                            u = 800
                                            o = 0
                                            for peak in peaks:
                                                for p in peak:
                                                    pdf.drawString(40, u, "%35s: %d %.2f -> %d %d" % (
                                                                   do[o][(1 + do[o].rfind("/")):], int(p.NPeakCenter),
                                                                   p.NPeakCenterMin / 60., int(aligned[o][0]), aligned[o][1]))
                                                u -= 20
                                                o += 1

                                        maxGroup = max(aligned, key=itemgetter(1))[1]
                                        minGroup = min(aligned, key=itemgetter(1))[1]

                                        groupedChromPeaks = []
                                        groupedChromPeaksAVGMz = []
                                        groupedChromPeaksAVGTimes = []

                                        for i in range(maxGroup + 1):
                                            groupedChromPeaks.append({})
                                            groupedChromPeaksAVGMz.append([])
                                            groupedChromPeaksAVGTimes.append([])

                                        # calculate average values (RT and mz) for feature pairs in the sub-subclusters
                                        j = 0
                                        for k in partChromPeaks.keys():
                                            for i in range(len(partChromPeaks[k])):
                                                if not (groupedChromPeaks[aligned[j][1]].has_key(k)):
                                                    groupedChromPeaks[aligned[j][1]][k] = []
                                                groupedChromPeaks[aligned[j][1]][k].append((aligned[j], partChromPeaks[k][i]))
                                                groupedChromPeaksAVGMz[aligned[j][1]].append(partChromPeaks[k][i].mz)
                                                groupedChromPeaksAVGTimes[aligned[j][1]].append(partChromPeaks[k][i].NPeakCenterMin)

                                                j = j + 1

                                        assert (j == len(aligned))
                                        assert (len(groupedChromPeaks) == len(groupedChromPeaksAVGMz) == len(groupedChromPeaksAVGTimes))

                                        # write results to data matrix and SQLite DB
                                        for i in range(minGroup, maxGroup + 1):
                                            if len(groupedChromPeaks[i]) > 0:
                                                avgmz = sum(groupedChromPeaksAVGMz[i]) / len(groupedChromPeaksAVGMz[i])
                                                avgtime = sum(groupedChromPeaksAVGTimes[i]) / len(groupedChromPeaksAVGTimes[i])

                                                f.write(str(curNum))
                                                f.write("\t")
                                                f.write("")
                                                f.write("\t")
                                                f.write(str(avgmz))
                                                f.write("\t")
                                                f.write(str(avgmz + xCount * tracersDeltaMZ[tracer] / cLoading))
                                                f.write("\t")
                                                f.write(str(xCount * tracersDeltaMZ[tracer] / cLoading))
                                                f.write("\t")
                                                f.write("%.6f - %.6f"%(min(groupedChromPeaksAVGMz[i]), max(groupedChromPeaksAVGMz[i])))
                                                f.write("\t")
                                                f.write("%.2f" % (avgtime / 60.))
                                                f.write("\t")
                                                f.write("%.3f - %.3f" %(min(groupedChromPeaksAVGTimes[i])/60., max(groupedChromPeaksAVGTimes[i])/60.))
                                                f.write("\t")
                                                f.write("%d" % xCount)
                                                f.write("\t")
                                                f.write(str(cLoading))
                                                f.write("\t")
                                                f.write(scanEvent)
                                                f.write("\t")
                                                f.write(ionMode)
                                                f.write("\t")
                                                f.write(str(tracer))

                                                f.write("\t")
                                                f.write("") # OGroup
                                                f.write("\t")
                                                f.write("") # Ion
                                                f.write("\t")
                                                f.write("") # Loss
                                                f.write("\t")
                                                f.write("") # M


                                                SQLInsert(resDB.curs, "GroupResults", id=curNum, mz=avgmz, lmz=avgmz + xCount * tracersDeltaMZ[tracer] / cLoading,
                                                            dmz=xCount * tracersDeltaMZ[tracer] / cLoading, rt=avgtime, xn=xCount, charge=cLoading, ionisationMode=ionMode, scanEvent=scanEvent,
                                                            tracer=str(tracer))

                                                doublePeak = 0
                                                toapp = ""
                                                for j in range(len(results)):
                                                    res = results[j].filePath
                                                    f.write("\t")
                                                    toapp = toapp + "\t"
                                                    if groupedChromPeaks[i].has_key(res) and len(groupedChromPeaks[i][res]) > 0:

                                                        for peak in groupedChromPeaks[i][res]:
                                                            SQLInsert(resDB.curs, "FoundFeaturePairs", file=results[j].fileName, featurePairID=peak[1].id, featureGroupID=peak[1].fGroupID, resID=curNum, areaN=peak[1].NPeakArea, areaL=peak[1].LPeakArea, featureType="foundPattern")

                                                        appe = False
                                                        for peak in groupedChromPeaks[i][res]:
                                                            if appe:
                                                                toapp = toapp + ";"
                                                                doublePeak = doublePeak + 1
                                                            else:
                                                                appe = True
                                                            toapp = toapp + str(peak[1].NPeakArea / peak[1].LPeakArea)

                                                        # Area of native, monoisotopic isotopolog
                                                        appe = False
                                                        for peak in groupedChromPeaks[i][res]:
                                                            if appe:
                                                                f.write(";")
                                                            else:
                                                                appe = True
                                                            f.write(str(peak[1].NPeakArea))

                                                        f.write("\t")
                                                        # Area of labeled isotopolog
                                                        appe = False
                                                        for peak in groupedChromPeaks[i][res]:
                                                            if appe:
                                                                f.write(";")
                                                            else:
                                                                appe = True
                                                            f.write(str(peak[1].LPeakArea))

                                                        f.write("\t")
                                                        # Abundance of native, monoisotopic isotopolog
                                                        appe = False
                                                        for peak in groupedChromPeaks[i][res]:
                                                            if appe:
                                                                f.write(";")
                                                            else:
                                                                appe = True
                                                            f.write(str(peak[1].NPeakAbundance))

                                                        f.write("\t")
                                                        # Abundance of labeled isotopolog
                                                        appe = False
                                                        for peak in groupedChromPeaks[i][res]:
                                                            if appe:
                                                                f.write(";")
                                                            else:
                                                                appe = True
                                                            f.write(str(peak[1].LPeakAbundance))

                                                        f.write("\t")
                                                        # Feature ID
                                                        appe = False
                                                        for peak in groupedChromPeaks[i][res]:
                                                            if appe:
                                                                f.write(";")
                                                            else:
                                                                appe = True
                                                            f.write("%d" % peak[1].id)

                                                        f.write("\t")
                                                        # Group ID
                                                        appe = False
                                                        for peak in groupedChromPeaks[i][res]:
                                                            if appe:
                                                                f.write(";")
                                                            else:
                                                                appe = True
                                                            f.write("%d" % peak[1].fGroupID)

                                                        f.write("\t")
                                                        # CorrelatesTo
                                                        appe = False
                                                        for peak in groupedChromPeaks[i][res]:
                                                            if appe:
                                                                f.write(";")
                                                            else:
                                                                appe = True
                                                            f.write("%s" % str(peak[1].correlationsToOthers))

                                                    else:
                                                        f.write("\t\t\t\t\t\t")
                                                f.write(toapp)
                                                f.write("\t%d" % doublePeak)
                                                f.write("\n")

                                                curNum = curNum + 1

                                        if writePDF: pdf.showPage()

                                    except Exception as e:
                                        import traceback
                                        traceback.print_exc()
                                        logging.error(str(traceback))

                                        logging.error("Error", str(e))
                                    grp = grp + 1

                                    for res in results:
                                        if toDel.has_key(res.filePath):
                                            toDel[res.filePath] = [a for a in toDel[res.filePath]]
                                            toDel[res.filePath].sort()
                                            toDel[res.filePath].reverse()
                                            for i in toDel[res.filePath]:
                                                res.featurePairs.pop(i)
                                                chromPeaksAct = chromPeaksAct - 1
                                                totalChromPeaksProc = totalChromPeaksProc + 1

                                    if pwValSet is not None: pwValSet.put(Bunch(mes="value",val=totalChromPeaksProc))
                                    if pwTextSet is not None:

                                        # Log time used for bracketing
                                        elapsed = (time.time() - start) / 60.
                                        hours = ""
                                        if elapsed >= 60.:
                                            hours = "%d hours " % (elapsed // 60)
                                        mins = "%.2f mins" % (elapsed % 60.)
                                        pwTextSet.put(Bunch(mes="text",val="<p align='right' >%s%s elapsed</p>\n\n\n\nBracketing results\n%d feature pairs (%d individual pairs processed).."%(hours, mins, curNum, totalChromPeaksProc)))

            if writePDF: pdf.save()


            import uuid
            import platform
            import datetime

            identifier="%s_%s_%s"%(str(uuid.uuid1()), str(platform.node()), str(datetime.datetime.now()))

            f.write("## MetExtract II %s\n"%(Bunch(MetExtractVersion=meVersion, RVersion=rVersion, UUID_ext=identifier).dumpAsJSon().replace("\"", "'")))
            f.write("## Individual files processing parameters %s\n"%(generalProcessingParams.dumpAsJSon().replace("\"", "'")))
            processingParams=Bunch()
            processingParams.FPBracketing_xCounts=xCounts
            processingParams.FPBracketing_groupSizePPM=groupSizePPM
            processingParams.FPBracketing_positiveScanEvent=positiveScanEvent
            processingParams.FPBracketing_negativeScanEvent=negativeScanEvent
            processingParams.FPBracketing_maxTimeDeviation=maxTimeDeviation
            processingParams.FPBracketing_maxLoading=maxLoading
            processingParams.FPBracketing_resultsFile=file
            processingParams.FPBracketing_align=align
            if align: processingParams.FPBracketing_nPolynom=nPolynom
            f.write("## Bracketing files processing parameters %s\n"%(processingParams.dumpAsJSon().replace("\"", "'")))

        exportAsFeatureML.convertMEMatrixToFeatureMLSepPolarities(file)

    except Exception as ex:
        import traceback
        traceback.print_exc()
        print(str(traceback))

        logging.error("Error during bracketing of files")

    finally:
        resDB.conn.commit()
        resDB.curs.close()
        resDB.conn.close()






























########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################














# store used configuration to DB file
def writeMetaboliteGroupingConfigToDB(curs, minConnectionsInFiles, minConnectionRate, groups):
    SQLInsert(curs, "config", key="MEGROUP_groups", value=str(groups))
    SQLInsert(curs, "config", key="MEGROUP_minConnectionsInFiles", value=str(minConnectionsInFiles))
    SQLInsert(curs, "config", key="MEGROUP_minConnectionRate", value=str(minConnectionRate))




def calculateMetaboliteGroups(file="./results.tsv", groups=[],
                              minConnectionsInFiles=1, minConnectionRate=0.4,
                              runIdentificationInstance=None):

    resDB=Bunch(conn=None, curs=None)
    resDB.conn=connect(file+getDBSuffix())
    resDB.curs=resDB.conn.cursor()

    try:
        # Select only those groups that shall be used for metabolite group grouping
        useGroups=[]
        useGroupsForConfig=[]
        for group in groups:
            if group.useForMetaboliteGrouping:
                useGroups.append(group)
                useGroupsForConfig.append(str(group.name)+":"+str(group.files))

        # connect to results db

        writeMetaboliteGroupingConfigToDB(resDB.curs, minConnectionsInFiles, minConnectionRate, str(useGroupsForConfig).replace("'", "").replace("\"", ""))


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

        doublePeaks=0
        for row in table.getData(cols=["COUNT(*)"], where="doublePeak>0"):
            doublePeaks=int(row)

        if doublePeaks>0:
            print "found double peaks:", doublePeaks
        writeCols=[]
        for column in table.getColumns():
            if not column.name.endswith("_CorrelatesTo"):
                writeCols.append(column.name)
        TableUtils.saveFile(table, file.replace(".tsv", ".doublePeaks.tsv"), cols=writeCols, order="OGroup, MZ, Xn", where="doublePeak>0")

        table.executeSQL("DELETE FROM :table: WHERE doublePeak>0")
        ## file specific columns

        for group in useGroups:
            for f in group.files:
                fShort=f[f.rfind("/")+1 : f.rfind(".")]
                if fShort[0].isdigit():
                    fShort="_"+fShort

                assert "%s_FID"%fShort in cols
                assert "%s_GroupID"%fShort in cols


        # fetch all feature pairs from the results file
        nodes={}
        nums=[featurePairNum for featurePairNum in table.getData(cols=["Num"])]
        for featurePairNum in nums:
            nodes[featurePairNum]={}

        # get overall feature pair clusters

        # add 1 to the connection density if a connection was detected in a file
        totalFilesToUse=0
        done=set()
        for group in useGroups:
            for f in group.files:
                if f not in done:
                    totalFilesToUse+=1
                    done.add(f)
                    fShort=f[f.rfind("/")+1 : f.rfind(".")]
                    if fShort[0].isdigit():
                        fShort="_"+fShort


                    fileNumMapping={}
                    for featurePairNum, featurePairNumInFile in table.getData(cols=["Num", "%s_FID"%fShort]):
                        if featurePairNumInFile!="":
                            featurePairNumInFile=float(featurePairNumInFile)
                            fileNumMapping[featurePairNumInFile]=featurePairNum
                    for featurePairNum, groupID, correlationsToOthers in table.getData(cols=["Num", "%s_GroupID"%fShort, "%s_CorrelatesTo"%fShort]):

                        if groupID!="":
                            correlationsToOthers=ast.literal_eval(str(correlationsToOthers))

                            for corWith in correlationsToOthers:
                                if corWith in fileNumMapping.keys():
                                    a=featurePairNum
                                    b=fileNumMapping[corWith]
                                    if a not in nodes.keys():
                                        nodes[a]={}
                                    if b not in nodes[a].keys():
                                        nodes[a][b]=Bunch(posConns=0, negConns=0)
                                    nodes[a][b].posConns+=1

                                    b=featurePairNum
                                    a=fileNumMapping[corWith]
                                    if a not in nodes.keys():
                                        nodes[a]={}
                                    if b not in nodes[a].keys():
                                        nodes[a][b]=Bunch(posConns=0, negConns=0)
                                    nodes[a][b].posConns+=1


        toRemove=set()
        for peakA in nodes.keys():
            for peakB in nodes[peakA].keys():
                node=nodes[peakA][peakB]
                if node.negConns==0 and node.posConns==0:
                    if peakA<peakB:
                        toRemove.add((peakA, peakB))
                    else:
                        toRemove.add((peakB, peakA))

        for peakA, peakB in toRemove:
            del nodes[peakA][peakB]
            if peakB!=peakA:
                del nodes[peakB][peakA]


        # remove edges with no connection
        toRemove=set()
        for peakA in nodes.keys():
            for peakB in nodes[peakA].keys():
                conns=nodes[peakA][peakB]
                if conns.posConns>=minConnectionsInFiles and ( conns.negConns==0 or ( conns.negConns>0 and (1.*conns.posConns/(conns.posConns+conns.negConns))>=minConnectionRate )):
                    pass
                else:
                    if peakA<peakB:
                        toRemove.add((peakA, peakB))
                    else:
                        toRemove.add((peakB, peakA))

        for peakA, peakB in toRemove:
            del nodes[peakA][peakB]
            if peakB!=peakA:
                del nodes[peakB][peakA]



        # Separate feature pair clusters; crude
        tGroups = getSubGraphsFromDictDict(nodes)

        # Separate feature pair clusters; softer
        curGroup=1
        for tGroup in tGroups:
            if len(tGroup)==1:
                # if only one feature pair is in the group, do nothing
                table.setData(cols=["OGroup"], vals=[curGroup], where="Num=%d"%tGroup[0])
                resDB.curs.execute("UPDATE GroupResults SET OGroup=%d WHERE id = %d"%(curGroup, tGroup[0]))
            else:
                table.setData(cols=["OGroup"], vals=[curGroup], where="Num IN (%s)"%(",".join([str(t) for t in tGroup])))
                resDB.curs.execute("UPDATE GroupResults SET OGroup=%d WHERE id IN (%s)"%(curGroup, ",".join([str(t) for t in tGroup])))

            curGroup+=1

        '''
        # 2nd separation
        for oGroupID in range(curGroup):
            nodes={}
            nums=[featurePairNum for featurePairNum in table.getData(cols=["Num"], where="OGroup=%d"%oGroupID)]
            if len(nums)>1:
                for featurePairNum in nums:
                    nodes[featurePairNum]={}
                    for featurePairNum2 in nums:
                        nodes[featurePairNum][featurePairNum2]=Bunch(posConns=0, negConns=0)

                # add 1 to the connection density if a connection was detected in a file
                totalFilesToUse=0
                done=set()
                for group in useGroups:
                    for f in group.files:
                        if f not in done:
                            totalFilesToUse+=1
                            done.add(f)
                            fShort=f[f.rfind("/")+1 : f.rfind(".")]
                            if fShort[0].isdigit():
                                fShort="_"+fShort

                            # get feature pair clusters from current file
                            fileGroupIDs=defaultdict(list)

                            fileNumMapping={}
                            for featurePairNum, featurePairNumInFile in table.getData(cols=["Num", "%s_FID"%fShort], where="OGroup=%d"%oGroupID):
                                if featurePairNumInFile!="":
                                    fileNumMapping[featurePairNumInFile]=featurePairNum

                            for featurePairNum, groupID, correlationsToOthers in table.getData(cols=["Num", "%s_GroupID"%fShort, "%s_CorrelatesTo"%fShort], where="OGroup=%d"%oGroupID):
                                if groupID!="":
                                    correlationsToOthers=ast.literal_eval(str(correlationsToOthers))
                                    fileGroupIDs[groupID].append(featurePairNum)

                                    for corWith in correlationsToOthers:
                                        nodes[featurePairNum][fileNumMapping[corWith]].posConns+=1
                print nums
        '''


        #tGroups=groups



        # Annotate groups with common adducts and in-source fragments
        if runIdentificationInstance is not None:
            groups=defaultdict(list)
            for row in table.getData(cols=["Num", "OGroup", "MZ", "Ionisation_Mode", "Charge", "Ion", "Xn", "Loss", "M"]):
                num, ogrp, mz, ionMode, loading, adducts, xCount, fDesc, ms=row
                groups[ogrp].append(ChromPeakPair(id=num, fGroupID=ogrp, mz=mz, ionMode=ionMode, loading=loading, adducts=[], heteroAtomsFeaturePairs=[],
                                                  xCount=xCount, fDesc=[]))

            for ogrp in groups.keys():
                chromPeaks={}
                for fp in groups[ogrp]:
                    chromPeaks[fp.id]=fp

                runIdentificationInstance.annotateChromPeaks(chromPeaks.keys(), chromPeaks)

            for ogrp in groups.keys():
                for fp in groups[ogrp]:
                    table.setData(cols=["Ion", "Loss", "M"], vals=[",".join([str(a) for a in fp.adducts]),
                                                              "",
                                                              ",".join([str(a) for a in fp.Ms])],
                                  where="Num = %d"%(fp.id))
                    resDB.curs.execute("UPDATE GroupResults SET Ion='%s', Loss='%s', M='%s' WHERE id = %d"%(",".join([str(a) for a in fp.adducts]), "", ",".join([str(a) for a in fp.Ms]), fp.id))

        resDB.curs.execute("DELETE FROM GroupResults WHERE id NOT IN (%s)"%(",".join([str(num) for num in table.getData(cols=["Num"])])))


        ## reassign feature group ids
        resDB.curs.execute("UPDATE GroupResults SET OGroup='X'||OGroup")
        table.executeSQL("UPDATE :table: SET OGroup='X'||OGroup")
        oGrps=[]
        curGrp=1
        curFP=1
        for row in resDB.curs.execute("SELECT OGroup, AVG(rt) FROM GroupResults GROUP BY OGroup ORDER BY AVG(rt)"):
            oGrps.append(str(row[0]))
        for tgrp in oGrps:
            resDB.curs.execute("UPDATE GroupResults SET OGroup='%d' WHERE OGroup='%s'"%(curGrp, tgrp))
            table.setData(cols=["OGroup"], vals=[curGrp], where="OGroup='%s'"%(tgrp))
            curGrp=curGrp+1




        processingParams=Bunch()
        processingParams.MEConvoluting_groups=str(useGroupsForConfig).replace("'","").replace("\"","")
        processingParams.MEConvoluting_connThreshold=minConnectionsInFiles
        table.addComment("## Convolution FPs processing parameters %s"%(processingParams.dumpAsJSon().replace("\"", "'")))

        writeCols=[]
        for column in table.getColumns():
            if not column.name.endswith("_CorrelatesTo"):
                writeCols.append(column.name)

        TableUtils.saveFile(table, file, cols=writeCols, order="OGroup, MZ, Xn")

    except Exception as ex:
        import traceback
        traceback.print_exc()

        raise ex

    finally:
        resDB.conn.commit()
        resDB.curs.close()
        resDB.conn.close()
