from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.charts.lineplots import ScatterPlot
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.widgets.markers import makeMarker
from reportlab.lib import pagesizes
from reportlab.lib.colors import Color, HexColor
from reportlab.lib.units import mm
from reportlab.platypus import Paragraph, Table
from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF
from reportlab.lib.styles import getSampleStyleSheet

import sqlite3

import matplotlib


from math import sqrt
from os.path import join
import sys
sys.path.append("D:/PyMetExtract/PyMetExtract")
import pickle
import base64
from TableUtils import TableUtils
from utils import smoothDataSeries, Bunch
from Chromatogram import Chromatogram

from DesignPatterns.SingletonDecorator import Singleton

def mean(x):
    if len(x)==0: return 0
    return sum(x)/len(x)
def std(x):
    if len(x)<2: return 0
    me=mean(x)

    return sqrt(sum([pow(y-me, 2) for y in x])/len(x))
def rstd(x):
    if len(x)<2: return 0
    st=std(x)
    me=mean(x)

    return st/me

def scientificFormatter(val):
    return "%.0e"%val
def noLabel(a):
    return ""


class ResDB:
    def __init__(self, file):
        self.conn = sqlite3.connect(file)
        self.curs = self.conn.cursor()

    def close(self):
        self.curs.close()
        self.conn.close()





























########################################################################################################################
########################################################################################################################
########################################################################################################################


def generatePageFor(feature, mzXMLs, experimentalGroups, pdf, ppm=5., plotLabeld=True, normEICs=False, showGroupRTShift=True):
    totalFiles=sum([len(e.files) for e in experimentalGroups])
    filesDone=0

    pdf.drawString(20, 820, "Num: %s, OGroup: %s,  MZ: %.5f, LMZ: %.5f, Xn: %d,  RT: %.2f, ScanEvent: %s, Ionisation mode: %s" % (feature.num, feature.ogroup, feature.mz, feature.lmz, feature.xn, feature.rt, feature.filterLine, feature.ionisationMode))
    pdf.drawString(20, 800, "Comment: %s"%(feature.comment))

    maxEICValue=0

    ## EICs
    eicsDD = []
    eicsSeparatedDD = []
    eicsCols = []
    eicsLDD = []
    eicsLSeparatedDD = []
    eicsLCols = []
    scanDD = []
    scanCols = []
    scanWidths = []
    scanBackgroundDD = []
    scanBackgroundCols = []
    scanBackgroundWidths = []
    AreaDD = []
    AreaCols = []
    AreaLDD = []
    AreaLCols = []


    AreaDD.append([((feature.rt - feature.rtBorders) -.1, -2.75 * ppm)])
    AreaCols.append("Black")
    AreaDD.append([((feature.rt - feature.rtBorders) -.1, 2.75 * ppm)])
    AreaCols.append("Black")
    AreaDD.append([((feature.rt + feature.rtBorders) +.1, 2.75 * ppm)])
    AreaCols.append("Black")
    AreaDD.append([((feature.rt + feature.rtBorders) +.1, -2.75 * ppm)])
    AreaCols.append("Black")

    AreaLDD.append([((feature.rt - feature.rtBorders) -.1, -2.75 * ppm)])
    AreaLCols.append("Black")
    AreaLDD.append([((feature.rt - feature.rtBorders) -.1, 2.75 * ppm)])
    AreaLCols.append("Black")
    AreaLDD.append([((feature.rt + feature.rtBorders) +.1, 2.75 * ppm)])
    AreaLCols.append("Black")
    AreaLDD.append([((feature.rt + feature.rtBorders) +.1, -2.75 * ppm)])
    AreaLCols.append("Black")

    for i, expGroup in enumerate(experimentalGroups):
        for file in expGroup.files:

            useFilterLine=""
            for fl in mzXMLs[file].getFilterLinesExtended():
                if fl in feature.filterLine:
                    useFilterLine=fl

            try:
                #print(("\rProcessing.. |%%-%ds| (%%d/%%d)"%(totalFiles)) % ("*"*filesDone, filesDone, totalFiles),)
                filesDone+=1

                ## Native EICs
                eic, times, scanIds, mzs = mzXMLs[file].getEIC(mz=feature.mz, ppm=ppm, filterLine=useFilterLine,
                                                              removeSingles=False, intThreshold=0, useMS1=True)

                if normEICs:
                    mval = max(maxEICValue, max([eic[j] for j in range(len(eic)) if
                                                        (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)]))
                    eic=[e/mval for e in eic]

                maxEICValue=max(maxEICValue, max([eic[j] for j in range(len(eic)) if
                               (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)]))

                ## plot overlaied EICs to drawing
                eicsDD.append([(times[j] / 60., eic[j]) for j in range(len(eic)) if
                               (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)])

                ## plot separated EICs to drawing
                eicsSeparatedDD.append([((times[j] + showGroupRTShift*(i+3) * 60) / 60. - feature.rt, eic[j]) for j in range(len(eic)) if
                                        (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)])

                eicsCols.append(expGroup.color)

            except Exception as ex:
                print("-- skipping for", file, "(message: %s)" % str(ex))


            try:
                ## Labeled EICs
                eic, times, scanIds, mzsL = mzXMLs[file].getEIC(mz=feature.lmz, ppm=ppm, filterLine=useFilterLine,
                                                              removeSingles=False, intThreshold=0, useMS1=True)

                if normEICs:
                    mval = max(maxEICValue, max([eic[j] for j in range(len(eic)) if
                                                        (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)]))
                    eic=[e/mval for e in eic]

                maxEICValue=max(maxEICValue, max([eic[j] for j in range(len(eic)) if
                               (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)]))

                ## plot overlaied EICs to drawing
                eicsLDD.append([(times[j] / 60., eic[j]) for j in range(len(eic)) if
                               (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)])

                ## plot separated EICs to drawing
                eicsLSeparatedDD.append([((times[j] + showGroupRTShift*(i+3) * 60) / 60. - feature.rt, eic[j]) for j in range(len(eic)) if
                                        (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)])

                eicsLCols.append(expGroup.color)

            except Exception as ex:
                print("-- skipping for", file, "(message: %s)" % str(ex))

            try:
                ## find best scan in EIC peak
                bestScanIndex=max([(j, ab) for j, ab in enumerate(eic) if
                                           (feature.rt - feature.rtBorders/3) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders/3)],
                                  key=lambda x: x[1])[0]

                scan=mzXMLs[file].getIthMS1Scan(index=bestScanIndex, filterLine=useFilterLine)
                uInds=[j for j, mz in enumerate(scan.mz_list) if
                         feature.mz-3 <= mz <= feature.lmz+3]
                f=[]
                for uInd in uInds:
                    mz=scan.mz_list[uInd]
                    ab=scan.intensity_list[uInd]

                    if abs(mz-feature.mz)*1E6/feature.mz <= 5. or \
                       abs(mz-feature.lmz)*1E6/feature.lmz <= 5.:
                        scanDD.append([(mz, 0), (mz, ab)])
                        scanCols.append(matplotlib.colors.cnames[expGroup.color.lower()])
                        scanWidths.append(1.)
                    elif abs(mz-feature.mz-1.00335484)*1E6/feature.mz <= 5. or \
                         abs(mz+1.00335484-feature.lmz)*1E6/feature.lmz <= 5.:
                        scanDD.append([(mz, 0), (mz, ab)])
                        scanCols.append(matplotlib.colors.cnames[expGroup.color.lower()])
                        scanWidths.append(.5)
                    else:
                        scanBackgroundDD.append([(mz, 0), (mz, ab)])
                        scanBackgroundCols.append("#396FAC")
                        scanBackgroundWidths.append(.1)
            except Exception as ex:
                print("-- skipping for", file, "(message: %s)"%str(ex))


            ## get chromatogram area
            scans, times, scanIDs= \
                mzXMLs[file].getSpecificArea((feature.rt - feature.rtBorders)*60,
                                             (feature.rt + feature.rtBorders)*60,
                                             feature.mz*(1E6-ppm*2.5)/1E6,
                                             feature.mz*(1E6+ppm*2.5)/1E6,
                                             filterLine=useFilterLine)

            for j, scan in enumerate(scans):
                rt=times[j]
                for k in range(len(scan)):
                    mz, ab=scan[k]

                    AreaDD.append([(rt/60., (mz-feature.mz)*1E6/feature.mz)])
                    AreaCols.append(expGroup.color)


            ## get chromatogram area for lmz
            scans, times, scanIDs= \
                mzXMLs[file].getSpecificArea((feature.rt - feature.rtBorders)*60,
                                             (feature.rt + feature.rtBorders)*60,
                                             feature.lmz*(1E6-ppm*2.5)/1E6,
                                             feature.lmz*(1E6+ppm*2.5)/1E6,
                                             filterLine=useFilterLine)

            for j, scan in enumerate(scans):
                rt=times[j]
                for k in range(len(scan)):
                    mz, ab=scan[k]

                    AreaLDD.append([(rt/60., (mz-feature.lmz)*1E6/feature.lmz)])
                    AreaLCols.append(expGroup.color)






    ## plot overlaied EICs
    drawingEICs = Drawing(500, 350)
    lp = LinePlot()
    lp.x = 20
    lp.y = 20
    lp.height = 230
    lp.width = 100
    lp.data = eicsDD
    lp.joinedLines = 1
    for i, col in enumerate(eicsCols):
        lp.lines[i].strokeColor = HexColor(matplotlib.colors.cnames[col.lower()])
        lp.lines[i].strokeWidth = 0.1
    lp.yValueAxis.labelTextFormat = scientificFormatter
    lp.yValueAxis.valueMax=maxEICValue
    drawingEICs.add(lp)
    renderPDF.draw(drawingEICs, pdf, 20, 510)
    pdf.drawString(60, 770, "Native EICs")

    ## plot separated EICs
    drawingEICs = Drawing(500, 500)
    lp = LinePlot()
    lp.x = 20
    lp.y = 20
    lp.height = 230
    lp.width = 400
    lp.data = eicsSeparatedDD
    lp.joinedLines = 1
    for i, col in enumerate(eicsCols):
        lp.lines[i].strokeColor = HexColor(matplotlib.colors.cnames[col.lower()])
        lp.lines[i].strokeWidth = 0.1
    lp.yValueAxis.labelTextFormat = scientificFormatter
    lp.yValueAxis.valueMax=maxEICValue
    drawingEICs.add(lp)
    renderPDF.draw(drawingEICs, pdf, 170, 510)
    pdf.drawString(240, 770, "Native EICs (artificially separated by assigned group)")

    drawingEICs = Drawing(500, 500)
    lp = LinePlot()
    lp.x = 20
    lp.y = 20
    lp.height = 230
    lp.width = 400
    lp.data = eicsSeparatedDD
    lp.joinedLines = 1
    for i, col in enumerate(eicsCols):
        lp.lines[i].strokeColor = HexColor(matplotlib.colors.cnames[col.lower()])
        lp.lines[i].strokeWidth = 0.1
    lp.yValueAxis.labelTextFormat = scientificFormatter
    drawingEICs.add(lp)
    renderPDF.draw(drawingEICs, pdf, 620, 510)
    pdf.drawString(240, 770, "Native EICs (artificially separated by assigned group)")



    if plotLabeld:
        ## plot overlaied EICs
        drawingEICs = Drawing(500, 350)
        lp = LinePlot()
        lp.x = 20
        lp.y = 20
        lp.height = 230
        lp.width = 100
        lp.data = eicsLDD
        lp.joinedLines = 1
        for i, col in enumerate(eicsLCols):
            lp.lines[i].strokeColor = HexColor(matplotlib.colors.cnames[col.lower()])
            lp.lines[i].strokeWidth = 0.1
        lp.yValueAxis.labelTextFormat = scientificFormatter
        lp.yValueAxis.valueMax=maxEICValue
        drawingEICs.add(lp)
        renderPDF.draw(drawingEICs, pdf, 20, 230)
        pdf.drawString(60, 490, "Labeled EICs")

        ## plot separated EICs
        drawingEICs = Drawing(500, 500)
        lp = LinePlot()
        lp.x = 20
        lp.y = 20
        lp.height = 230
        lp.width = 400
        lp.data = eicsLSeparatedDD
        lp.joinedLines = 1
        for i, col in enumerate(eicsLCols):
            lp.lines[i].strokeColor = HexColor(matplotlib.colors.cnames[col.lower()])
            lp.lines[i].strokeWidth = 0.1
        lp.yValueAxis.labelTextFormat = scientificFormatter
        lp.yValueAxis.valueMax=maxEICValue
        drawingEICs.add(lp)
        renderPDF.draw(drawingEICs, pdf, 170, 230)
        pdf.drawString(240, 490, "Labeled EICs (artificially separated by assigned group)")

        drawingEICs = Drawing(500, 500)
        lp = LinePlot()
        lp.x = 20
        lp.y = 20
        lp.height = 230
        lp.width = 400
        lp.data = eicsLSeparatedDD
        lp.joinedLines = 1
        for i, col in enumerate(eicsLCols):
            lp.lines[i].strokeColor = HexColor(matplotlib.colors.cnames[col.lower()])
            lp.lines[i].strokeWidth = 0.1
        lp.yValueAxis.labelTextFormat = scientificFormatter
        drawingEICs.add(lp)
        renderPDF.draw(drawingEICs, pdf, 620, 230)
        pdf.drawString(240, 490, "Labeled EICs (artificially separated by assigned group)")







    try:
        ## plot MS scans
        drawingEICs = Drawing(500, 310)
        lp = LinePlot()
        lp.x = 20
        lp.y = 20
        lp.height = 160
        lp.width = 400
        lp.data = scanBackgroundDD+scanDD
        lp.joinedLines = 1
        for i, col in enumerate(scanBackgroundCols):
            lp.lines[i].strokeColor = HexColor(col)
            lp.lines[i].strokeWidth = scanBackgroundWidths[i]
        for i, col in enumerate(scanCols):
            lp.lines[i+len(scanBackgroundCols)].strokeColor = HexColor(col)
            lp.lines[i+len(scanBackgroundCols)].strokeWidth = scanWidths[i]
        lp.yValueAxis.labelTextFormat = scientificFormatter
        drawingEICs.add(lp)
        renderPDF.draw(drawingEICs, pdf, 20, 20)
        pdf.drawString(240, 210, "MS scans (all samples)")

        # helper function for ReportLab
        def noLabel(a):
            return ""
        ## plot chromatogram area
        drawingEICs = Drawing(500, 310)
        sp = ScatterPlot()
        sp.x = 20
        sp.y = 20
        sp.height = 160
        sp.width = 250
        sp.data = AreaDD
        for i, col in enumerate(AreaCols):
            col=HexColor(matplotlib.colors.cnames[col.lower()])
            col.alpha=0.25
            sp.lines[i].strokeColor = col
        sp.lineLabelFormat=noLabel
        drawingEICs.add(sp)
        renderPDF.draw(drawingEICs, pdf, 620, 20)


        ## plot chromatogram area
        drawingEICs = Drawing(500, 310)
        sp = ScatterPlot()
        sp.x = 20
        sp.y = 20
        sp.height = 160
        sp.width = 250
        sp.data = AreaLDD
        for i, col in enumerate(AreaLCols):
            col=HexColor(matplotlib.colors.cnames[col.lower()])
            col.alpha=0.25
            sp.lines[i].strokeColor = col
        sp.lineLabelFormat=noLabel
        drawingEICs.add(sp)
        renderPDF.draw(drawingEICs, pdf, 900, 20)

    except Exception as ex:
        print(ex)

    pdf.drawString(3, 3, "Notes: Raw-data is shown, no processing has been performed, axis labels missing (bug in library), order/color of separated EICs corresponds to group table, names for groups in separated EICs cannot be set")

    data=[["", "G", "Name"]]
    style=[('FONTSIZE', (0, 0), (-1, -1), 8)]
    for i, expGroup in enumerate(experimentalGroups):
        data.append(["", str(i+1), expGroup.name])
        style.append(('BACKGROUND', (0, len(data)-1), (0, len(data)-1), expGroup.color))

    table = Table(data, style=style)
    w, h = table.wrapOn(pdf, 150, len(data))
    table.drawOn(pdf, 1050, 220)





def generatePDF(experimentalGroups, metabolitesToPlot, saveTo, ppm=5., mzXMLs=None, pw=None, drawMetaboliteIntroPage=True,
                intensityCutoff=-1, plotLabeld=True, normEICs=False, showGroupRTShift=True):

    importFilter=[]
    for metabolite in metabolitesToPlot:
        for feature in metabolite.features:
            importFilter.append((feature.mz*(1.-ppm*5/1E6), feature.mz*(1.+ppm*5/1E6)))
            importFilter.append((feature.lmz * (1. - ppm * 5 / 1E6), feature.lmz * (1. + ppm * 5 / 1E6)))

    if mzXMLs is None:
        mzXMLs={}
        print("Importing mzXMLs")
        for i, expGroup in enumerate(experimentalGroups):
            print("   Group: %s"%(expGroup.name))
            for file in expGroup.files:
                print("      File: %s"%(file))

                if file not in mzXMLs.keys():
                    mzXML = Chromatogram()
                    mzXML.parse_file(file, intensityCutoff=intensityCutoff)#, mzFilter=importFilter)
                    print("         --> FilterLines:", mzXML.getFilterLines())
                    mzXMLs[file] = mzXML



    done=0
    ## geneate PDF page for each feature
    print("\nGenerating PDF pages for features")

    pdf=None

    if pw!=None:
        pw.setMax(sum([len(metabolite.features) for metabolite in metabolitesToPlot]))
        pw.setValue(0)

    done=0
    f=0
    for metabolite in metabolitesToPlot:
        print(("\rProcessing.. |%%-%ds| (%%d/%%d)" % (len(metabolitesToPlot))) % ("*" * done, done, len(metabolitesToPlot)),)
        if not done%100:
            if pdf is not None:
                pdf.save()
            pdf = canvas.Canvas(saveTo.replace(".pdf", "_%d.pdf"%(done//100)), pagesize=pagesizes.landscape(pagesizes.A3))
        done+=1

        if drawMetaboliteIntroPage:
            pdf.drawString(20, 820, "Name: %s OGroup: %s Number of features: %d" % (metabolite.name, metabolite.ogroup, len(metabolite.features)))
            pdf.showPage()

        if pw!=None: pw.setText("   Metabolite: %s"%(metabolite.ogroup))
        for feature in metabolite.features:
            if pw!=None: pw.setText(" Metabolite: %s      Num: %s, mz: %.5f, rt: %.2f"%(metabolite.ogroup, feature.num, feature.mz, feature.rt))

            try:
                generatePageFor(feature, mzXMLs, experimentalGroups, pdf, ppm=ppm, plotLabeld=plotLabeld, normEICs=normEICs, showGroupRTShift=showGroupRTShift)
            except Exception as ex:
                print(ex.message)
                pass

            pdf.showPage()

            if pw!=None: pw.setValueu(f)
            f=f+1


    ## save PDF
    pdf.save()




























########################################################################################################################
########################################################################################################################
########################################################################################################################


















if __name__=="__main__" and True:


    ########################################################################################################################
    ######   Negative mode
    ########################################################################################################################
    if False:
        experimentalGroups=[]
        path="E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs"
        experimentalGroups.append(Bunch(name="NHC dmalQ",
                                        files=[path+"/NHC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_26_1AH_Cells_n.mzXML",
                                               path+"/NHC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_28_1BH_Cells_n.mzXML"
                                              ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="NHC dmalQDmaa",
                                        files=[path+"/NHC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_30_2AH_Cells_n.mzXML",
                                               path+"/NHC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_32_2BH_Cells_n.mzXML"
                                              ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="NHC dmnaT",
                                        files=[path+"/NHC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_10_B1H_Cells_n.mzXML",
                                               path+"/NHC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_11_B2H_Cells_n.mzXML",
                                               path+"/NHC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_12_B3H_Cells_n.mzXML"
                                               ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="NHC wildtype",
                                        files=[path+"/NHC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_7_A1H_Cells_n.mzXML",
                                               path+"/NHC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_8_A2H_Cells_n.mzXML",
                                               path+"/NHC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_9_A3H_Cells_n.mzXML"
                                              ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="NNC dmalQ",
                                        files=[path+"/NNC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_25_1A_Cells_n.mzXML",
                                               path+"/NNC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_27_1B_Cells_n.mzXML"
                                              ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="NNC dmalQDmaa",
                                        files=[path+"/NNC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_29_2A_Cells_n.mzXML",
                                               path+"/NNC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_31_2B_Cells_n.mzXML"
                                              ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="NNC dmnaT",
                                        files=[path+"/NNC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_4_B1_Cells_n.mzXML",
                                               path+"/NNC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_B2_Cells_n.mzXML",
                                               path+"/NNC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_6_B3_Cells_n.mzXML"
                                              ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="NNC wildtype",
                                        files=[path+"/NNC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_1_A1_Cells_n.mzXML",
                                               path+"/NNC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_2_A2_Cells_n.mzXML",
                                               path+"/NNC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_3_A3_Cells_n.mzXML"
                                              ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="NHM dmalQ",
                                        files=[path+"/NHM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_34_1AH_Media_n.mzXML",
                                               path+"/NHM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_36_1BH_Media_n.mzXML"
                                              ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="NHM dmalQDmaa",
                                        files=[path+"/NHM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_38_2AH_Media_n.mzXML",
                                               path+"/NHM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_40_2BH_Media_n.mzXML"
                                              ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="NHM dmnaT",
                                        files=[path+"/NHM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_22_B1H_Media_n.mzXML",
                                               path+"/NHM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_23_B2H_Media_n.mzXML",
                                               path+"/NHM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_24_B3H_Media_n.mzXML"
                                              ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="NHM wildtype",
                                        files=[path+"/NHM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_19_A1H_Media_n.mzXML",
                                               path+"/NHM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_20_A2H_Media_n.mzXML",
                                               path+"/NHM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_21_A3H_Media_n.mzXML"
                                              ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="NNM dmalQ",
                                        files=[path+"/NNM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_33_1A_Media_n.mzXML",
                                               path+"/NNM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_35_1B_Media_n.mzXML"
                                              ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="NNM dmalQDmaa",
                                        files=[path+"/NNM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_37_2A_Media_n.mzXML",
                                               path+"/NNM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_39_2B_Media_n.mzXML"
                                              ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="NNM dmnaT",
                                        files=[path+"/NNM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_16_B1_Media_n.mzXML",
                                               path+"/NNM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_17_B2_Media_n.mzXML",
                                               path+"/NNM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_18_B3_Media_n.mzXML"
                                              ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="NNM wildtype",
                                        files=[path+"/NNM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_13_A1_Media_n.mzXML",
                                               path+"/NNM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_14_A2_Media_n.mzXML",
                                               path+"/NNM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_15_A3_Media_n.mzXML"
                                              ],
                                        color="Slategrey"))

        # experimentalGroups.append(Bunch(name="N_blanks",
        #                                 files=[path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_1_Blank_n.mzXML",
        #                                        path + "/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_2_Blank_n.mzXML",
        #                                        path + "/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_3_Blank_n.mzXML",
        #                                        path + "/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_4_Blank_n.mzXML",
        #                                        path + "/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_Blank_n.mzXML",
        #                                        path + "/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_Blank_n_170222024019.mzXML",
        #                                        path + "/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_Blank_n_170222055313.mzXML",
        #                                       ],
        #                                 color="Black"))
        #
        # experimentalGroups.append(Bunch(name="O_acetylserine",
        #                                 files=[path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_1_O_acetylserine_n.mzXML",
        #                                        path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_1_O_acetylserine_n_repeat.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylaspartate",
        #                                 files=[path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_2_N_acetylaspartate_n.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylglutamate",
        #                                 files=[path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_3_N_acetylglutamate_n.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="acetylglycine",
        #                                 files=[path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_4_N_acetylglycine_n.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylgalactosamine",
        #                                 files=[path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_N_acetylgalactosamine_n.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylglucosamine",
        #                                 files=[path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_6_N_acetylglucosamine_n.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylmannosamine",
        #                                 files=[path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_7_N_acetylmannosamine_n.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylserine",
        #                                 files=[path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_8_N_acetylserine_n.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="acetylurea",
        #                                 files=[path+"/N_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_9_acetylurea_n.mzXML"
        #                                       ],
        #                                 color="Orange"))



        fps=TableUtils.readFile("E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/results.txt", sep="\t")


        ogroups=set()
        for row in fps.getData(cols=["OGroup"], getResultsAsBunchObjects=True):
            ogroups.add(row.OGroup)

        metabolitesToPlot=[]

        for ogroup in ogroups:
            features=[]
            for row in fps.getData(cols=["Num", "MZ", "L_MZ", "RT", "ScanEvent", "Type", "Comment", "Ionisation_Mode"], where="OGroup=%d"%ogroup, getResultsAsBunchObjects=True):
                if row.Ionisation_Mode=="-":
                    features.append(Bunch(num=row.Num,
                                          ogroup=ogroup,
                                          mz=row.MZ,
                                          lmz=row.L_MZ,
                                          rt=row.RT,
                                          rtBorders=.33,
                                          filterLine=row.ScanEvent,
                                          comment=row.Type+"; "+row.Comment))
            if len(features)>0:
                metabolitesToPlot.append(Bunch(name=ogroup, ogroup=ogroup, features=features))

        generatePDF(experimentalGroups, metabolitesToPlot, saveTo="E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/featurePairs_negMode.pdf")

    ########################################################################################################################
    ######   Positive mode
    ########################################################################################################################
    if False:
        experimentalGroups=[]
        path = "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs"
        experimentalGroups.append(Bunch(name="PHC dmalQ",
                                        files=[
                                            path + "/PHC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_26_1AH_Cells_p.mzXML",
                                            path + "/PHC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_28_1BH_Cells_p.mzXML"
                                            ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="PHC dmalQDmaa",
                                        files=[
                                            path + "/PHC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_30_2AH_Cells_p.mzXML",
                                            path + "/PHC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_32_2BH_Cells_p.mzXML"
                                            ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="PHC dmnaT",
                                        files=[
                                            path + "/PHC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_10_B1H_Cells_p.mzXML",
                                            path + "/PHC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_11_B2H_Cells_p.mzXML",
                                            path + "/PHC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_12_B3H_Cells_p.mzXML"
                                            ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="PHC wildtype",
                                        files=[
                                            path + "/PHC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_7_A1H_Cells_p.mzXML",
                                            path + "/PHC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_8_A2H_Cells_p.mzXML",
                                            path + "/PHC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_9_A3H_Cells_p.mzXML"
                                            ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="PNC dmalQ",
                                        files=[path + "/PNC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_25_1A_Cells_p.mzXML",
                                               path + "/PNC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_27_1B_Cells_p.mzXML"
                                               ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="PNC dmalQDmaa",
                                        files=[
                                            path + "/PNC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_29_2A_Cells_p.mzXML",
                                            path + "/PNC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_31_2B_Cells_p.mzXML"
                                            ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="PNC dmnaT",
                                        files=[path + "/PNC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_4_B1_Cells_p.mzXML",
                                               path + "/PNC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_B2_Cells_p.mzXML",
                                               path + "/PNC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_6_B3_Cells_p.mzXML"
                                               ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="PNC wildtype",
                                        files=[
                                            path + "/PNC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_1_A1_Cells_p.mzXML",
                                            path + "/PNC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_2_A2_Cells_p.mzXML",
                                            path + "/PNC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_3_A3_Cells_p.mzXML"
                                            ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="PHM dmalQ",
                                        files=[
                                            path + "/PHM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_34_1AH_Media_p.mzXML",
                                            path + "/PHM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_36_1BH_Media_p.mzXML"
                                            ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="PHM dmalQDmaa",
                                        files=[
                                            path + "/PHM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_38_2AH_Media_p.mzXML",
                                            path + "/PHM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_40_2BH_Media_p.mzXML"
                                            ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="PHM dmnaT",
                                        files=[
                                            path + "/PHM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_22_B1H_Media_p.mzXML",
                                            path + "/PHM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_23_B2H_Media_p.mzXML",
                                            path + "/PHM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_24_B3H_Media_p.mzXML"
                                            ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="PHM wildtype",
                                        files=[
                                            path + "/PHM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_19_A1H_Media_p.mzXML",
                                            path + "/PHM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_20_A2H_Media_p.mzXML",
                                            path + "/PHM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_21_A3H_Media_p.mzXML"
                                            ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="PNM dmalQ",
                                        files=[path + "/PNM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_33_1A_Media_p.mzXML",
                                               path + "/PNM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_35_1B_Media_p.mzXML"
                                               ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="PNM dmalQDmaa",
                                        files=[
                                            path + "/PNM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_37_2A_Media_p.mzXML",
                                            path + "/PNM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_39_2B_Media_p.mzXML"
                                            ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="PNM dmnaT",
                                        files=[path + "/PNM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_16_B1_Media_p.mzXML",
                                               path + "/PNM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_17_B2_Media_p.mzXML",
                                               path + "/PNM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_18_B3_Media_p.mzXML"
                                               ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="PNM wildtype",
                                        files=[
                                            path + "/PNM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_13_A1_Media_p.mzXML",
                                            path + "/PNM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_14_A2_Media_p.mzXML",
                                            path + "/PNM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_15_A3_Media_p.mzXML"
                                            ],
                                        color="Slategrey"))

        # experimentalGroups.append(Bunch(name="P_blanks",
        #                                 files=[path+"/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_1_Blank_p.mzXML",
        #                                        path + "/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_1_Blank_p_170221095225.mzXML",
        #                                        path + "/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_2_Blank_p.mzXML",
        #                                        path + "/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_3_Blank_p.mzXML",
        #                                        path + "/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_4_Blank_p.mzXML",
        #                                        path + "/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_Blank_p.mzXML",
        #                                        path + "/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_Blank_p_170221011813.mzXML",
        #                                        path + "/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_Blank_p_170221042956.mzXML",
        #                                       ],
        #                                 color="Black"))
        #
        # experimentalGroups.append(Bunch(name="O_acetylserine",
        #                                 files=[path+"/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_1_O_acetylserine_p.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylaspartate",
        #                                 files=[path+"/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_2_N_acetylaspartate_p.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylglutamate",
        #                                 files=[path+"/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_3_N_acetylglutamate_p.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="acetylglycine",
        #                                 files=[path+"/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_4_N_acetylglycine_p.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylgalactosamine",
        #                                 files=[path+"/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_N_acetylgalactosamine_p.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylglucosamine",
        #                                 files=[path+"/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_6_N_acetylglucosamine_p.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylmannosamine",
        #                                 files=[path+"/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_7_N_acetylmannosamine_p.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="N_acetylserine",
        #                                 files=[path+"/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_8_N_acetylserine_p.mzXML"
        #                                       ],
        #                                 color="Orange"))
        # experimentalGroups.append(Bunch(name="acetylurea",
        #                                 files=[path+"/P_blanks_stds/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_9_acetylurea_p.mzXML"
        #                                       ],
        #                                 color="Orange"))

        fps = TableUtils.readFile("E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/results.txt", sep="\t")

        ogroups = set()
        for row in fps.getData(cols=["OGroup"], getResultsAsBunchObjects=True):
            ogroups.add(row.OGroup)

        metabolitesToPlot = []

        for ogroup in ogroups:
            features = []
            for row in fps.getData(cols=["Num", "MZ", "L_MZ", "RT", "ScanEvent", "Type", "Comment", "Ionisation_Mode"],
                                   where="OGroup=%d" % ogroup, getResultsAsBunchObjects=True):
                if row.Ionisation_Mode == "+":
                    features.append(Bunch(num=row.Num,
                                          ogroup=ogroup,
                                          mz=row.MZ,
                                          lmz=row.L_MZ,
                                          rt=row.RT,
                                          rtBorders=.33,
                                          filterLine=row.ScanEvent,
                                          comment=row.Type + "; " + row.Comment))
            if len(features) > 0:
                metabolitesToPlot.append(Bunch(name=ogroup, ogroup=ogroup, features=features))

        generatePDF(experimentalGroups, metabolitesToPlot,
                    saveTo="E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/featurePairs_posMode.pdf")

    ########################################################################################################################
    ######   Experiment Benedikt Warth
    ########################################################################################################################
    if False:
        experimentalGroups = []

        experimentalGroups.append(Bunch(name="HEP DON",
                                        files=[
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON3_SUP_034_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON_INCUB_SOL_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON_INCUB_SOL_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON1_LYS_011_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON1_LYS_011_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON1_SUP_026_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON1_SUP_026_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON2_LYS_015_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON2_LYS_015_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON2_SUP_030_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON2_SUP_030_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON3_LYS_019_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON3_LYS_019_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_DON/HEP_DON3_SUP_034_neg.mzXML"
                                        ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="HEP SOL",
                                        files=[
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL1_LYS_012_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL1_LYS_012_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL1_SUP_027_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL1_SUP_027_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL2_LYS_016_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL2_LYS_016_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL2_SUP_031_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL2_SUP_031_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL3_LYS_020_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL3_LYS_020_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL3_SUP_035_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL3_SUP_035_pos.mzXML"
                                        ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="HEP SOL",
                                        files=[
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL1_LYS_012_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL1_LYS_012_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL1_SUP_027_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL1_SUP_027_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL2_LYS_016_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL2_LYS_016_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL2_SUP_031_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL2_SUP_031_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL3_LYS_020_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL3_LYS_020_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL3_SUP_035_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HEP_SOL/HEP_SOL3_SUP_035_pos.mzXML"
                                        ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="HT DON",
                                        files=[
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON3_SUP_032_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON_INBUB_SOL_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON_INBUB_SOL_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON1_LYS_007b_equil_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON1_LYS_007c_equil_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON1_LYS_009_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON1_LYS_009_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON1_SUP_024_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON1_SUP_024_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON2_LYS_013_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON2_LYS_013_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON2_SUP_028_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON2_SUP_028_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON3_LYS_017_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON3_LYS_017_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_DON/HT_DON3_SUP_032_neg.mzXML"
                                        ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="HT SOL",
                                        files=[
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL3_SUP_033_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL1_LYS_010_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL1_LYS_010_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL1_SUP_025_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL1_SUP_025_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL2_LYS_014_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL2_LYS_014_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL2_SUP_029_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL2_SUP_029_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL3_LYS_018_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL3_LYS_018_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/HT_SOL/HT_SOL3_SUP_033_neg.mzXML"
                                        ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="Standards",
                                        files=[
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD01_DONmix_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD1_DONmix_004_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD1_DONmix_004_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD1_DONmix_038_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD1_DONmix_038_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD1_DONmix_ACN_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD1_DONmix_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD10_DONmix_005_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD10_DONmix_005_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD10_DONmix_039_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD10_DONmix_039_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD10_DONmix_ACN_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD10_DONmix_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD100_DONmix_006_022_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD100_DONmix_006_022_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD100_DONmix_006_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD100_DONmix_006_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD100_DONmix_007_023_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD100_DONmix_007_023_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD100_DONmix_040_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD100_DONmix_040_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD100_DONmix_ACN_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD100_DONmix_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD1000_DONmix__ACN_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD1000_DONmix_ACN_2_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD1000_DONmix_ACN_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD1000_DONmix_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/2B4_DON_1_180320_042_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/2B4_DON_1_180320_042_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/2B7_DON_3_180320_043_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/2B7_DON_3_180320_043_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/3ADON_100ppb_045_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/3ADON_100ppb_045_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/15ADON_aprox50ppm_046_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/15ADON_aprox50ppm_046_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/D15GlcA_Met_cage_3_2011_044_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/D15GlcA_Met_cage_3_2011_044_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD01_DONmix_003_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD01_DONmix_003_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD01_DONmix_037_neg.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD01_DONmix_037_pos.mzXML",
                                            "E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/data/Standards/STD01_DONmix_ACN_neg.mzXML"
                                        ],
                                        color="Olivedrab"))



        fps = TableUtils.readFile("E:/180523_433_DON_Warth/DONConjugatesDB.tsv", sep="\t")

        metabolitesToPlot = []
        for row in fps.getData(cols=["Num", "Name", "Rt_min", "MonoisotopicMass"], getResultsAsBunchObjects=True):

            border=.33
            try:
                row.Rt_min=float(row.Rt_min)
            except:
                row.Rt_min=7
                border=7

            mass=float(row.MonoisotopicMass)
            features=[]
            features.append(Bunch(num="%s_%d"%(row.Num, 1),
                                  ogroup=row.Name,
                                  mz=mass+1.007276,
                                  lmz=mass+1.007276+15*1.00335484,
                                  rt=row.Rt_min,
                                  ionisationMode="+",
                                  xn=15,
                                  rtBorders=border,
                                  filterLine="Q Exactive (MS lvl: 1, pol: +)",
                                  comment="[M+H]+"))
            features.append(Bunch(num="%s_%d"%(row.Num, 2),
                                  ogroup=row.Name,
                                  mz=mass+22.989218,
                                  lmz=mass+22.989218+15*1.00335484,
                                  rt=row.Rt_min,
                                  ionisationMode="+",
                                  xn=15,
                                  rtBorders=border,
                                  filterLine="Q Exactive (MS lvl: 1, pol: +)",
                                  comment="[M+Na]+"))
            features.append(Bunch(num="%s_%d"%(row.Num, 3),
                                  ogroup=row.Name,
                                  mz=mass-1.007276,
                                  lmz=mass-1.007276+15*1.00335484,
                                  rt=row.Rt_min,
                                  ionisationMode="-",
                                  xn=15,
                                  rtBorders=border,
                                  filterLine="Q Exactive (MS lvl: 1, pol: -)",
                                  comment="[M-H]-"))

            features.append(Bunch(num="%s_%d"%(row.Num, 1),
                                  ogroup=row.Name,
                                  mz=mass+1.007276,
                                  lmz=mass+1.007276+15*1.00335484,
                                  rt=7,
                                  ionisationMode="+",
                                  xn=15,
                                  rtBorders=7,
                                  filterLine="Q Exactive (MS lvl: 1, pol: +)",
                                  comment="[M+H]+"))
            features.append(Bunch(num="%s_%d"%(row.Num, 2),
                                  ogroup=row.Name,
                                  mz=mass+22.989218,
                                  lmz=mass+22.989218+15*1.00335484,
                                  rt=7,
                                  ionisationMode="+",
                                  xn=15,
                                  rtBorders=7,
                                  filterLine="Q Exactive (MS lvl: 1, pol: +)",
                                  comment="[M+Na]+"))
            features.append(Bunch(num="%s_%d"%(row.Num, 3),
                                  ogroup=row.Name,
                                  mz=mass-1.007276,
                                  lmz=mass-1.007276+15*1.00335484,
                                  rt=7,
                                  ionisationMode="-",
                                  xn=15,
                                  rtBorders=7,
                                  filterLine="Q Exactive (MS lvl: 1, pol: -)",
                                  comment="[M-H]-"))
            metabolitesToPlot.append(Bunch(name=row.Name, ogroup=row.Name, features=features))

        generatePDF(experimentalGroups, metabolitesToPlot,
                    saveTo="E:/180523_433_DON_Warth/EVAL_LowAbundanceSettings/featurePairs.pdf")

    ########################################################################################################################
    ######   Orbitrap problem July 2018
    ########################################################################################################################
    if False:
        experimentalGroups = []

        experimentalGroups.append(Bunch(name="170508_350_MethylatedCompounds",
                                        files=[
                                            "D:/20180708_OrbitrapProblem/170515_350_QC_D.mzXML",
                                            "D:/20180708_OrbitrapProblem/170508_350_QC_A.mzXML",
                                            "D:/20180708_OrbitrapProblem/170508_350_QC_B.mzXML",
                                            "D:/20180708_OrbitrapProblem/170508_350_QC_C.mzXML",
                                            "D:/20180708_OrbitrapProblem/170508_350_QC_D.mzXML",
                                            "D:/20180708_OrbitrapProblem/170509_350_QC_B.mzXML",
                                            "D:/20180708_OrbitrapProblem/170509_350_QC_C.mzXML",
                                            "D:/20180708_OrbitrapProblem/170509_350_QC_D.mzXML",
                                            "D:/20180708_OrbitrapProblem/170510_350_QC_B.mzXML",
                                            "D:/20180708_OrbitrapProblem/170510_350_QC_C.mzXML",
                                            "D:/20180708_OrbitrapProblem/170510_350_QC_D.mzXML",
                                            "D:/20180708_OrbitrapProblem/170511_350_QC_B.mzXML",
                                            "D:/20180708_OrbitrapProblem/170511_350_QC_C.mzXML",
                                            "D:/20180708_OrbitrapProblem/170511_350_QC_D.mzXML",
                                            "D:/20180708_OrbitrapProblem/170512_350_QC_B.mzXML",
                                            "D:/20180708_OrbitrapProblem/170512_350_QC_C.mzXML",
                                            "D:/20180708_OrbitrapProblem/170512_350_QC_D.mzXML",
                                            "D:/20180708_OrbitrapProblem/170513_350_QC_B.mzXML",
                                            "D:/20180708_OrbitrapProblem/170513_350_QC_C.mzXML",
                                            "D:/20180708_OrbitrapProblem/170513_350_QC_D.mzXML",
                                            "D:/20180708_OrbitrapProblem/170514_350_QC_B.mzXML",
                                            "D:/20180708_OrbitrapProblem/170514_350_QC_C.mzXML",
                                            "D:/20180708_OrbitrapProblem/170514_350_QC_D.mzXML",
                                            "D:/20180708_OrbitrapProblem/170515_350_QC_B.mzXML",
                                            "D:/20180708_OrbitrapProblem/170515_350_QC_C.mzXML"
                                        ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="180420_422_MSMS_for_Papers",
                                        files=[
                                            "D:/20180708_OrbitrapProblem/180420_422_QC_5_neg.mzXML",
                                            "D:/20180708_OrbitrapProblem/180420_422_QC_1_neg.mzXML",
                                            "D:/20180708_OrbitrapProblem/180420_422_QC_1_pos.mzXML",
                                            "D:/20180708_OrbitrapProblem/180420_422_QC_2_neg.mzXML",
                                            "D:/20180708_OrbitrapProblem/180420_422_QC_2_pos.mzXML",
                                            "D:/20180708_OrbitrapProblem/180420_422_QC_3_neg.mzXML",
                                            "D:/20180708_OrbitrapProblem/180420_422_QC_3_pos.mzXML",
                                            "D:/20180708_OrbitrapProblem/180420_422_QC_4_neg.mzXML",
                                            "D:/20180708_OrbitrapProblem/180420_422_QC_4_pos.mzXML"
                                        ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="180511_417_Trichoderma_Susanne_Innsbruck_Vorversuch",
                                        files=[
                                            "D:/20180708_OrbitrapProblem/180511_417_QC_4.mzXML",
                                            "D:/20180708_OrbitrapProblem/180511_417_QC_5.mzXML",
                                            "D:/20180708_OrbitrapProblem/180511_417_QC_6.mzXML",
                                            "D:/20180708_OrbitrapProblem/180511_417_QC_1.mzXML",
                                            "D:/20180708_OrbitrapProblem/180511_417_QC_2.mzXML",
                                            "D:/20180708_OrbitrapProblem/180511_417_QC_3.mzXML"
                                        ],
                                        color="Dodgerblue"))

        experimentalGroups.append(Bunch(name="180601_434_Sorbicillinoids",
                                        files=[
                                            "D:/20180708_OrbitrapProblem/180601_434_QC_4.mzXML",
                                            "D:/20180708_OrbitrapProblem/180601_434_QC_1.mzXML",
                                            "D:/20180708_OrbitrapProblem/180601_434_QC_2.mzXML",
                                            "D:/20180708_OrbitrapProblem/180601_434_QC_3.mzXML"
                                        ],
                                        color="Slategrey"))

        experimentalGroups.append(Bunch(name="180615_438_Labelboxpaper",
                                        files=[
                                            "D:/20180708_OrbitrapProblem/180615_438_pos_QC_Std_5.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438_pos_QC_Std_6.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438_pos_QC_Std_1.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438_pos_QC_Std_2.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438_pos_QC_Std_3.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438_pos_QC_Std_4.mzXML"
                                        ],
                                        color="Chocolate"))

        experimentalGroups.append(Bunch(name="180615_438180622_KMT6",
                                        files=[
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_6.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_7.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_8.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_9.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_1.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_2.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_3.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_4.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_5.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_14.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_15.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_16.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_10.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_11.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_12.mzXML",
                                            "D:/20180708_OrbitrapProblem/180615_438180622_QC_13.mzXML"
                                        ],
                                        color="Orange"))

        experimentalGroups.append(Bunch(name="180702_440_sorbicillinoids_posneg",
                                        files=[
                                            "D:/20180708_OrbitrapProblem/180702_440_QC_2.mzXML",
                                            "D:/20180708_OrbitrapProblem/180702_440_QC_3.mzXML",
                                            "D:/20180708_OrbitrapProblem/180702_440_QC_5.mzXML",
                                            "D:/20180708_OrbitrapProblem/180702_440_QC_1.mzXML"
                                        ],
                                        color="Yellowgreen"))


        # t=experimentalGroups
        # experimentalGroups=[]
        # for g in t:
        #     g.files=g.files[1:2]
        #     experimentalGroups.append(g)

        fps = TableUtils.readFile("D:/20180708_OrbitrapProblem/Standards.txt", sep="\t")

        metabolitesToPlot = []
        for row in fps.getData(cols=["Num", "Name", "Rt_min", "MonoisotopicMass"], getResultsAsBunchObjects=True):


            # if len(metabolitesToPlot)>3:
            #     continue

            border = 1.5
            try:
                row.Rt_min = float(row.Rt_min)
            except:
                row.Rt_min = 7
                border = 7

            mass = float(row.MonoisotopicMass)
            features = []
            features.append(Bunch(num="%s_%d" % (row.Num, 1),
                                  ogroup=row.Name,
                                  mz=mass + 1.007276,
                                  lmz=mass + 1.007276 + 15 * 1.00335484,
                                  rt=row.Rt_min,
                                  ionisationMode="+",
                                  xn=15,
                                  rtBorders=border,
                                  filterLine="Q Exactive (MS lvl: 1, pol: +)",
                                  comment="[M+H]+"))
            features.append(Bunch(num="%s_%d" % (row.Num, 2),
                                  ogroup=row.Name,
                                  mz=mass + 22.989218,
                                  lmz=mass + 22.989218 + 15 * 1.00335484,
                                  rt=row.Rt_min,
                                  ionisationMode="+",
                                  xn=15,
                                  rtBorders=border,
                                  filterLine="Q Exactive (MS lvl: 1, pol: +)",
                                  comment="[M+Na]+"))
            features.append(Bunch(num="%s_%d" % (row.Num, 3),
                                  ogroup=row.Name,
                                  mz=mass - 1.007276,
                                  lmz=mass - 1.007276 + 15 * 1.00335484,
                                  rt=row.Rt_min,
                                  ionisationMode="-",
                                  xn=15,
                                  rtBorders=border,
                                  filterLine="Q Exactive (MS lvl: 1, pol: -)",
                                  comment="[M-H]-"))

            metabolitesToPlot.append(Bunch(name=row.Name, ogroup=row.Name, features=features))

        generatePDF(experimentalGroups, metabolitesToPlot,
                    saveTo="D:/20180708_OrbitrapProblem/featurePairs.pdf",
                    intensityCutoff=1E4, plotLabeld=False)

    ########################################################################################################################
    ######   Biological data Maria Phe paper
    ########################################################################################################################
    if False:
        experimentalGroups = []

        experimentalGroups.append(Bunch(name="GML samples",
                                        files=[
                                            "D:/Maria_20171219/160112_posneg_236_12C13C_fullyLab_Remus_untreated_1.mzXML",
                                            "D:/Maria_20171219/160112_posneg_236_12C13C_fullyLab_Remus_untreated_2.mzXML",
                                            "D:/Maria_20171219/160112_posneg_238_12C13C_fullyLab_Remus_untreated_3.mzXML",
                                            "D:/Maria_20171219/160112_posneg_238_12C13C_fullyLab_Remus_untreated_4.mzXML",
                                            "D:/Maria_20171219/160112_posneg_236_12C13C_fullyLab_Remus_untreated_5.mzXML"
                                        ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="PHE samples",
                                        files=[
                                            "D:/Maria_20171219/160112_posneg_236_13C_Phe_1.mzXML",
                                            "D:/Maria_20171219/160112_posneg_236_13C_Phe_2.mzXML",
                                            "D:/Maria_20171219/160112_posneg_236_13C_Phe_3.mzXML"
                                        ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="TRP samples",
                                        files=[
                                            "D:/Maria_20171219/160112_posneg_236_13C_Trp_1.mzXML",
                                            "D:/Maria_20171219/160112_posneg_236_13C_Trp_2.mzXML",
                                            "D:/Maria_20171219/160112_posneg_236_13C_Trp_3.mzXML"
                                        ],
                                        color="Dodgerblue"))

        experimentalGroups.append(Bunch(name="blank samples",
                                        files=[
                                            "D:/Maria_20171219/blank_1.mzXML",
                                            "D:/Maria_20171219/blank_2.mzXML"
                                        ],
                                        color="Slategrey"))


        metabolitesToPlot = []

        features=[]
        group="515_1"
        features.append(Bunch(num="%s" % (group),
                              ogroup=group,
                              mz=177.0550676,
                              lmz=177.055676+9*1.00335484,
                              rt=15.96,
                              ionisationMode="-",
                              xn=9,
                              rtBorders=1,
                              filterLine="Exactive (MS lvl: 1, pol: -)",
                              comment=""))
        metabolitesToPlot.append(Bunch(name=group, ogroup=group, features=features))

        features=[]
        group="518_1"
        features.append(Bunch(num="%s" % (group),
                              ogroup=group,
                              mz=209.0810254,
                              lmz=209.0810254+9*1.00335484,
                              rt=15.96,
                              ionisationMode="+",
                              xn=9,
                              rtBorders=1,
                              filterLine="Exactive (MS lvl: 1, pol: +)",
                              comment=""))
        metabolitesToPlot.append(Bunch(name=group, ogroup=group, features=features))

        features=[]
        group="798_1"
        features.append(Bunch(num="%s" % (group),
                              ogroup=group,
                              mz=303.0863463,
                              lmz=303.0863463+9*1.00335484,
                              rt=21.39,
                              ionisationMode="+",
                              xn=9,
                              rtBorders=1,
                              filterLine="Exactive (MS lvl: 1, pol: +)",
                              comment=""))
        metabolitesToPlot.append(Bunch(name=group, ogroup=group, features=features))

        features=[]
        group="926_1"
        features.append(Bunch(num="%s" % (group),
                              ogroup=group,
                              mz=301.070708,
                              lmz=301.070708+9*1.00335484,
                              rt=23.22,
                              ionisationMode="+",
                              xn=9,
                              rtBorders=1,
                              filterLine="Exactive (MS lvl: 1, pol: +)",
                              comment=""))
        metabolitesToPlot.append(Bunch(name=group, ogroup=group, features=features))

        generatePDF(experimentalGroups, metabolitesToPlot,
                    saveTo="D:/Maria_20171219/featurePairs_missingInEitherExperiment.pdf",
                    intensityCutoff=1E4, plotLabeld=True)

    ########################################################################################################################
    ######   Orbitrap problem July 2018 - 180710_443
    ########################################################################################################################
    if False:
        experimentalGroups = []

        experimentalGroups.append(Bunch(name="ColBIMMA5_QCStds",
                                        files=[
                                            "D:/180710_443_posneg_OrbitrapTests/ColBiMMA5_QC_Std_0k5ppm_1.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/ColBiMMA5_QC_Std_0k5ppm_2.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/ColBiMMA5_QC_Std_0k5ppm_3.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/ColBiMMA5_QC_Std_0k5ppm_4.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/ColBiMMA5_QC_Std_0k5ppm_5.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/ColBiMMA5_QC_Std_0k5ppm_6.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/ColBiMMA5_QC_Std_0k5ppm_7.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/ColBiMMA5_QC_Std_0k5ppm_8.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/ColBiMMA5_QC_Std_0k5ppm_9.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/ColBiMMA5_QC_Std_0k5ppm_10.mzXML"
                                        ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="Col195_QCStds",
                                        files=[
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_QC_Std_0k5ppm_1.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_QC_Std_0k5ppm_2.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_QC_Std_0k5ppm_3.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_QC_Std_0k5ppm_4.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_QC_Std_0k5ppm_5.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_QC_Std_0k5ppm_6.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_QC_Std_0k5ppm_7.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_QC_Std_0k5ppm_8.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_QC_Std_0k5ppm_9.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_QC_Std_0k5ppm_10.mzXML"
                                        ],
                                        color="Olivedrab"))

        experimentalGroups.append(Bunch(name="Col195_12C13CWheat",
                                        files=[
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_12C13CWheat_1.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_12C13CWheat_2.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_12C13CWheat_3.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_12C13CWheat_4.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_12C13CWheat_5.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_12C13CWheat_6.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_12C13CWheat_7.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_12C13CWheat_8.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_12C13CWheat_9.mzXML",
                                            "D:/180710_443_posneg_OrbitrapTests/Col195_12C13CWheat_10.mzXML"
                                        ],
                                        color="Dodgerblue"))

        if False:
            t=experimentalGroups
            experimentalGroups=[]
            for g in t:
                g.files=g.files[1:2]
                experimentalGroups.append(g)

        fps = TableUtils.readFile("D:/180710_443_posneg_OrbitrapTests/Standards.txt", sep="\t")

        metabolitesToPlot = []
        for row in fps.getData(cols=["Num", "Name", "Rt_min", "MonoisotopicMass"], getResultsAsBunchObjects=True):

            # if len(metabolitesToPlot)>3:
            #     continue

            border = 1.5
            try:
                row.Rt_min = float(row.Rt_min)
            except:
                row.Rt_min = 7
                border = 7

            mass = float(row.MonoisotopicMass)
            features = []
            features.append(Bunch(num="%s_%d" % (row.Num, 1),
                                  ogroup=row.Name,
                                  mz=mass + 1.007276,
                                  lmz=mass + 1.007276 + 15 * 1.00335484,
                                  rt=row.Rt_min,
                                  ionisationMode="+",
                                  xn=15,
                                  rtBorders=border,
                                  filterLine="Q Exactive (MS lvl: 1, pol: +)",
                                  comment="[M+H]+"))
            features.append(Bunch(num="%s_%d" % (row.Num, 2),
                                  ogroup=row.Name,
                                  mz=mass + 22.989218,
                                  lmz=mass + 22.989218,
                                  rt=row.Rt_min,
                                  ionisationMode="+",
                                  xn=15,
                                  rtBorders=border,
                                  filterLine="Q Exactive (MS lvl: 1, pol: +)",
                                  comment="[M+Na]+"))
            features.append(Bunch(num="%s_%d" % (row.Num, 3),
                                  ogroup=row.Name,
                                  mz=mass - 1.007276,
                                  lmz=mass - 1.007276,
                                  rt=row.Rt_min,
                                  ionisationMode="-",
                                  xn=15,
                                  rtBorders=border,
                                  filterLine="Q Exactive (MS lvl: 1, pol: -)",
                                  comment="[M-H]-"))

            metabolitesToPlot.append(Bunch(name=row.Name, ogroup=row.Name, features=features))

        generatePDF(experimentalGroups, metabolitesToPlot,
                    saveTo="D:/180710_443_posneg_OrbitrapTests//featurePairs_QCStds.pdf",
                    intensityCutoff=1E4, plotLabeld=False)

    ########################################################################################################################
    ######   Biologcial publication Maria
    ########################################################################################################################
    if True:
        experimentalGroups = []

        experimentalGroups.append(Bunch(name="422",
                                        files=[
                                            "E:/180420_422_MSMS_for_Papers/RawData/EVAL_BiologicalPublication_GMLPHE/MatchMSMSDataToDataMatrix/BiologicalPaper/exampleFiles/422_GML_Remus_DON_MSMS_Rep1_pos.mzXML"
                                        ],
                                        color="Firebrick"))

        experimentalGroups.append(Bunch(name="Reference (243,244) experiment",
                                        files=[
                                            "E:/180420_422_MSMS_for_Papers/RawData/EVAL_BiologicalPublication_GMLPHE/MatchMSMSDataToDataMatrix/BiologicalPaper/exampleFiles/160125_pos_243_13C_Phe_Remus_3.mzXML",
                                            "E:/180420_422_MSMS_for_Papers/RawData/EVAL_BiologicalPublication_GMLPHE/MatchMSMSDataToDataMatrix/BiologicalPaper/exampleFiles/160125_pos_243_13C_Phe_Remus_1.mzXML",
                                            "E:/180420_422_MSMS_for_Papers/RawData/EVAL_BiologicalPublication_GMLPHE/MatchMSMSDataToDataMatrix/BiologicalPaper/exampleFiles/160125_pos_243_13C_Phe_Remus_2.mzXML"
                                        ],
                                        color="Dodgerblue"))

        experimentalGroups.append(Bunch(name="Marc experiment",
                                        files=[
                                            "E:/180420_422_MSMS_for_Papers/RawData/EVAL_BiologicalPublication_GMLPHE/MatchMSMSDataToDataMatrix/BiologicalPaper/exampleFiles/Remus_DON_0h_pooled_12C13C.mzXML"
                                        ],
                                        color="Olivedrab"))



        fps = TableUtils.readFile("E:/180420_422_MSMS_for_Papers/RawData/EVAL_BiologicalPublication_GMLPHE/MatchMSMSDataToDataMatrix/BiologicalPaper/exampleFiles/matchedFeatures.txt", sep="\t")

        metabolitesToPlot = []
        for row in fps.getData(cols=["Num", "Name", "Rt_min", "MonoisotopicMass"], getResultsAsBunchObjects=True):

            # if len(metabolitesToPlot)>3:
            #     continue

            border = 3
            try:
                row.Rt_min = float(row.Rt_min)
            except:
                row.Rt_min = 7
                border = 7

            mass = float(row.MonoisotopicMass)
            features = []
            features.append(Bunch(num="%s_%d" % (row.Num, 1),
                                  ogroup=row.Name,
                                  mz=mass,
                                  lmz=mass,
                                  rt=row.Rt_min,
                                  ionisationMode="+",
                                  xn=15,
                                  rtBorders=border,
                                  filterLine=["Q Exactive (MS lvl: 1, pol: +)", "LTQ Orbitrap XL (MS lvl: 1, pol: +)"],
                                  comment="[M+H]+"))

            metabolitesToPlot.append(Bunch(name=row.Name, ogroup=row.Name, features=features))

        generatePDF(experimentalGroups, metabolitesToPlot,
                    saveTo="E:/180420_422_MSMS_for_Papers/RawData/EVAL_BiologicalPublication_GMLPHE/MatchMSMSDataToDataMatrix/BiologicalPaper/exampleFiles/matchedFeatures.pdf",
                    intensityCutoff=1E4, plotLabeld=False, ppm=10., normEICs=True, showGroupRTShift=False)
