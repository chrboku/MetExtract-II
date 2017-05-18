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


def generatePageFor(feature, mzXMLs, experimentalGroups, pdf, ppm=5.):
    totalFiles=sum([len(e.files) for e in experimentalGroups])
    filesDone=0

    pdf.drawString(20, 820, "Num: %s OGroup: %s  MZ: %.5f LMZ: %.5f  RT: %.2f ScanEvent: %s" % (feature.num, feature.ogroup, feature.mz, feature.lmz, feature.rt, feature.filterLine))
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
    for i, expGroup in enumerate(experimentalGroups):
        for file in expGroup.files:
            try:
                print ("\rProcessing.. |%%-%ds| (%%d/%%d)"%(totalFiles)) % ("*"*filesDone, filesDone, totalFiles),
                filesDone+=1

                ## Native EICs
                eic, times, scanIds, mzs = mzXMLs[file].getEIC(mz=feature.mz, ppm=ppm, filterLine=feature.filterLine,
                                                              removeSingles=False, intThreshold=0, useMS1=True)
                maxEICValue=max(maxEICValue, max([eic[j] for j in range(len(eic)) if
                               (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)]))

                ## plot overlaied EICs to drawing
                eicsDD.append([(times[j] / 60., eic[j]) for j in range(len(eic)) if
                               (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)])

                ## plot separated EICs to drawing
                eicsSeparatedDD.append([((times[j] + (i+1) * 60) / 60. - feature.rt, eic[j]) for j in range(len(eic)) if
                                        (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (
                                        feature.rt + feature.rtBorders)])

                eicsCols.append(expGroup.color)



                ## Labeled EICs
                eic, times, scanIds, mzsL = mzXMLs[file].getEIC(mz=feature.lmz, ppm=ppm, filterLine=feature.filterLine,
                                                              removeSingles=False, intThreshold=0, useMS1=True)
                maxEICValue=max(maxEICValue, max([eic[j] for j in range(len(eic)) if
                               (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)]))

                ## plot overlaied EICs to drawing
                eicsLDD.append([(times[j] / 60., eic[j]) for j in range(len(eic)) if
                               (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders)])

                ## plot separated EICs to drawing
                eicsLSeparatedDD.append([((times[j] + (i+1) * 60) / 60. - feature.rt, eic[j]) for j in range(len(eic)) if
                                        (feature.rt - feature.rtBorders) <= (times[j] / 60.) <= (
                                        feature.rt + feature.rtBorders)])

                eicsLCols.append(expGroup.color)


                ## find best scan in EIC peak
                bestScanIndex=max([(j, ab) for j, ab in enumerate(eic) if
                                           (feature.rt - feature.rtBorders/3) <= (times[j] / 60.) <= (feature.rt + feature.rtBorders/3)],
                                  key=lambda x: x[1])[0]

                scan=mzXMLs[file].getIthMS1Scan(index=bestScanIndex, filterLine=feature.filterLine)
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
            except Exception:
                print "-- skipping for", file


            # ## get chromatogram area
            # scans, times, scanIDs=\
            #     MzXMLProxy.Instance().getSpecificArea(file, (feature.rt - feature.rtBorders)*60, (feature.rt + feature.rtBorders)*60,
            #                                           feature.mz-3, feature.lmz+3, filterLine=feature.filterLine)
            #
            # for j, scan in enumerate(scans):
            #     rt=times[j]
            #     for k in range(len(scan)):
            #         mz, ab=scan[k]
            #
            #         AreaDD.append([(rt, mz)])
            #         AreaCols.append(expGroup.color)






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
    renderPDF.draw(drawingEICs, pdf, 170, 20)
    pdf.drawString(240, 210, "MS scans (all samples)")


    pdf.drawString(3, 3, "Notes: Raw-data is shown, no processing has been performed, axis labels missing (bug in library), order/color of separated EICs corresponds to group table, names for groups in separated EICs cannot be set")
    # ## Draw chromatogram area
    # drawingSP = Drawing(900, 450)
    # sc = ScatterPlot()
    # sc.x = 10
    # sc.y = 10
    # sc.height = 150
    # sc.width = 180
    # sc.data = AreaDD
    # for i, col in enumerate(AreaCols):
    #     sc.lines[i].symbol = makeMarker("Circle", size=2)
    #     sc.lines[i].symbol.fillColor = HexColor(matplotlib.colors.cnames[col.lower()])
    #     sc.lines[i].symbol.strokeWidth = 0
    #     sc.lines[i].symbol.strokeColor = HexColor(matplotlib.colors.cnames[col.lower()])
    #
    # sc.yLabel = "MZ"
    # sc.xLabel = "RT"
    # sc.lineLabelFormat = noLabel
    # sc.xValueAxis.labelTextFormat = noLabel
    # drawingSP.add(sc)
    # renderPDF.draw(drawingSP, pdf, 820, 400)








    data=[["", "Group number", "Groupname / order"]]
    style=[('FONTSIZE', (0, 0), (-1, -1), 8)]
    for i, expGroup in enumerate(experimentalGroups):
        data.append(["", str(i+1), expGroup.name])
        style.append(('BACKGROUND', (0, len(data)-1), (0, len(data)-1), expGroup.color))

    table = Table(data, style=style)
    w, h = table.wrapOn(pdf, 150, len(data))
    table.drawOn(pdf, 620, 50)


    pdf.showPage()
    print "\r",



def generatePDF(experimentalGroups, metabolitesToPlot, saveTo, mzXMLs=None):

    if mzXMLs is None:
        mzXMLs={}
        print "Importing mzXMLs"
        for i, expGroup in enumerate(experimentalGroups):
            print "   Group: %s"%(expGroup.name)
            for file in expGroup.files:
                print "      File: %s"%(file)

                if file not in mzXMLs.keys():
                    mzXML = Chromatogram()
                    mzXML.parse_file(file)
                    # print mzXML.getFilterLines()
                    mzXMLs[file] = mzXML



    done=0
    ## geneate PDF page for each feature
    print "\nGenerating PDF pages for features"

    pdf=None
    for metabolite in metabolitesToPlot:
        if not done%200:
            if pdf is not None:
                pdf.save()
            pdf = canvas.Canvas(saveTo.replace(".pdf", "_%d.pdf"%(done//50)), pagesize=pagesizes.landscape(pagesizes.A3))
            # pdf = canvas.Canvas(saveTo, pagesize=pagesizes.A4)
        done+=1

        pdf.drawString(20, 820, "Name: %s OGroup: %s  Number of features: %d" % (metabolite.name, metabolite.ogroup, len(metabolite.features)))
        pdf.showPage()

        print "   Metabolite: %s"%(metabolite.ogroup)
        for feature in metabolite.features:
            print "      Num: %s, mz: %.5f, rt: %.2f"%(feature.num, feature.mz, feature.rt)
            generatePageFor(feature, mzXMLs, experimentalGroups, pdf)

    ## save PDF
    pdf.save()




























########################################################################################################################
########################################################################################################################
########################################################################################################################


















if __name__=="__main__" and True:


    ########################################################################################################################
    ######   Negative mode
    ########################################################################################################################
    if True:
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
    if True:
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

