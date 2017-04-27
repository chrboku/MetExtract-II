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
            eicsSeparatedDD.append([((times[j] + i * 60) / 60., eic[j]) for j in range(len(eic)) if
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
            eicsLSeparatedDD.append([((times[j] + i * 60) / 60., eic[j]) for j in range(len(eic)) if
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








    data=[["", "Groupname / order"]]
    style=[('FONTSIZE', (0, 0), (-1, -1), 8)]
    for i, expGroup in enumerate(experimentalGroups):
        data.append(["", expGroup.name])
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

    for metabolite in metabolitesToPlot:
        if not done%20:
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
    currentWorkspace="E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData"
    experimentalGroups=[]
    # experimentalGroups.append(Bunch(name="Exp2 negMode CellsHeavy dmalQ",
    #                                 files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Cell1AH1.mzXML",
    #                                        "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Cell1BH1.mzXML",
    #                                        ],
    #                                 color="#B22222"))
    # experimentalGroups.append(Bunch(name="Exp2 negMode CellsHeavy dmalQDmaa",
    #                                 files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Cell2BH.mzXML",
    #                                        "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Cell2BH_150917210015.mzXML",
    #                                        ],
    #                                 color="#B22222"))
    # experimentalGroups.append(Bunch(name="Exp2 negMode CellsNative dmalQ",
    #                                 files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Cell1A1.mzXML",
    #                                        "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Cell1B1.mzXML",
    #                                        ],
    #                                 color="#2F4F4F"))
    # experimentalGroups.append(Bunch(name="Exp2 negMode CellsNative dmalQDmaa",
    #                                 files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Cell1A2.mzXML",
    #                                        "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Cell2A.mzXML",
    #                                        "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Cell2B.mzXML"
    #                                        ],
    #                                 color="#2F4F4F"))
    # experimentalGroups.append(Bunch(name="Exp2 negMode SupHeavy dmalQ",
    #                                 files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Sup1AHM.mzXML",
    #                                        "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Sup1BHM.mzXML",
    #                                        ],
    #                                 color="#B22222"))
    # experimentalGroups.append(Bunch(name="Exp2 negMode SupHeavy dmalQDmaa",
    #                                 files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Sup2AHM.mzXML",
    #                                        "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Sup2BHM.mzXML",
    #                                        "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Sup2BM.mzXML"
    #                                        ],
    #                                 color="#B22222"))
    # experimentalGroups.append(Bunch(name="Exp2 negMode SupNative dmalQ",
    #                                 files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Sup1AM.mzXML",
    #                                        "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Sup1BM.mzXML",
    #                                        ],
    #                                 color="#2F4F4F"))
    # experimentalGroups.append(Bunch(name="Exp2 negMode SupNative dmalQDmaa",
    #                                 files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment2/heavy acetate data_2ndProcessing_correctlyConverted/QE1_TG_Sup2AM.mzXML",
    #                                        ],
    #                                 color="#2F4F4F"))




    experimentalGroups.append(Bunch(name="Exp3 negMode CellsHeavy dmalQ",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_26_1AH_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_28_1BH_Cells_n.mzXML"
                                          ],
                                    color="#B22222"))

    experimentalGroups.append(Bunch(name="Exp3 negMode CellsHeavy dmalQDmaa",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_30_2AH_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_32_2BH_Cells_n.mzXML"
                                          ],
                                    color="#B22222"))

    experimentalGroups.append(Bunch(name="Exp3 negMode CellsHeavy dmnaT",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_10_B1H_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_11_B2H_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_12_B3H_Cells_n.mzXML"
                                           ],
                                    color="#B22222"))

    experimentalGroups.append(Bunch(name="Exp3 negMode CellsHeavy wildtype",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_7_A1H_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_8_A2H_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_9_A3H_Cells_n.mzXML"
                                          ],
                                    color="#B22222"))

    experimentalGroups.append(Bunch(name="Exp3 negMode CellsNative dmalQ",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_25_1A_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNC_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_27_1B_Cells_n.mzXML"
                                          ],
                                    color="#2F4F4F"))

    experimentalGroups.append(Bunch(name="Exp3 negMode CellsNative dmalQDmaa",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_29_2A_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNC_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_31_2B_Cells_n.mzXML"
                                          ],
                                    color="#2F4F4F"))

    experimentalGroups.append(Bunch(name="Exp3 negMode CellsNative dmnaT",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_4_B1_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_5_B2_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNC_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_6_B3_Cells_n.mzXML"
                                          ],
                                    color="#2F4F4F"))

    experimentalGroups.append(Bunch(name="Exp3 negMode CellsNative wildtype",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_1_A1_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_2_A2_Cells_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNC_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_3_A3_Cells_n.mzXML"
                                          ],
                                    color="#2F4F4F"))

    experimentalGroups.append(Bunch(name="Exp3 negMode MediumHeavy dmalQ",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_34_1AH_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_36_1BH_Media_n.mzXML"
                                          ],
                                    color="#9ACD32"))

    experimentalGroups.append(Bunch(name="Exp3 negMode MediumHeavy dmalQDmaa",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_38_2AH_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_40_2BH_Media_n.mzXML"
                                          ],
                                    color="#9ACD32"))

    experimentalGroups.append(Bunch(name="Exp3 negMode MediumHeavy dmnaT",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_22_B1H_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_23_B2H_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_24_B3H_Media_n.mzXML"
                                          ],
                                    color="#9ACD32"))

    experimentalGroups.append(Bunch(name="Exp3 negMode MediumHeavy wildtype",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_19_A1H_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_20_A2H_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NHM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_21_A3H_Media_n.mzXML"
                                          ],
                                    color="#9ACD32"))

    experimentalGroups.append(Bunch(name="Exp3 negMode MediumNative dmalQ",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_33_1A_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNM_dmalQ/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_35_1B_Media_n.mzXML"
                                          ],
                                    color="#2F4F4F"))

    experimentalGroups.append(Bunch(name="Exp3 negMode MediumNative dmalQDmaa",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_37_2A_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNM_dmalQDmaa/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_39_2B_Media_n.mzXML"
                                          ],
                                    color="#2F4F4F"))

    experimentalGroups.append(Bunch(name="Exp3 negMode MediumNative dmnaT",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_16_B1_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_17_B2_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNM_dmnaT/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_18_B3_Media_n.mzXML"
                                          ],
                                    color="#2F4F4F"))

    experimentalGroups.append(Bunch(name="Exp3 negMode MediumNative wildtype",
                                    files=["E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_13_A1_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_14_A2_Media_n.mzXML",
                                           "E:/Hanson_MetaboliteRepair_AcetylatedSamples/Experiment3/RawData/mzXMLs/NNM_wildtype/QE2_jdg_62_Hanson_Niehaus_Hanson_Niehaus_15_A3_Media_n.mzXML"
                                          ],
                                    color="#2F4F4F"))


    metabolitesToPlot=[]

    features=[]
    features.append(Bunch(num="a", ogroup="1",
                            mz=419.09528,
                            lmz=419.09528 +2*1.00335484+3*1.00628,
                            rt=3.83,
                            rtBorders=.33,
                            filterLine="Q Exactive (MS lvl: 1, pol: -)",
                            comment="Annotated as [Acetylmaltose+Cl]-"))
    features.append(Bunch(num="a_37Cl", ogroup="1",
                            mz=421.09321,
                            lmz=421.09321 +2*1.00335484+3*1.00628,
                            rt=3.83,
                            rtBorders=.33,
                            filterLine="Q Exactive (MS lvl: 1, pol: -)",
                            comment="Annotated as [Acetylmaltose+37Cl]-"))
    features.append(Bunch(num="b", ogroup="1",
                            mz=429.124744925,
                            lmz=429.124744925 +2*1.00335484+3*1.00628,
                            rt=3.83,
                            rtBorders=.33,
                            filterLine="Q Exactive (MS lvl: 1, pol: -)",
                            comment="Annotated as [Acetylmaltose+FA-H]-"))
    metabolitesToPlot.append(Bunch(name="unknown1", ogroup="1",
                                   features=features))

    features=[]
    features.append(Bunch(num="a", ogroup="2",
                            mz=256.059379,
                            lmz=256.059379 +2*1.00335484+3*1.00628,
                            rt=4.78,
                            rtBorders=.25,
                            filterLine="Q Exactive (MS lvl: 1, pol: -)",
                            comment=""))
    features.append(Bunch(num="b", ogroup="2",
                            mz=258.0556897,
                            lmz=258.0556897 +2*1.00335484+3*1.00628,
                            rt=4.78,
                            rtBorders=.25,
                            filterLine="Q Exactive (MS lvl: 1, pol: -)",
                            comment=""))
    metabolitesToPlot.append(Bunch(name="unknown2", ogroup="2",
                                   features=features))

    features=[]
    features.append(Bunch(num="a", ogroup="3",
                            mz=274.093205981,
                            lmz=274.093205981 +2*1.00335484+3*1.00628,
                            rt=6.73,
                            rtBorders=.33,
                            filterLine="Q Exactive (MS lvl: 1, pol: -)",
                            comment=""))
    features.append(Bunch(num="b", ogroup="3",
                            mz=310.069820134,
                            lmz=310.069820134 +2*1.00335484+3*1.00628,
                            rt=6.73,
                            rtBorders=.33,
                            filterLine="Q Exactive (MS lvl: 1, pol: -)",
                            comment=""))
    features.append(Bunch(num="c", ogroup="3",
                            mz=312.066038469,
                            lmz=312.066038469 +2*1.00335484+3*1.00628,
                            rt=6.73,
                            rtBorders=.33,
                            filterLine="Q Exactive (MS lvl: 1, pol: -)",
                            comment=""))
    metabolitesToPlot.append(Bunch(name="unknown3", ogroup="3",
                                   features=features))

    generatePDF(experimentalGroups, metabolitesToPlot, saveTo=currentWorkspace + "/featurePlot.pdf")


