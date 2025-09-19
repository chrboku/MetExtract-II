
import FTICR

from mePyGuis.TSVLoaderEdit import *

import logging
import LoggingSetup
LoggingSetup.LoggingSetup.Instance().initLogging()

import os


from formulaTools import formulaTools
from mePyGuis.FTICRWindow import Ui_MainWindow
from PySide6 import QtCore, QtGui, QtWidgets
from utils import Bunch, mean
from mePyGuis.ProgressWrapper import ProgressWrapper



#<editor-fold desc="### MatPlotLib imports and setup">
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

matplotlib.rcParams['savefig.dpi'] = 300
font = {'size': 16}
matplotlib.rc('font', **font)
globAlpha = 0.05
predefinedColors = ["FireBrick", "YellowGreen", "SteelBlue", "DarkOrange", "Teal", "CornflowerBlue","DarkOliveGreen", "SlateGrey", "CadetBlue", "DarkCyan", "Black", "DarkSeaGreen", "DimGray","GoldenRod", "LightBlue", "MediumTurquoise", "RoyalBlue"]

# show only bottom and left axes in a matplotlib graph
def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

# remove all axes in a matplotlib graph
def noaxis(ax):
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axis('off')

#from mpltools import style
#style.use('ggplot')
#</editor-fold>



class FTICRModuleWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, initDir=None):
        logging.info("MetExtract II (module FTICRExtract)")

        self.fT=formulaTools()
        self.fticr=FTICR.FTICRProcessing()

        QtWidgets.QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.setWindowTitle("MetExtract II - FTICRExtract")

        self.initDir="."

        self.samplesModel=QtGui.QStandardItemModel()
        self.loadedSamples.setModel(self.samplesModel)

        self.loadedSamples.selectionModel().selectionChanged.connect(self.draw)
        self.results.itemSelectionChanged.connect(self.draw)

        self.actionNew_processing.triggered.connect(self.newProcessing)
        self.actionExit.triggered.connect(self.exitApp)
        self.loadSamplesTSVButton.clicked.connect(self.loadSamplesTSV)
        self.loadSamplesMZXMLButton.clicked.connect(self.loadSamplesMZXML)
        self.process.clicked.connect(self.processSamples)
        self.stackedWidget.setCurrentIndex(0)

        self.dbs=[]
        self.addDB.clicked.connect(self.addSFDB)

        #Setup plot
        #http://eli.thegreenplace.net/2009/01/20/matplotlib-with-pyqt-guis/
        self.pl = QtCore.QObject()
        self.pl.dpi = 50
        self.pl.fig = Figure((5.0, 4.0), dpi=self.pl.dpi, facecolor='white')
        self.pl.fig.subplots_adjust(left=0.05, bottom=0.1, right=0.99, top=0.95)
        self.pl.canvas = FigureCanvas(self.pl.fig)
        self.pl.canvas.setParent(self.result_MassSpectrum)
        self.pl.mpl_toolbar = NavigationToolbar(self.pl.canvas, self.result_MassSpectrum)
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.pl.canvas)
        vbox.addWidget(self.pl.mpl_toolbar)
        self.result_MassSpectrum.setLayout(vbox)
        self.clearPlot()



    def newProcessing(self):
        if self.samplesModel.rowCount() > 0 and \
                        QtWidgets.QMessageBox.question(self, "MetExtract",
                                                   "Do you want to discard the already loaded and processed samples? ",
                                                   QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.Yes:

            self.fticr=FTICR.FTICRProcessing()

            self.samplesModel = QtGui.QStandardItemModel()
            self.loadedSamples.setModel(self.samplesModel)
            self.loadedSamples.selectionModel().selectionChanged.connect(self.draw)

            if QtWidgets.QMessageBox.question(self, "MetExtract",
                                                   "Do you want to discard the loaded databases? ",
                                                   QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.Yes:
                model = QtGui.QStandardItemModel()
                self.DBs.setModel(model)
                self.dbs=[]

            self.clearPlot(twinxs=1)
            self.drawCanvas()

            self.stackedWidget.setCurrentIndex(0)



    def addSFDB(self):
        dialog = TSVLoaderEdit(
            mapping={"Name": "Name", "Sum formula": "SumFormula"},
            order=["Name", "Sum formula"])

        dialog.setTitle("Select sum formula database file(s)")
        dialog.setDescription("Select sum formula database file(s)")
        dialog.setWindowTitle("Select sum formula database file(s)")
        dialog.resize(560, 80)

        if dialog.executeDialog() == QtWidgets.QDialog.Accepted:
            mappingInfo = dialog.getUserSelectedMapping()

            model = QtGui.QStandardItemModel()
            self.DBs.setModel(model)

            for filename in dialog.getUserSelectedFile():
                filename = str(filename).replace("\\", "/")
                f=filename[filename.rfind("/") + 1:]

                imported=0

                with open(filename, "rb") as fin:
                    headers={}
                    for rowi, row in enumerate(fin):
                        cells=row.split("\t")
                        if rowi==0:
                            for celli, cell in enumerate(cells):
                                headers[cell]=celli
                        else:
                            b=Bunch(DBName=f,
                                    sumFormula=cells[headers[mappingInfo["Sum formula"]]],
                                    name=cells[headers[mappingInfo["Name"]]])
                            self.dbs.append(b)
                            imported=imported+1

                item = QtGui.QStandardItem("%s: %d entries"%(f, imported))
                model.appendRow(item)


    def loadSamplesTSV(self):

        dialog = TSVLoaderEdit(
            mapping={"mz": "MZ", "Intensity": "Int"},
            order=["mz", "Intensity"])

        dialog.setTitle("Select FT-ICR MS file(s)")
        dialog.setDescription("Select FT-ICR MS file(s)")
        dialog.setWindowTitle("Select FT-ICR MS file(s)")
        dialog.resize(560, 80)

        texts=[]

        if dialog.executeDialog() == QtWidgets.QDialog.Accepted:
            mappingInfo = dialog.getUserSelectedMapping()

            imported=0

            pw = ProgressWrapper(1, parent=self, showIndProgress=False)
            pw.show()
            pw.getCallingFunction()("max")(len(dialog.getUserSelectedFile()))
            pw.getCallingFunction()("value")(0)

            for filename in dialog.getUserSelectedFile():
                filename = str(filename).replace("\\", "/")
                pw.getCallingFunction()("text")("Importing %s"%(filename))

                self.initDir = filename
                self.initDir = self.initDir[:self.initDir.rfind("/")]

                msScan = self.fticr.importMSScan(filename, colNameMZ=mappingInfo["mz"],
                                                 colNameIntensity=mappingInfo["Intensity"])
                mes="Read file '%s' with %d signals" % (filename, len(msScan.mz_list))
                logging.info(mes)
                texts.append(mes)

                b = Bunch(fileName=filename[filename.rfind("/") + 1:], path=filename, color="Firebrick", msScan=msScan)
                item = QtGui.QStandardItem()
                item.setData(b)
                item.setText(b.fileName)

                self.samplesModel.appendRow(item)

                imported=imported+1
                pw.getCallingFunction()("value")(imported)

            pw.hide()

            self.stackedWidget.setCurrentIndex(1)

    def loadSamplesMZXML(self):

        texts=[]

        filenames = QtWidgets.QFileDialog.getOpenFileNames(self, caption="Select mzXML files", dir=self.initDir,
                                                       filter="mzXML (*.mzxml);;mzML (*.mzml);;All files (*.*)")
        imported=0

        pw = ProgressWrapper(1, parent=self, showIndProgress=False)
        pw.show()
        pw.getCallingFunction()("max")(len(filenames))
        pw.getCallingFunction()("value")(0)

        for filename in filenames:
            filename=str(filename).replace("\\", "/")
            pw.getCallingFunction()("text")("Importing %s"%filename)

            msScan = self.fticr.importMSScan(filename)
            mes = "Read file '%s' with %d signals" % (filename, len(msScan.mz_list))
            logging.info(mes)
            texts.append(mes)

            b = Bunch(fileName=filename[filename.rfind("/") + 1:], path=filename, color="Firebrick", msScan=msScan)
            item = QtGui.QStandardItem()
            item.setData(b)
            item.setText(b.fileName)

            self.samplesModel.appendRow(item)

            imported=imported+1
            pw.getCallingFunction()("value")(imported)

        pw.hide()

        if imported>0:
            self.stackedWidget.setCurrentIndex(1)

    def processSamples(self):

        msScans = {}
        for rowi in range(self.samplesModel.rowCount()):
            item = self.samplesModel.item(rowi)
            data=item.data()
            msScans[data.fileName]=data.msScan

        if len(msScans)==0:
            QtWidgets.QMessageBox.warning(self, "MetExtract - FT-ICR modul", "Error: FT-ICR spectra are required for processing. Please laod some samples and then rerun the analysis", QtWidgets.QMessageBox.Ok)
            return

        saveFile=QtWidgets.QFileDialog.getSaveFileName(caption="Select results file", dir="./results.tsv",
                                                    filter="TSV file (*.tsv)", )
        saveFile=str(saveFile)

        if saveFile is None or saveFile=="":
            QtWidgets.QMessageBox.warning(self, "MetExtract - FT-ICR modul", "Error: Results must be saved to a file. Please rerun the analysis and save it to a file", QtWidgets.QMessageBox.Ok)
            return

        enr12C=self.enr12C.value()
        enr13C=self.enr13C.value()
        maxEnrDeviation=self.maxEnrDeviation.value()

        matchPPM=self.matchPPM.value()
        clusteringPPM=self.bracketingPPM.value()
        annotationPPM=self.annotationPPM.value()

        adducts={}
        if self.adduct_MmHm.isChecked():
            adducts["[M-H]-"]=  -1.007276
        if self.adduct_MpClm.isChecked():
            adducts["[M+Cl]-"]= 34.969402
        if self.adduct_MKm2Hm.isChecked():
            adducts["[M+K-2H]-"]= 36.948606
        if self.adduct_MNam2Hm.isChecked():
            adducts["[M+Na-2H]-"]= 20.974666
        if self.adduct_MpFAmHm.isChecked():
            adducts["[M+FA-H]-"]= 44.998201
        if self.adduct_MpBrm.isChecked():
            adducts["[M+Br]-"]= 78.918885

        elemsMin=self.fT.parseFormula(str(self.SGRElements_min.text()))
        elemsMax=self.fT.parseFormula(str(self.SGRElements_max.text()))

        atoms=[]
        atomsRange=[]
        for elem in ["C", "H", "N", "O", "S", "P"]:
            if elem in elemsMax.keys():
                atoms.append(elem)
                atomsRange.append((elemsMin[elem] if elem in elemsMin.keys() else 0, elemsMax[elem]))

        for elem in elemsMax.keys():
            if elem not in ["C", "H", "N", "O", "S", "P"]:
                atoms.append(elem)
                atomsRange.append((elemsMin[elem] if elem in elemsMin.keys() else 0, elemsMax[elem]))

        searchCAtomsRange=str(self.cAtomsCounts.text())
        a=searchCAtomsRange.replace(" ", "").replace(";",",").split(",")
        searchCAtomsRange=[]
        for j in a:
            if "-" in j:
                searchCAtomsRange.extend(range(int(j[0:j.find("-")]), int(j[j.find("-")+1:len(j)])+1))
            elif j!="":
                searchCAtomsRange.append(int(j))


        pw=ProgressWrapper()
        pw.show()
        foundSFs, averageMZErrors = self.fticr.processMSScans(
                                        msScans,
                                        enr12C=enr12C, enr13C=enr13C, maxEnrDeviation=maxEnrDeviation, checkIsotopologPattern=self.checkIsotopologPattern.isChecked(),
                                        intensityThreshold=self.intensityThreshold.value(), searchCAtomsRange=searchCAtomsRange,
                                        ppm=matchPPM, clusteringPPM=clusteringPPM, annotationPPM=annotationPPM,
                                        atoms=atoms, atomsRange=atomsRange, adducts=adducts,
                                        useSevenGoldenRules=self.useSevenGoldenRules.checkState()==QtCore.Qt.Checked,
                                        recalibrateWithUniqueSumFormulas=self.calibrate.checkState()==QtCore.Qt.Checked,
                                        dbs=self.dbs,
                                        setFunctionMax=pw.getCallingFunction()("max"), setFunctionValue=pw.getCallingFunction()("value"), setFunctionText=pw.getCallingFunction()("text"))
        logging.info("%d signal pairs were annotated with unique sum formulas. %d non-unique sum formulas" % (len(set(foundSFs)), len([f for f in foundSFs if len(f.sfs) > 1])))

        logging.info("Generating data matrix")
        self.fticr.generateDataMatrix(msScans, foundSFs, ppm=matchPPM)

        logging.info("Saving data matrix (%s)"%saveFile)
        self.fticr.writeMatrixToFile(saveFile, foundSFs, msScans, ppm=matchPPM)
        pw.hide()

        self.results.clear()

        item=QtWidgets.QTreeWidgetItem(["Calibration"])
        item.data=Bunch(type="other")
        for file in sorted(msScans.keys()):
            item2=QtWidgets.QTreeWidgetItem([file, "%.3f"%(averageMZErrors[file]) if file in averageMZErrors.keys() else "-"])
            item2.data=Bunch(type="other")
            item.addChild(item2)
        self.results.addTopLevelItem(item)

        item=QtWidgets.QTreeWidgetItem(["Distribution PPM Error"])
        item.data=Bunch(type="PPMErrorGenerated")
        self.results.addTopLevelItem(item)

        item=QtWidgets.QTreeWidgetItem(["Enrichment variation"])
        item.data=Bunch(type="EnrichmentVariation")
        self.results.addTopLevelItem(item)

        for i, re in enumerate(sorted(foundSFs, key=lambda sf:sf.meanMZ)):
            ## re=Bunch(sfs=sfs, dbs=dbs, meanMZ=meanmz, cn=cn, files={})
            re.type="SignalPair"
            sfShow=""

            if len(re.dbs) == 1:
                sfShow=re.dbs[0].name
            elif len(re.dbs)>1:
                sfShow=">1"

            if len(re.sfs) == 1:
                x=re.sfs[0].ion.replace("M", re.sfs[0].sf)
                sfShow=x if len(re.dbs)==0 else "%s / %s" % (sfShow, x)
            elif len(re.sfs) > 1:
                x=">1"
                sfShow="> 1" if len(re.dbs)==0 else "%s / %s" % (sfShow, x)

            if len(re.dbs)==0 and len(re.sfs)==0:
                sfShow="-"

            item=QtWidgets.QTreeWidgetItem(["%.5f"%(re.meanMZ), "%d"%(re.cn), sfShow])
            item.data=re

            for sf in re.sfs:
                ## re.sfs=[Bunch(sf=t, ion="[M+Cl]-", theoMZ=mz, ppmError=ppmError), ...]
                itemKid=QtWidgets.QTreeWidgetItem([sf.ion, sf.sf, "%.3f"%sf.ppmError])
                itemKid.data=re
                item.addChild(itemKid)

            for db in re.dbs:
                ## Bunch(ion=adductName, dbName=db.DBName, name=db.name, ppmError=ppmError)
                itemKid=QtWidgets.QTreeWidgetItem([db.ion, db.name, "%.3f"%db.ppmError])
                itemKid.data=re
                item.addChild(itemKid)

            self.results.addTopLevelItem(item)

        self.stackedWidget.setCurrentIndex(2)






    def exitApp(self):
        self.close()

    def clearPlot(self, twinxs=1, shareX=False, shareY=False):
        if twinxs==0:
            twinxs=1
        try:
            ylim = self.pl.twinxs[0].get_ylim()
            xlim = self.pl.twinxs[0].get_xlim()
            for ax in self.pl.twinxs:
                self.pl.fig.delaxes(ax)
        except Exception:
            ylim = [0, 1]
            xlim = [0, 1]

        self.pl.twinxs=[]
        for i in range(1,twinxs+1):
            ax = self.pl.fig.add_subplot(twinxs*100+10+i)#, sharex=self.pl.twinxs[0] if shareX and i>1 else "None", sharey=self.pl.twinxs[0] if shareY and i>1 else "None")

            simpleaxis(ax)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            self.pl.twinxs.append(ax)

        return (xlim, ylim)

    def setLims(self, xlim=None, ylim=None, twinxsind=1):
        ax = self.pl.twinxs[twinxsind]
        if ylim!=None:
            if len(ylim) == 1:
                ax.set_ylim(ylim[0])
            elif len(ylim) == 2:
                ax.set_ylim(ylim[0], ylim[1])
        if xlim!=None:
            if len(xlim) == 1:
                ax.set_xlim(xlim[0])
            else:
                ax.set_xlim(xlim[0], xlim[1])

    def drawCanvas(self):
        self.pl.fig.tight_layout()
        self.pl.canvas.draw()

    def drawPlot(self, x=range(10), y=range(1, 11), color="lightgrey", marker=".", linestyle="solid",
                 xlab="MZ", ylab="Intensity [counts]", title="",
                 twinxsindex=0):
        try:
            ax = self.pl.twinxs[twinxsindex]
            ax.set_title(title)
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
            ax.plot(x, y, color=color, marker=marker, linestyle=linestyle)

        except Exception as ex:
            logging.warning("Exception caught: %s"%(ex.message))


    def draw(self):
        self.clearPlot()

        matchppm=self.matchPPM.value()

        x2all = []
        x3all = []

        selIndsSamples=[index.row() for index in self.loadedSamples.selectedIndexes()]
        selectedResults= self.results.selectedItems()


        types=set()
        for selectedResult in selectedResults:
            types.add(selectedResult.data.type)

        if len(types)>1:
            QtWidgets.QMessageBox.warning(self, "MetExtract - FT-ICR modul",
                                      "Error: Different results cannot be visualized at the same time",
                                      QtWidgets.QMessageBox.Ok)
        elif len(types)==0:
            return
        elif "SignalPair" in types:
                            self.clearPlot(twinxs=max(1, min(len(selIndsSamples), 5)), shareX=True)
                            shownSample=0
                            zoom = True
                            for rowi in range(self.samplesModel.rowCount()):
                                if rowi in selIndsSamples:
                                    x1 = []
                                    y1 = []
                                    x1a = []
                                    y1a = []
                                    x2 = []
                                    y2 = []
                                    x2out = []
                                    y2out = []
                                    x2heteroIsotope = []
                                    y2heteroIsotope = []
                                    x3 = []
                                    y3 = []

                                    item = self.samplesModel.item(rowi)
                                    data=item.data()

                                    for j in range(len(data.msScan.mz_list)):
                                        mz=data.msScan.mz_list[j]

                                        used=False
                                        for itemSF in selectedResults:
                                            b=itemSF.data

                                            for cn in range(0, b.cn+1):
                                                if abs(mz-b.meanMZ-cn*1.00335484)*1E6/mz < matchppm:
                                                    used=True
                                                    x2.append(mz);x2.append(mz);x2.append(mz)
                                                    y2.append(0);y2.append(data.msScan.intensity_list[j]);y2.append(0)
                                            for cn in range(-50, 0)+range(b.cn+1, b.cn+50):
                                                if abs(mz-b.meanMZ-cn*1.00335484)*1E6/mz < matchppm:
                                                    used=True
                                                    x2out.append(mz);x2out.append(mz);x2out.append(mz)
                                                    y2out.append(0);y2out.append(data.msScan.intensity_list[j]);y2out.append(0)
                                            for delta in [ 2.014102- 1.007825, # H
                                                          15.000109-14.003074, # N
                                                          17.999160-15.994915, # O
                                                          33.967867-31.972071, # S
                                                          36.965903-34.968853, # Cl
                                                          ]:
                                                if abs(mz-b.meanMZ-delta)*1E6/mz < matchppm:
                                                    used=True
                                                    x2heteroIsotope.append(mz);x2heteroIsotope.append(mz);x2heteroIsotope.append(mz)
                                                    y2heteroIsotope.append(0);y2heteroIsotope.append(data.msScan.intensity_list[j]);y2heteroIsotope.append(0)
                                                if abs(mz-b.meanMZ-cn*1.00335484-delta)*1E6/mz < matchppm:
                                                    used=True
                                                    x2heteroIsotope.append(mz);x2heteroIsotope.append(mz);x2heteroIsotope.append(mz)
                                                    y2heteroIsotope.append(0);y2heteroIsotope.append(data.msScan.intensity_list[j]);y2heteroIsotope.append(0)

                                        if False:
                                            for tlii in range(self.results.topLevelItemCount()):
                                                tli=self.results.topLevelItem(tlii)
                                                if tli.data.type=="SignalPair":
                                                    bl=tli.data
                                                    for cn in [-2,-1,0,1,2,bl.cn-2,bl.cn-1,bl.cn,bl.cn+1,bl.cn+2]:
                                                        if abs(mz-bl.meanMZ-cn*1.00335484)*1E6/mz < matchppm:
                                                            used=True
                                                            x1a.append(mz);x1a.append(mz);x1a.append(mz)
                                                            y1a.append(0);y1a.append(data.msScan.intensity_list[j]);y1a.append(0)
                                        if not used:
                                            x1.append(mz);x1.append(mz);x1.append(mz)
                                            y1.append(0);y1.append(data.msScan.intensity_list[j]);y1.append(0)

                                    for itemSF in selectedResults:
                                        b = itemSF.data
                                        x3.append(b.meanMZ);x3.append(b.meanMZ);x3.append(b.meanMZ)
                                        y3.append(0);y3.append(1.1*max(y2 if len(y2)>0 else y1));y3.append(0)
                                        x3.append(b.meanMZ+b.cn*1.00335484);x3.append(b.meanMZ+b.cn*1.00335484);x3.append(b.meanMZ+b.cn*1.00335484)
                                        y3.append(0);y3.append(1.1*max(y2 if len(y2)>0 else y1));y3.append(0)
                                    self.drawPlot(x1,y1, color="whitesmoke" if len(x2)>0 else "lightgrey", title=data.fileName, twinxsindex=shownSample%5)
                                    self.drawPlot(x1a,y1a, color="Slategrey" if len(x2)>0 else "lightgrey", title=data.fileName, twinxsindex=shownSample%5)
                                    self.drawPlot(x3, y3, color="Dodgerblue", title=data.fileName, twinxsindex=shownSample%5)
                                    self.drawPlot(x2, y2, color="Firebrick", title=data.fileName, twinxsindex=shownSample%5)
                                    self.drawPlot(x2out, y2out, color="Orange", title=data.fileName, twinxsindex=shownSample%5)
                                    self.drawPlot(x2heteroIsotope, y2heteroIsotope, color="Olivedrab", title=data.fileName, twinxsindex=shownSample%5)

                                    if zoom and len(x2) > 0:
                                        ylim = (max(y2 + y2out + y3) * -0.05, max(y2 + y2out + y3) * 1.05)
                                        self.setLims(ylim=ylim, twinxsind=shownSample%5)
                                    x2all.extend(x2)
                                    x2all.extend(x2out)
                                    x3all.extend(x3)
                                    shownSample=shownSample+1

                            if zoom and len(x2all) > 0:
                                for i in range(len(self.pl.twinxs)):
                                    xlim = (min(x2all + x3all) - 3.5, max(x2all + x3all) + 3.5)
                                    self.setLims(xlim=xlim, twinxsind=i)
                            self.drawCanvas()


        elif "PPMErrorGenerated" in types:
                            self.clearPlot(twinxs=2)

                            xu=[]
                            yu=[]
                            xnu=[]
                            ynu=[]

                            xei={}
                            yei={}

                            for tlii in range(self.results.topLevelItemCount()):
                                tli = self.results.topLevelItem(tlii)
                                if tli.data.type == "SignalPair":
                                    b = tli.data
                                    for sf in b.sfs:
                                        if len(b.sfs)==1:
                                            xu.append(b.meanMZ)
                                            yu.append(sf.ppmError)
                                        else:
                                            xnu.append(b.meanMZ)
                                            ynu.append(sf.ppmError)

                            self.drawPlot(xu, yu, color="Firebrick", linestyle="None",
                                          title="PPM Error generated sum formulas\n(mean from all samples)", xlab="MZ", ylab="PPM Error",
                                          twinxsindex=0)
                            self.drawPlot(xnu, ynu, color="Dodgerblue", linestyle="None",
                                          title="PPM Error generated sum formulas\n(mean from all samples)", xlab="MZ", ylab="PPM Error",
                                          twinxsindex=0)
                            self.pl.twinxs[0].axhline(y=  matchppm, color="black")
                            self.pl.twinxs[0].axhline(y= -matchppm, color="black")

                            for tlii in range(self.results.topLevelItemCount()):
                                tli = self.results.topLevelItem(tlii)
                                if tli.data.type == "SignalPair":
                                    b = tli.data

                                    for file in b.isotopologRatios["mzError"].keys():
                                        if file not in xei.keys():
                                            xei[file]=[]
                                            yei[file]=[]

                                        xei[file].append(b.meanMZ)
                                        yei[file].append(b.isotopologRatios["mzError"][file])

                            for file in xei.keys():
                                color="Firebrick"
                                for selInd in selIndsSamples:
                                    s=str(self.samplesModel.item(selInd).text())
                                    if s.find(file)!=-1:
                                        color="Dodgerblue"
                                self.drawPlot(xei[file], yei[file], color=color, linestyle="solid",
                                              title="PPM Error generated sum formulas\n(individual samples)", xlab="MZ", ylab="PPM Error",
                                              twinxsindex=1)
                            self.pl.twinxs[1].axhline(y=  matchppm, color="black")
                            self.pl.twinxs[1].axhline(y= -matchppm, color="black")

                            self.drawCanvas()


        elif "EnrichmentVariation" in types:
                            self.clearPlot(twinxs=2, shareX=True, shareY=True)

                            xn=[]
                            yn=[]
                            xl=[]
                            yl=[]

                            for tlii in range(self.results.topLevelItemCount()):
                                tli = self.results.topLevelItem(tlii)
                                if tli.data.type == "SignalPair":
                                    b = tli.data

                                    for file in b.isotopologRatios["native"].keys():
                                        xn.append(b.meanMZ)
                                        yn.append(b.isotopologRatios["native"][file])

                                    for file in b.isotopologRatios["labeled"].keys():
                                        xl.append(b.meanMZ)
                                        yl.append(b.isotopologRatios["labeled"][file])

                            self.drawPlot(xn, yn, color="Firebrick", linestyle="None",
                                          title="Error between theoretical ratio for native carbon isotopologs", xlab="MZ", ylab="Delta observed - theoretical (%)",
                                          twinxsindex=0)
                            self.drawPlot(xl, yl, color="Firebrick", linestyle="None",
                                          title="Error between theoretical ratio for 13C-labeled carbon isotopologs", xlab="MZ", ylab="Delta observed - theoretical (%)",
                                          twinxsindex=1)

                            self.pl.twinxs[0].axhline(y=0, color="black")
                            self.pl.twinxs[0].axhline(y= 0.05, color="slategrey")
                            self.pl.twinxs[0].axhline(y=-0.05, color="slategrey")
                            self.pl.twinxs[0].axhline(y=mean(yn, skipExtremes=0.1), color="Firebrick")
                            #self.pl.twinxs[0].text(min(xn)*1.1, mean(yn, skipExtremes=0.1)+0.1, "Mean enrichment %.2f%%"%(100*mean(yn, skipExtremes=0.1)))


                            self.pl.twinxs[1].axhline(y=0, color="black")
                            self.pl.twinxs[1].axhline(y= 0.05, color="slategrey")
                            self.pl.twinxs[1].axhline(y=-0.05, color="slategrey")
                            self.pl.twinxs[1].axhline(y=mean(yl, skipExtremes=0.1), color="Firebrick")
                            #self.pl.twinxs[1].text(min(xl) * 1.1, mean(yl, skipExtremes=0.1) + 0.1,"Mean enrichment %.2f%%" % (100 * mean(yl, skipExtremes=0.1)))
                            self.drawCanvas()





if __name__=="__main__":

    import sys

    app = QtWidgets.QApplication(sys.argv)
    Dialog = FTICRModuleWindow()

    Dialog.show()
    x = app.exec()

    sys.exit(x)
