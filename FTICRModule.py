
import FTICR

from mePyGuis.TSVLoaderEdit import *

import logging
import LoggingSetup
LoggingSetup.LoggingSetup.Instance().initLogging()


from formulaTools import formulaTools
from mePyGuis.FTICRWindow import Ui_MainWindow
from PyQt4 import QtGui, QtCore
from utils import Bunch, get_main_dir, is_int, is_float
from multiprocessing import cpu_count
from mePyGuis.ProgressWrapper import ProgressWrapper
from multiprocessing import Pool, freeze_support, cpu_count, Manager



#<editor-fold desc="### MatPlotLib imports and setup">
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.cm import get_cmap

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

def clearPlot(plt, setXtoZero=False):
    for ax in plt.twinxs:
        plt.fig.delaxes(ax)
    plt.axes = plt.fig.add_subplot(111)

    simpleaxis(plt.axes)
    if setXtoZero:
        plt.axes.spines['bottom'].set_position('zero')

    plt.axes.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.twinxs = [plt.axes]
def drawCanvas(plt, ylim=None, xlim=None):

    for ax in plt.twinxs:
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

    plt.canvas.draw()


class FTICRModuleWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, initDir=None):
        logging.info("MetExtract II (module FTICRExtract)")

        self.fT=formulaTools()

        QtGui.QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.setWindowTitle("MetExtract II - FTICRExtract")

        self.initDir="."

        self.samplesModel=QtGui.QStandardItemModel()
        self.loadedSamples.setModel(self.samplesModel)

        self.loadedSamples.selectionModel().selectionChanged.connect(self.draw)
        self.results.itemSelectionChanged.connect(self.draw)

        self.actionExit.triggered.connect(self.exitApp)
        self.actionLoad_samples.triggered.connect(self.loadSamples)

        self.process.clicked.connect(self.processSamples)

        #Setup plot
        #http://eli.thegreenplace.net/2009/01/20/matplotlib-with-pyqt-guis/
        self.pl = QtCore.QObject()
        self.pl.dpi = 50
        self.pl.fig = Figure((5.0, 4.0), dpi=self.pl.dpi, facecolor='white')
        self.pl.fig.subplots_adjust(left=0.05, bottom=0.1, right=0.99, top=0.95)
        self.pl.canvas = FigureCanvas(self.pl.fig)
        self.pl.canvas.setParent(self.result_MassSpectrum)
        self.pl.mpl_toolbar = NavigationToolbar(self.pl.canvas, self.result_MassSpectrum)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.pl.canvas)
        vbox.addWidget(self.pl.mpl_toolbar)
        self.result_MassSpectrum.setLayout(vbox)
        self.clearPlot()


    def processSamples(self):

        msScans = {}
        for rowi in range(self.samplesModel.rowCount()):
            item = self.samplesModel.item(rowi)
            data=item.data().toPyObject()
            msScans[data.fileName]=data.msScan

        if len(msScans)==0:
            QtGui.QMessageBox.warning(self, "MetExtract - FT-ICR modul", "Error: FT-ICR spectra are required for processing. Please laod some samples and then rerun the analysis", QtGui.QMessageBox.Ok)
            return

        saveFile=QtGui.QFileDialog.getSaveFileName(caption="Select results file", directory="./results.tsv",
                                                    filter="TSV file (*.tsv)", )
        saveFile=str(saveFile)

        if saveFile is None or saveFile=="":
            QtGui.QMessageBox.warning(self, "MetExtract - FT-ICR modul", "Error: Results must be saved to a file. Please rerun the analysis and save it to a file", QtGui.QMessageBox.Ok)
            return

        matchPPM=self.matchPPM.value()
        clusteringPPM=self.bracketingPPM.value()
        annotationPPM=self.annotationPPM.value()


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

        pw=ProgressWrapper()
        pw.show()
        foundSFs = FTICR.processMSScans(msScans, ppm=matchPPM, clusteringPPM=clusteringPPM, annotationPPM=annotationPPM,
                                        atoms=atoms, atomsRange=atomsRange,
                                        useSevenGoldenRules=self.useSevenGoldenRules.checkState()==QtCore.Qt.Checked,
                                        setFunctionMax=pw.getCallingFunction()("max"), setFunctionValue=pw.getCallingFunction()("value"), setFunctionText=pw.getCallingFunction()("text"))
        logging.info("%d signal pairs were annotated with unique sum formulas. %d non-unique sum formulas" % (
        len(set(foundSFs)), len([f for f in foundSFs if len(f.sfs) > 1])))

        logging.info("Generating data matrix")
        FTICR.generateDataMatrix(msScans, foundSFs, ppm=matchPPM)

        logging.info("Saving data matrix (%s)"%saveFile)
        FTICR.writeMatrixToFile(saveFile, foundSFs, msScans, ppm=matchPPM)
        pw.hide()

        self.results.clear()

        for i, re in enumerate(sorted(foundSFs, key=lambda sf:sf.meanMZ)):
            sfShow=""
            if len(re.sfs) == 1:
                sfShow=re.sfs[0].ion.replace("M", re.sfs[0].sf)
            elif len(re.sfs) > 1:
                sfShow="> 1"
            item=QtGui.QTreeWidgetItem(["%.5f"%(re.meanMZ), "%d"%(re.cn), sfShow])
            item.data=re

            for sf in re.sfs:
                ## sf=t, ion="[M+Cl]-", theoMZ=mz, ppmError=ppmError
                itemKid=QtGui.QTreeWidgetItem([sf.ion, sf.sf, "%.3f"%sf.ppmError])
                itemKid.data=re
                item.addChild(itemKid)

            self.results.addTopLevelItem(item)






    def exitApp(self):
        self.close()

    def loadSamples(self):
        if self.samplesModel.rowCount()>0 and \
                        QtGui.QMessageBox.question(self, "MetExtract", "Do you want to keep the already loaded samples? ",
                                      QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.No:
            self.samplesModel=QtGui.QStandardItemModel()
            self.loadedSamples.setModel(self.samplesModel)
            self.loadedSamples.selectionModel().selectionChanged.connect(self.draw)

        dialog = TSVLoaderEdit(
            mapping={"mz": "MZ", "Intensity": "Int"},
            order=["mz", "Intensity"])

        dialog.setTitle("Select FT-ICR MS file(s)")
        dialog.setDescription("Select FT-ICR MS file(s)")
        dialog.setWindowTitle("Select FT-ICR MS file(s)")
        dialog.resize(560, 80)

        if dialog.executeDialog() == QtGui.QDialog.Accepted:
            mappingInfo = dialog.getUserSelectedMapping()

            for filename in dialog.getUserSelectedFile():
                filename=str(filename).replace("\\","/")

                self.initDir = filename
                self.initDir = self.initDir[:self.initDir.rfind("/")]

                msScan=FTICR.importMSScan(filename, colNameMZ=mappingInfo["mz"], colNameIntensity=mappingInfo["Intensity"])
                logging.info("Read file '%s' with %d signals"%(filename, len(msScan.mz_list)))

                b=Bunch(fileName=filename[filename.rfind("/")+1:], path=filename, color="Firebrick", msScan=msScan)
                item=QtGui.QStandardItem()
                item.setData(b)
                item.setText(b.fileName)

                self.samplesModel.appendRow(item)


    def clearPlot(self, twinxs=1):
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
            ax = self.pl.fig.add_subplot(twinxs*100+10+i)

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

    def drawPlot(self, x=range(10), y=range(1, 11), color="lightgrey",
                 xlab="MZ", ylab="Intensity [counts]", title="",
                 twinxsindex=1):
        try:
            ax = self.pl.twinxs[twinxsindex]
            ax.set_title(title)
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)
            ax.plot(x, y, color=color)

        except Exception as ex:
            logging.warning("Exception caught", ex)


    def draw(self):
        xlim, ylim=self.clearPlot()
        if xlim[0]==0 and xlim[1]==1 and ylim[0]==0 and ylim[1]==1:
            xlim=None
            ylim=None

        matchppm=self.matchPPM.value()

        x2all = []
        x3all = []

        selIndsSamples=[index.row() for index in self.loadedSamples.selectedIndexes()]
        selectedResults= self.results.selectedItems()

        self.clearPlot(twinxs=max(1, min(len(selIndsSamples), 5)))
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
                x3 = []
                y3 = []

                item = self.samplesModel.item(rowi)
                data=item.data().toPyObject()

                for j in range(len(data.msScan.mz_list)):
                    mz=data.msScan.mz_list[j]

                    used=False
                    for itemSF in selectedResults:
                        b=itemSF.data
                        for cn in range(0-50, b.cn+51):
                            if abs(mz-b.meanMZ-cn*1.00335484)*1E6/mz < matchppm:
                                used=True
                                x2.append(mz);x2.append(mz);x2.append(mz)
                                y2.append(0);y2.append(data.msScan.intensity_list[j]);y2.append(0)

                    if False:
                        for tlii in range(self.results.topLevelItemCount()):
                            tli=self.results.topLevelItem(tlii)
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

                if zoom and len(x2) > 0:
                    ylim = (max(y2 + y3) * -0.05, max(y2 + y3) * 1.05)
                    self.setLims(ylim=ylim, twinxsind=shownSample%5)
                x2all.extend(x2)
                x3all.extend(x3)
                shownSample=shownSample+1

        if zoom and len(x2all) > 0:
            for i in range(len(self.pl.twinxs)):
                xlim = (min(x2all + x3all) - 3.5, max(x2all + x3all) + 3.5)
                self.setLims(xlim=xlim, twinxsind=i)
        self.drawCanvas()






if __name__=="__main__":
    # add freeze support (required for multiprocessing)
    freeze_support()

    import sys

    app = QtGui.QApplication(sys.argv)
    Dialog = FTICRModuleWindow()

    Dialog.show()
    x = app.exec_()

    sys.exit(x)

    window=FTICRModuleWindow()
    window.show()