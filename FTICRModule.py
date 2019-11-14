
import FTICR

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

        self.resultsModel=QtGui.QStandardItemModel()
        self.results.setModel(self.resultsModel)

        self.loadedSamples.selectionModel().selectionChanged.connect(self.draw)
        self.results.selectionModel().selectionChanged.connect(self.draw)

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
        self.pl.axes = self.pl.fig.add_subplot(111)
        simpleaxis(self.pl.axes)
        self.pl.twinxs = [self.pl.axes]
        self.pl.mpl_toolbar = NavigationToolbar(self.pl.canvas, self.result_MassSpectrum)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.pl.canvas)
        vbox.addWidget(self.pl.mpl_toolbar)
        self.result_MassSpectrum.setLayout(vbox)




    def processSamples(self):
        matchPPM=self.matchPPM.value()
        clusteringPPM=self.bracketingPPM.value()

        msScans = {}
        for rowi in range(self.samplesModel.rowCount()):
            item = self.samplesModel.item(rowi)
            data=item.data().toPyObject()
            msScans[data.fileName]=data.msScan

        foundSFs = FTICR.processMSScans(msScans, ppm=matchPPM, clusteringPPM=clusteringPPM)
        logging.info("%d signal pairs were annotated with unique sum formulas. %d non-unique sum formulas" % (len(set(foundSFs)), len([f for f in foundSFs if len(f.sfs) > 1])))

        self.resultsModel = QtGui.QStandardItemModel()
        self.results.setModel(self.resultsModel)
        self.results.selectionModel().selectionChanged.connect(self.draw)


        for i, sf in enumerate(sorted(foundSFs, key=lambda sf:sf.meanMZ)):
            item=QtGui.QStandardItem()
            item.setData(sf) ##sfs=sfs, meanMZ=meanmz, cn=cn, files={}
            if len(sf.sfs)==1:
                item.setText(sf.sfs[0] + "    (MZ %.5f)"%sf.meanMZ)
            else:
                item.setText(", ".join(sf.sfs))

            self.resultsModel.appendRow(item)






    def exitApp(self):
        self.close()

    def loadSamples(self):
        if self.samplesModel.rowCount()>0 and \
                        QtGui.QMessageBox.question(self, "MetExtract", "Do you want to keep the already loaded samples? ",
                                      QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.No:
            self.samplesModel=QtGui.QStandardItemModel()
            self.loadedSamples.setModel(self.samplesModel)
            self.loadedSamples.selectionModel().selectionChanged.connect(self.draw)

        filenames = QtGui.QFileDialog.getOpenFileNames(self, caption="Select FT-ICR mass spectra tables",
                                                       directory=self.initDir, filter="tsv (*.tsv);;All files (*.*)")
        for filename in filenames:
            filename=str(filename).replace("\\","/")

            self.initDir = filename
            self.initDir = self.initDir[:self.initDir.rfind("/")]

            msScan=FTICR.importMSScan(filename, colNameMZ="MZ", colNameIntensity="Int")
            logging.info("Read file '%s' with %d signals"%(filename, len(msScan.mz_list)))

            b=Bunch(fileName=filename[filename.rfind("/")+1:], path=filename, color="Firebrick", msScan=msScan)
            item=QtGui.QStandardItem()
            item.setData(b)
            item.setText(b.fileName)

            self.samplesModel.appendRow(item)


    def clearPlot(self):
        ylim = self.pl.twinxs[0].get_ylim()
        xlim = self.pl.twinxs[0].get_xlim()

        for ax in self.pl.twinxs:
            self.pl.fig.delaxes(ax)
        self.pl.axes = self.pl.fig.add_subplot(111)

        simpleaxis(self.pl.axes)

        self.pl.axes.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        self.pl.twinxs = [self.pl.axes]

        return (xlim, ylim)

    def drawCanvas(self, xlim=None, ylim=None):
        if ylim is None:
            ylim = self.pl.twinxs[0].get_ylim()
        if xlim is None:
            xlim = self.pl.twinxs[0].get_xlim()

        for ax in self.pl.twinxs:
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

        self.pl.canvas.draw()

    def drawPlot(self, x=range(10), y=range(1, 11), color="lightgrey",
                 xlab="MZ", ylab="Intensity [counts]"):
        try:
            ax = self.pl.twinxs[0]

            ax.set_title("")
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

        x1 = []
        y1 = []
        x2 = []
        y2 = []

        selIndsSamples=[index.row() for index in self.loadedSamples.selectedIndexes()]
        selIndsResults = [index.row() for index in self.results.selectedIndexes()]
        for rowi in range(self.samplesModel.rowCount()):
            if rowi in selIndsSamples:
                item = self.samplesModel.item(rowi)
                data=item.data().toPyObject()

                for j in range(len(data.msScan.mz_list)):
                    mz=data.msScan.mz_list[j]

                    used=False
                    for rowa in range(self.resultsModel.rowCount()):
                        if rowa in selIndsResults:
                            b=self.resultsModel.item(rowa)
                            b=b.data().toPyObject()

                            for cn in range(0-2, b.cn+3):
                                if abs(mz-b.meanMZ-cn*1.00335484)*1E6/mz < matchppm:
                                    used=True
                                    x2.append(mz);x2.append(mz);x2.append(mz)
                                    y2.append(0);y2.append(data.msScan.intensity_list[j]);y2.append(0)
                    if not used:
                        x1.append(mz);x1.append(mz);x1.append(mz)
                        y1.append(0);y1.append(data.msScan.intensity_list[j]);y1.append(0)


                self.drawPlot(x1,y1, color="whitesmoke" if len(x2)>0 else "lightgrey")
                self.drawPlot(x2, y2, color="Firebrick")

        zoom=True
        if zoom and len(x2)>0:
            xlim=(min(x2)-3.5, max(x2)+3.5)
            ylim=(max(y2)*-0.05, max(y2)*1.05)
        self.drawCanvas(xlim=xlim, ylim=ylim)






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