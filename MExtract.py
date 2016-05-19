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


import logging

import LoggingSetup

LoggingSetup.LoggingSetup.Instance().initLogging()



# EXPERIMENTAL: alternative peak picking algorithm. Currently under development

import sys, pprint
sys.displayhook=pprint.pprint

from mePyGuis.groupEdit import groupEdit
from mePyGuis.adductsEdit import adductsEdit
from mePyGuis.heteroAtomEdit import heteroAtomEdit
from mePyGuis.DependenciesDialog import DependenciesDialog
from mePyGuis.ProgressWrapper import ProgressWrapper

from utils import USEGRADIENTDESCENDPEAKPICKING


from MetExtractII_Main import MetExtractVersion

# experiment types
TRACER=object()
METABOLOME=object()

#<editor-fold desc="### check if R is installed and accessible">
def checkR():
    try:
        import rpy2.robjects as ro              # import RPy2 module
        r = ro.r                                # make R globally accessible

        v = r("R.Version()$version.string")     # if this R-commando is executed, the RPy2 connection to the
                                                # R subprocess has been established
        return True
    except:
        # The R subprocess could not be started / accessed successfully
        return False

def loadRConfFile(path):
    import os
    if os.path.isfile(path+"/RPATH.conf"):
        with open(path+"/RPATH.conf", "rb") as rconf:
            line=rconf.readline()
            os.environ["R_HOME"]=line
            return True
    else:
        return False


__RHOMEENVVAR=""
import os
from utils import get_main_dir
if "R_HOME" in os.environ.keys():
    __RHOMEENVVAR=os.environ["R_HOME"]


# try to load r configuration file (does not require any environment variables or registry keys)
if not loadRConfFile(path=get_main_dir()) or not checkR():
    os.environ["R_HOME"]=get_main_dir()+"/R"

    if checkR():
        with open("RPATH.conf", "wb") as rconf:
            rconf.write(get_main_dir()+"/R")
            tryLoad=False

            # Show a dialog box to the user that R could not be started
            from os import sys
            from PyQt4 import QtGui, QtCore

            app = QtGui.QApplication(sys.argv)

            QtGui.QMessageBox.information(None, "MetExtract",
                      "R successfully configured\nUsing MetExtract R-Installation\nPlease restart",
                      QtGui.QMessageBox.Ok)
            sys.exit(0)
    else:

        os.environ["R_HOME"]=__RHOMEENVVAR
        os.environ["R_HOME_FROM"]="RPATH environment variable"
        if not checkR():

            logging.error("Error: R could not be loaded correctly (No RPATH.conf file or R_HOME environment variable found)\nPlease make sure it is installed and accessible")

            # Show a dialog box to the user that R could not be started
            from os import sys
            from PyQt4 import QtGui, QtCore

            app = QtGui.QApplication(sys.argv)

            if QtGui.QMessageBox.warning(None, "MetExtract",
                                      "Error: R could not be loaded\nPlease make sure it is installed and accessible\n"
                                      "The default installation path is C:\\Program Files\\R\n"
                                      "Do you want to specify the folder?",
                                      QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)==QtGui.QMessageBox.Yes:
                tryLoad=True
                from utils import get_main_dir
                lastDir=get_main_dir()
                while tryLoad:
                    folder = str(QtGui.QFileDialog.getExistingDirectory(None, "Select R-directory (not bin folder)", directory=lastDir))
                    if folder=="":
                        sys.exit(1)
                    else:
                        lastDir=folder
                        os.environ["R_HOME"]=folder
                        if checkR():
                            with open("RPATH.conf", "wb") as rconf:
                                rconf.write(folder)
                                tryLoad=False

                                QtGui.QMessageBox.information(None, "MetExtract",
                                          "R successfully configured\nPlease restart",
                                          QtGui.QMessageBox.Ok)
                                sys.exit(0)
                        else:
                            if QtGui.QMessageBox.warning(None, "MetExtract",
                                          "Error: R could not be loaded from the specified location\n"
                                          "%s\n\n"
                                          "Please make sure it is installed and accessible\n"
                                          "The default installation path is C:\\Program Files\\R\n"
                                          "Do you want to specify the folder?"%folder,
                                          QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)==QtGui.QMessageBox.Yes:
                                pass
                            else:
                                sys.exit(1)
            else:
                sys.exit(1)
else:
    os.environ["R_HOME_FROM"]="RPATH.conf of MetExtract II"
#</editor-fold>
#<editor-fold desc="### Check if R-dependencies are installed. If not try to fetch them from CRAN and Bioconductor">


# Returns version of r subprocess
def getRVersion():
    try:

        import rpy2.robjects as ro              # import RPy2 module
        r = ro.r                                # make R globally accessible

        v = r("R.Version()$version.string")     # get R-Version
        return v[0]
    except:
        logging.error("Error: R could not be loaded..")

# Checks, if necessary R dependencies are installed
def checkRDependencies(r):
    # R function to check, if package is installed
    r("is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])")

    # Dialog showing if the necessary R packages are installed / could be installed successfully
    RPackageAvailable=0
    RPackageCheck=1
    RPackageError=2
    dependencies = {"R": RPackageCheck, "waveslim": RPackageCheck, "signal": RPackageCheck, "ptw": RPackageCheck, "MASS": RPackageCheck, "MassSpecWavelet": RPackageCheck}
    dialog = DependenciesDialog(dependencies=dependencies, statuss={RPackageAvailable: "olivedrab", RPackageCheck: "orange", RPackageError: "firebrick"},
                                dependencyOrder=["R", "waveslim", "signal", "ptw", "MASS", "MassSpecWavelet"])
    dialog.show()
    dialog.setDependencyStatus("R", RPackageAvailable)

    # Check each dependency, if it is installed and update the dialog accordingly
    missingDependency = False
    for dep in ["waveslim", "signal", "ptw", "MASS"]:
        if str(r("is.installed(\"%s\")" % dep)[0]).lower() == "true":
            dialog.setDependencyStatus(dep, RPackageAvailable)
        else:
            logging.info("installing %s.." % dep)
            r("install.packages(\"%s\", repos='http://cran.us.r-project.org')" % dep)
            if str(r("is.installed(\"%s\")" % dep)[0]).lower() == "true":
                dialog.setDependencyStatus(dep, RPackageAvailable)
            else:
                dialog.setDependencyStatus(dep, RPackageError)
                missingDependency = True

    for dep in ["MassSpecWavelet"]:
        if str(r("is.installed(\"%s\")" % dep)[0]).lower() == "true":
            dialog.setDependencyStatus(dep, RPackageAvailable)
        else:
            logging.info("installing %s.." % dep)
            r("source(\"http://bioconductor.org/biocLite.R\")")
            r("biocLite(\"%s\", suppressUpdates=TRUE)" % dep)
            if str(r("is.installed(\"%s\")" % dep)[0]).lower() == "true":
                dialog.setDependencyStatus(dep, RPackageAvailable)
            else:
                dialog.setDependencyStatus(dep, RPackageError)
                missingDependency = True

    dialog.hide()

    return missingDependency
#</editor-fold>


#<editor-fold desc="### Python Standard Library imports">
import sys
import gc
from pickle import loads, dumps
from collections import defaultdict
import base64
import re
import functools
from operator import itemgetter
from multiprocessing import Pool, freeze_support, cpu_count, Manager
from sqlite3 import *
from copy import copy, deepcopy
from xml.parsers.expat import ExpatError
from optparse import OptionParser
#from hashlib import sha256


#</editor-fold>
#<editor-fold desc="### PyQT 4 Imports">
from PyQt4 import QtGui, QtCore
from PyQt4.QtGui import QColor, QListWidgetItem
#</editor-fold>
#<editor-fold desc="### MatPlotLib imports and setup">
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.cm import get_cmap

matplotlib.rcParams['savefig.dpi'] = 300
font = {'size': 16}
matplotlib.rc('font', **font)
globAlpha = 0.05
predefinedColors = ["FireBrick", "YellowGreen", "SteelBlue", "DarkOrange", "Teal", "DarkSlateGrey", "CornflowerBlue","DarkOliveGreen", "SlateGrey", "CadetBlue", "DarkCyan", "Black", "DarkSeaGreen", "DimGray","GoldenRod", "LightBlue", "MediumTurquoise", "RoyalBlue"]

# required for Plotting table
def convertXaToX(x, a):
    return ((1 - a) * 1) + (a * x)

# taken from http://www.scipy.org/Cookbook/Matplotlib/Transformations
# New enough versions have offset_copy (by Eric Firing):
if 'offset_copy' in dir(matplotlib.transforms):
    from matplotlib.transforms import offset_copy

    def offset(ax, x, y):
        return offset_copy(ax.transData, x=x, y=y, units='dots')
else:  # Without offset_copy we have to do some black transform magic
    from matplotlib.transforms import blend_xy_sep_transform, identity_transform

    def offset(ax, x, y):
        # This trick makes a shallow copy of ax.transData (but fails for polar plots):
        trans = blend_xy_sep_transform(ax.transData, ax.transData)
        # Now we set the offset in pixels
        trans.set_offset((x, y), identity_transform())
        return trans

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

#<editor-fold desc="### MZXML">
from Chromatogram import Chromatogram
#</editor-fold>
#<editor-fold desc="### MassSpecWavelet Processing Class Import">
from chromPeakPicking.MassSpecWavelet import MassSpecWavelet
#</editor-fold>
#<editor-fold desc="### RunIdentification Import">
from runIdentification import RunIdentification
#</editor-fold>
#<editor-fold desc="### Group Results Import">
from bracketResults import bracketResults, calculateMetaboliteGroups, getDBSuffix
from annotateResultMatrix import addGroup as grpAdd
from annotateResultMatrix import addStatsColumnToResults
from annotateResultMatrix import performGroupOmit as grpOmit
#</editor-fold>
#<editor-fold desc="### UI Window Imports">
from mePyGuis.mainWindow import Ui_MainWindow
from mePyGuis.adductsEdit import ConfiguredAdduct, ConfiguredElement
from mePyGuis.heteroAtomEdit import ConfiguredHeteroAtom
from mePyGuis.TracerEdit import tracerEdit
from formulaTools import getIsotopeMass
#</editor-fold>
#<editor-fold desc="### Various Imports">
from utils import natSort, ChromPeakPair, getNormRatio, mean, SampleGroup, Bunch, SQLInsert, SQLSelectAsObject, get_main_dir, smoothDataSeries
from utils import FuncProcess, CallBackMethod
import HCA_general

from TableUtils import TableUtils

import pyperclip
#</editor-fold>

#<editor-fold desc="### debug imports">
import pprint
pp = pprint.PrettyPrinter(indent=1)
#</editor-fold>


# helper method for runIdentification and multiprocessing (multi core support)
def runFile(rI):
    rI.identify()


#<editor-fold desc="Re-Integrate Function definitions">
class procAreaInFile:
    # initialise re-integration for one LC-HRMS file
    def __init__(self, forFile, chromPeakFile=None, offset=0):
        if chromPeakFile is None:
            from utils import get_main_dir
            cpf = get_main_dir()+ "./chromPeakPicking/MassSpecWaveletIdentification.r"
            chromPeakFile=cpf

        self.chromPeakFile = chromPeakFile
        self.done = 0
        self.offset = offset
        self.newData = []
        self.forFile = forFile

    # find chromatographic peak for one detected feature
    def processAreaFor(self, colToProc, mz, rt, ppm, scanEvent, maxRTShift, scales, snrTH, smoothingWindow, smoothingWindowSize):

        if colToProc == "":

            eic, times, scanids = self.t.getEIC(mz, ppm, filterLine=scanEvent)
            eicSmoothed = smoothDataSeries(times, eic, windowLen=smoothingWindowSize,window=smoothingWindow)
            ret = self.CP.getPeaksFor(times, eicSmoothed, scales=scales, snrTh=snrTH)

            best = None

            ri = 0
            for r in ret:
                if abs(times[r.peakIndex] / 60. - rt) <= maxRTShift and (
                                best is None or abs(times[r.peakIndex] / 60. - rt) < abs(times[ret[best].peakIndex] / 60. - rt)):
                    best = ri
                ri += 1

            if best is not None:
                return ret[best].peakArea

        return None

    # re-integrate one detected feature pair
    def processArea(self, oldData, colToProc, colmz, colrt, colLmz, colIonMode, positiveScanEvent, negativeScanEvent,
                    ppm, maxRTShift, scales, snrTH, smoothingWindow, smoothingWindowSize):

        scanEvent = ""
        if oldData[self.colInd[colIonMode]] == "+":
            scanEvent = positiveScanEvent
        elif oldData[self.colInd[colIonMode]] == "-":
            scanEvent = negativeScanEvent
        else:
            logging.error("undefined scan event", oldData[self.colInd[colIonMode]])



        r = self.processAreaFor(oldData[self.colInd[colToProc + "_Area_N"]], oldData[self.colInd[colmz]],
                                oldData[self.colInd[colrt]], ppm, scanEvent, maxRTShift, scales, snrTH, smoothingWindow, smoothingWindowSize)

        nFound = False
        if r is not None:
            oldData[self.colInd[colToProc + "_Area_N"]] = r
            nFound = True

        r = self.processAreaFor(oldData[self.colInd[colToProc + "_Area_L"]], oldData[self.colInd[colLmz]],
                                oldData[self.colInd[colrt]], ppm, scanEvent, maxRTShift, scales, snrTH, smoothingWindow, smoothingWindowSize)

        lFound = False
        if r is not None:
            oldData[self.colInd[colToProc + "_Area_L"]] = r
            lFound = True

        if nFound and lFound:
            oldData[self.colInd[colToProc + "_fold"]] = oldData[self.colInd[colToProc + "_Area_N"]] / oldData[self.colInd[colToProc + "_Area_L"]]


        self.done += 1

        return oldData

    # re-integrate all detected feature pairs in a given LC-HRMS data file
    def processFile(self, oldData, colToProc, colmz, colrt, colLmz, colIonMode, positiveScanEvent, negativeScanEvent,
                    ppm, maxRTShift, scales, snrTH, smoothingWindow, smoothingWindowSize):
        logging.info("   Reintegration started for file %s" % self.forFile)

        z = 0
        for oDat in oldData:

            if self.queue is not None: self.queue.put(Bunch(pid=self.pID, mes="value", val=z))
            z += 1

            nDat = [copy(d) for d in oDat]
            try:
                nDat = self.processArea(nDat, colToProc, colmz, colrt, colLmz, colIonMode, positiveScanEvent,
                                        negativeScanEvent, ppm, maxRTShift, scales, snrTH, smoothingWindow, smoothingWindowSize)
            except Exception as ex:
                print ex
                pass

            self.newData.append(nDat)

        logging.info("   Reintegration finished for file %s" % self.forFile)

    # set re-integration parameters for the current sub-process
    def setParams(self, oldDataFile, headers, colToProc, colmz, colrt, colxcount, colloading, colLmz, colIonMode,
                  positiveScanEvent, negativeScanEvent, colnum, ppm, maxRTShift, scales, reintegrateIntensityCutoff, snrTH,
                  smoothingWindow, smoothingWindowSize):
        self.oldDataFile = oldDataFile
        self.headers = headers

        self.colInd = {}
        for i in range(len(self.headers)):
            self.colInd[self.headers[i]] = i

        self.colToProc = colToProc
        self.colmz = colmz
        self.colrt = colrt
        self.colxcount = colxcount
        self.colloading = colloading
        self.colLmz = colLmz
        self.colIonMode = colIonMode
        self.positiveScanEvent = positiveScanEvent
        self.negativeScanEvent = negativeScanEvent
        self.colnum = colnum
        self.ppm = ppm
        self.maxRTShift = maxRTShift
        self.scales = scales
        self.reintegrateIntensityCutoff = reintegrateIntensityCutoff
        self.snrTH = snrTH

        self.smoothingWindow=smoothingWindow
        self.smoothingWindowSize=smoothingWindowSize

        self.queue = None
        self.pID = None

    def setMultiProcessingQueue(self, qu, pID):
        self.queue = qu
        self.pID = pID

    # re-integrate one LC-HRMS data file
    def processFileParamed(self):
        # read all detected and grouped feature pairs
        self.oldData = TableUtils.readFile(self.oldDataFile, delim="\t").getData()

        if self.queue is not None:
            self.queue.put(Bunch(pid=self.pID, mes="start"))
            self.queue.put(Bunch(pid=self.pID, mes="max", val=len(self.oldData)))
            self.queue.put(Bunch(pid=self.pID, mes="value", val=0))
            self.queue.put(Bunch(pid=self.pID, mes="text", val=str(self.forFile)))

        # import LC-HRMS data file
        logging.info("   Reading file %s" % self.forFile)
        self.t = Chromatogram()
        self.t.parse_file(self.forFile, intensityCutoff=self.reintegrateIntensityCutoff)

        # initialise peak picking algorithm
        if not USEGRADIENTDESCENDPEAKPICKING:
                self.CP = MassSpecWavelet(self.chromPeakFile)
        else:
            from chromPeakPicking.GradientPeaks import GradientPeaks
            self.CP=GradientPeaks()                                                                      ## generic gradient descend peak picking
            self.CP=GradientPeaks(minInt=10000, minIntFlanks=10, minIncreaseRatio=.05, expTime=[25, 250]) ## Swiss Orbitrap HF data

        # re-integrate all detected feature pairs from the grouped results in this LC-HRMS data file
        self.processFile(self.oldData, self.colToProc, self.colmz, self.colrt, self.colLmz, self.colIonMode,
                         self.positiveScanEvent, self.negativeScanEvent, self.ppm, self.maxRTShift, self.scales,
                         self.snrTH, self.smoothingWindow, self.smoothingWindowSize)

        # ask to release used memory
        self.t.freeMe()
        gc.collect()

        if self.queue is not None: self.queue.put(Bunch(pid=self.pID, mes="end"))

        return self.forFile, self.newData

    def getIntegratedData(self):
        return self.newData

    # re-integrated data needs to be matched to the correct columns in the grouped results
    def matchIntegratedDataToTable(self, oldData):
        newData={}
        if oldData[self.colToProc + "_Area_N"] == "":
            for r in self.newData:
                if r[self.colInd[self.colnum]] == oldData[self.colnum]:
                    newData[self.colToProc + "_Area_N"] = r[self.colInd[self.colToProc + "_Area_N"]]

        if oldData[self.colToProc + "_Area_L"] == "":
            for r in self.newData:
                if r[self.colInd[self.colnum]] == oldData[self.colnum]:
                    newData[self.colToProc + "_Area_L"] = r[self.colInd[self.colToProc + "_Area_L"]]

        if oldData[self.colToProc + "_fold"] == "":
            for r in self.newData:
                if r[self.colInd[self.colnum]] == oldData[self.colnum]:
                    newData[self.colToProc + "_fold"] = r[self.colInd[self.colToProc + "_fold"]]
        return newData

# helper method for the multiprocessing package
def integrateFile(s):
    return s.processFileParamed()

def interruptReIntegrateFilesProcessing(pool, selfObj):

    if QtGui.QMessageBox.question(selfObj, "MetExtract", "Are you sure you want to cancel?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)==QtGui.QMessageBox.Yes:
        pool.close()
        pool.terminate()
        pool.join()

        selfObj.terminateJobs=True
        logging.info("Processing stopped by user")

        return True
    else:
        return False # don't close progresswrapper and continue processing files

def integrateResultsFile(file, toF, colToProc, colmz, colrt, colxcount, colloading, colLmz, colIonMode, colnum,
                          ppm=5., maxRTShift=0.25, scales=[3,19], reintegrateIntensityCutoff=0, snrTH=1, smoothingWindow=None, smoothingWindowSize=0,
                          positiveScanEvent="None", negativeScanEvent="None", pw=None, selfObj=None, cpus=1):

    if pw is not None: pw.getCallingFunction()("max")(len(colToProc))

    pwOffset=0
    colToProcCur={}
    writeConfig=True
    for c in colToProc.keys():
        colToProcCur[c]=colToProc[c]

        if len(colToProcCur)>=cpus:
            _integrateResultsFile(file, toF, colToProcCur, colmz, colrt, colxcount, colloading, colLmz, colIonMode, colnum,
                                  ppm, maxRTShift, scales, reintegrateIntensityCutoff, snrTH, smoothingWindow, smoothingWindowSize,
                                  positiveScanEvent, negativeScanEvent, pw, selfObj, cpus, pwOffset=pwOffset,
                                  writeConfig=writeConfig)
            pwOffset+=len(colToProcCur)
            colToProcCur={}
            writeConfig=False

    if len(colToProcCur)>0:
        _integrateResultsFile(file, toF, colToProcCur, colmz, colrt, colxcount, colloading, colLmz, colIonMode, colnum,
                              ppm, maxRTShift, scales, reintegrateIntensityCutoff, snrTH, smoothingWindow, smoothingWindowSize,
                              positiveScanEvent, negativeScanEvent, pw, selfObj, cpus, pwOffset=pwOffset,
                              writeConfig=writeConfig)
        colToProcCur={}
        writeConfig=False

# re-integrate all LC-HRMS data files with the grouped feature pairs results
# uses multiprocessing to distribute each LC-HRMS data file to the available CPU cores
def writeReintegrateConfigToDB(curs, maxRTShift, negativeScanEvent, positiveScanEvent, ppm, reintegrateIntensityCutoff,
                               scales):
    SQLInsert(curs, "config", key="FPREINT_ppm", value=str(ppm))
    SQLInsert(curs, "config", key="FPREINT_maxRTShift", value=str(maxRTShift))
    SQLInsert(curs, "config", key="FPREINT_scales", value=str(scales))
    SQLInsert(curs, "config", key="FPREINT_reintegrateIntensityCutoff", value=str(reintegrateIntensityCutoff))
    SQLInsert(curs, "config", key="FPREINT_positiveScanEvent", value=str(positiveScanEvent))
    SQLInsert(curs, "config", key="FPREINT_negativeScanEvent", value=str(negativeScanEvent))


def _integrateResultsFile(file, toF, colToProc, colmz, colrt, colxcount, colloading, colLmz, colIonMode, colnum,
                         ppm=5.,maxRTShift=0.25, scales=[3, 19], reintegrateIntensityCutoff=0, snrTH=1, smoothingWindow=None, smoothingWindowSize=0,
                         positiveScanEvent="None", negativeScanEvent="None",pw=None, selfObj=None, cpus=1, pwOffset=0,
                         writeConfig=True):

    # read results file
    table = TableUtils.readFile(file, delim="\t")

    # connect to results db
    resDB=Bunch(conn=None, curs=None)
    resDB.conn=connect(file+getDBSuffix())
    resDB.curs=resDB.conn.cursor()

    if writeConfig:
        writeReintegrateConfigToDB(resDB.curs, maxRTShift, negativeScanEvent, positiveScanEvent, ppm, reintegrateIntensityCutoff, scales)

    if pw is not None: pw.getCallingFunction()("text")("Integrating missed peaks")
    if pw is not None: pw.getCallingFunction()("value")(0+pwOffset)


    # check if all expected columns are present in the grouped results
    cols = [col.getName() for col in table.getColumns()]
    for col in colToProc.values():
        if str.isdigit(col[0]):
            col = "_" + col

        assert col + "_Area_N" in cols
        assert col + "_Area_L" in cols

    assert colmz in cols
    assert colrt in cols
    assert colxcount in cols
    assert colloading in cols
    assert colLmz in cols
    assert colIonMode in cols
    assert colnum in cols

    # initialise multiprocessing queue
    p = Pool(processes=min(len(colToProc), cpus))  #, maxtasksperchild=1) # only in python >=2.7; experimental
    manager = Manager()
    queue = manager.Queue()

    if pw is not None: pw.setCloseCallback(closeCallBack=CallBackMethod(_target=interruptReIntegrateFilesProcessing, pool=p, selfObj=selfObj).getRunMethod())

    # generate a re-integration object for each supplied LC-HRMS data file
    toProcFiles = []
    filesDict = {}
    i = 1
    for col in colToProc:

        u = procAreaInFile(col)
        s = colToProc[col]
        if str.isdigit(s[0]):
            s = "_" + s

        u.setParams(file, [x.getName() for x in table.getColumns()], s, colmz, colrt, colxcount, colloading,
                    colLmz, colIonMode, positiveScanEvent, negativeScanEvent, colnum, ppm, maxRTShift, scales, reintegrateIntensityCutoff, snrTH,
                    smoothingWindow, smoothingWindowSize)
        u.setMultiProcessingQueue(queue, i)

        toProcFiles.append(u)
        filesDict[col] = u

        i += 1

    # start the multiprocessing pool
    res = p.imap_unordered(integrateFile, toProcFiles)

    # wait until all subprocesses have finished re-integrating their respective LC-HRMS data file
    start = time.time()
    loop = True
    freeSlots = range(min(len(colToProc), cpus))
    assignedThreads = {}
    while loop and not selfObj.terminateJobs:
        completed = res._index
        if completed == len(colToProc):
            loop = False
        else:
            mess = {}
            while not (queue.empty()):
                mes = queue.get(block=False, timeout=1)
                if mes.pid not in mess:
                    mess[mes.pid] = {}
                mess[mes.pid][mes.mes] = mes

            for v in mess.values():
                if "end" in v.keys():
                    mes = v["end"]
                    freeS = assignedThreads[mes.pid]
                    if freeS == -1:
                        logging.error("Something went wrong..")
                        logging.error("Progress bars do not work correctly, but files will be processed and \"finished..\" will be printed..")
                    else:
                        pw.getCallingFunction(assignedThreads[mes.pid] + 1)("text")("")
                        pw.getCallingFunction(assignedThreads[mes.pid] + 1)("value")(0)
                        assignedThreads[mes.pid] = -1
                        freeSlots.append(freeS)

            for v in mess.values():
                if "start" in v.keys():
                    mes = v["start"]
                    if len(freeSlots) > 0:
                        w = freeSlots.pop()
                        assignedThreads[mes.pid] = w
                    else:
                        logging.error("Something went wrong..")
                        logging.error("Progress bars do not work correctly, but files will be processed and \"finished..\" will be printed..")

            for v in mess.values():
                for mes in v.values():
                    if mes.mes in ["log", "max", "value", "text"]:
                        if assignedThreads.has_key(mes.pid):
                            pw.getCallingFunction(assignedThreads[mes.pid] + 1)(mes.mes)(mes.val)
                        else:
                            pass

            elapsed = (time.time() - start) / 60.
            hours = ""
            if elapsed >= 60.:
                if elapsed < 120.:
                    hours = "1 hour "
                else:
                    hours = "%d hours " % (elapsed // 60)

            if pw is not None: pw.getCallingFunction()("value")(completed+pwOffset)
            if pw is not None: pw.getCallingFunction()("text")("%s%.2f min elapsed %d / %d files done (%d parallel)" % (
                hours, elapsed % 60, completed, len(colToProc), min(cpus, len(colToProc))))

            QtGui.QApplication.processEvents();
            time.sleep(.5)

    p.close()
    p.terminate()
    p.join()

    if selfObj.terminateJobs:
        return

    # after all LC-HRMS data files have been re-integrated, match the individual results
    for re in res:
        u = filesDict[re[0]]
        u.newData = re[1]
        table.applyFunction(u.matchIntegratedDataToTable)

    if writeConfig:
        processingParams=Bunch()
        processingParams.FPReintegrate_ppm=ppm
        processingParams.FPReintegrate_maxRTShift=maxRTShift
        processingParams.FPReintegrate_scales=scales
        processingParams.FPReintegrate_cutoff=reintegrateIntensityCutoff
        processingParams.FPReintegrate_positiveScanEvent=positiveScanEvent
        processingParams.FPReintegrate_negativeScanEvent=negativeScanEvent
        processingParams.FPReintegrate_smoothingWindow=smoothingWindow
        processingParams.FPReintegrate_smoothingWindowSize=smoothingWindowSize
        table.addComment("## Reintegration files processing parameters %s"%(processingParams.dumpAsJSon().replace("\"", "'")))

    # save the re-integrated data matrix to the input file
    TableUtils.saveFile(table, toF)

    # close results db
    resDB.conn.commit()
    resDB.curs.close()
    resDB.conn.close()


#</editor-fold>








def interruptIndividualFilesProcessing(selfObj, pool):
    if QtGui.QMessageBox.question(selfObj, "MetExtract", "Are you sure you want to cancel?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)==QtGui.QMessageBox.Yes:
        pool.close()
        pool.terminate()
        pool.join()

        selfObj.terminateJobs=True
        logging.info("Processing stopped by user")

        return True
    else:
        return False # don't close progresswrapper and continue processing files

def interruptBracketingOfFeaturePairs(selfObj, funcProc):
    if QtGui.QMessageBox.question(selfObj, "MetExtract", "Are you sure you want to cancel?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)==QtGui.QMessageBox.Yes:
        funcProc.terminate()
        funcProc.join()

        selfObj.terminateJobs=True
        logging.info("Processing stopped by user")
    else:
        return False # don't close progresswrapper and continue processing files


class mainWindow(QtGui.QMainWindow, Ui_MainWindow):

    #<editor-fold desc="### check group functions">
    def forceUpdateFL(self):
        self.updateLCMSSampleSettings(force=True)


    # check, if a specified mzxml file can be loaded successfully
    def checkFileImport(self, file):

        #fhash="%s_%s"%(file, sha256(open(file, 'rb').read()).hexdigest())  #self.ckeckedLCMSFiles[fhash]=Bunch(parsed=parsed, fls=tm.getFilterLinesPerPolarity(), pols=tm.getPolarities(), tics=tics)
        fhash="%s_NOHash"%(file)  ## ignore hash, files are not likely to change

        if fhash not in self.checkedLCMSFiles.keys():
            f=file.replace("\\","/")
            f=f[f.rfind("/")+1:f.rfind(".")]

            if not(re.match("^[a-zA-Z0-9_]*$", f)):
                self.checkedLCMSFiles[fhash]=Bunch(parsed="Invalid characters in file name")
            else:
                parsed=False
                try:
                    tm = Chromatogram()
                    tm.parse_file(file, ignoreCharacterData=True)

                    parsed=True
                    pols=tm.getPolarities()

                    fls=tm.getFilterLinesPerPolarity()

                    tics={}
                    if '+' in pols:
                        tic, times, scanIDs=tm.getTIC(filterLine=fls['+'].pop())
                        tics['+']=Bunch(tic=tic, times=times)
                    if '-' in pols:
                        tic, times, scanIDs=tm.getTIC(filterLine=fls['-'].pop())
                        tics['-']=Bunch(tic=tic, times=times)

                    self.checkedLCMSFiles[fhash]=Bunch(parsed=parsed, fls=tm.getFilterLinesPerPolarity(), pols=tm.getPolarities(), tics=tics)

                except ExpatError as ex:
                    self.checkedLCMSFiles[fhash]=Bunch(parsed="Parsing error")
                except Exception as ex:
                    self.checkedLCMSFiles[fhash]=Bunch(parsed="General error")

        return self.checkedLCMSFiles[fhash].parsed



    def updateLCMSSampleSettings(self, force=False):

        if not (self.ui.dontUpdateFilterLine.isChecked()) or force:

            # fetch LC-HRMS data (polarity, filter-lines)
            self.ui.processedFilesComboBox.clear()
            self.ui.processedFilesComboBox.addItem("--", Bunch(file=None))
            self.ui.res_ExtractedData.clear()
            self.ui.chromPeakName.setText("")

            filterLines = set()
            polarities = set()

            # find out, how many files are to be inspected
            all = 0
            done = []
            indGroups = {}
            filesDict = {}

            for group in [t.data(QListWidgetItem.UserType).toPyObject() for t in natSort(self.ui.groupsList.findItems('*', QtCore.Qt.MatchWildcard), key=lambda x: str(x.data(QListWidgetItem.UserType).toPyObject().name))]:
                indGroups[group.name] = natSort(group.files)
                for file in natSort(group.files):
                    if file not in done:
                        all += 1
                        done.append(file)

                    if file not in filesDict.keys():
                        filesDict[file]=set()
                    filesDict[file].add(group.name)


            multipleFiles = {}
            for file in filesDict.keys():
                if len(filesDict[file])>1:
                    multipleFiles[file]=filesDict[file]

            if len(multipleFiles)>0:
                fc=[]
                for f in multipleFiles.keys():
                    fc.append(f)
                    for group in multipleFiles[f]:
                        fc.append("  - group %s"%group)

                QtGui.QMessageBox.warning(self, "MetExtract II", "Some files were used more than once. Please check if all groups have been imported / created correctly\n\nThe files are:\n%s"%("\n".join(fc)),
                                          QtGui.QMessageBox.Ok)

            pw = ProgressWrapper(1, parent=self, showIndProgress=True, indGroups=indGroups)
            pw.show()
            pw.getCallingFunction()("max")(all)

            i = 0
            done = []
            commonFilterLines = {}
            color = False
            failed = defaultdict(list)
            self.clearPlot(self.ui.pl_tic)

            # try to import each LC-HRMS file and get its polarities and filter lines (MS only)
            # draw the TICs of the individual filter lines
            for group in [t.data(QListWidgetItem.UserType).toPyObject() for t in natSort(self.ui.groupsList.findItems('*', QtCore.Qt.MatchWildcard), key=lambda x: str(x.data(QListWidgetItem.UserType).toPyObject().name))]:
                for file in natSort(group.files):
                    if file not in done:

                        pw.getCallingFunction()("text")("Importing group " + group.name + "\n" + file)
                        pw.getCallingFunction()("statuscolor")(file, "Orange")
                        pw.getCallingFunction()("statustext")(file, text="File: %s\nStatus: %s" % (file, "importing"))

                        f=file.replace("\\","/")
                        f=f[f.rfind("/")+1:f.rfind(".")]

                        if not(re.match("^[a-zA-Z0-9_]*$", f)):
                            failed["Invalid characters"].append(file)
                            pw.getCallingFunction()("statuscolor")(file, "firebrick")
                            continue


                        if not(os.path.isfile(file)):
                            failed["File not found"].append(file)
                            pw.getCallingFunction()("statuscolor")(file, "firebrick")
                            continue

                        try:
                            if self.checkFileImport(file)==True:

                                #fhash="%s_%s"%(file, sha256(open(file, 'rb').read()).hexdigest())
                                fhash="%s_NOHash"%(file)  ## ignore hash, files are not likely to change
                                b=self.checkedLCMSFiles[fhash]

                                pols=b.pols
                                fls=b.fls

                                if '+' in pols:
                                    self.drawPlot(self.ui.pl_tic, 0, [t/60. for t in b.tics['+'].times], b.tics['+'].tic, useCol=group.color)
                                if '-' in pols:
                                    mult=1
                                    if '+' in pols and '-' in pols:
                                        mult=-1
                                    self.drawPlot(self.ui.pl_tic, 0, [t/60. for t in b.tics['-'].times], [e*mult for e in b.tics['-'].tic], useCol=group.color)

                                for pol in pols:
                                    if len(fls[pol])>0:
                                        if pol not in commonFilterLines:
                                            commonFilterLines[pol]=fls[pol]
                                        else:
                                            commonFilterLines[pol]=commonFilterLines[pol].intersection(fls[pol])

                                for pol in pols:
                                    polarities.add(pol)

                                # if file has been processed successfully (FileName.identified.sqlite DB exists) add it to the successfully processed list
                                if os.path.exists(file + ".identified.sqlite") and os.path.isfile(file + ".identified.sqlite"):
                                    fname=fname=file[(file.rfind("/") + 1):]
                                    self.ui.processedFilesComboBox.addItem(fname, userData=Bunch(file=file))

                                done.append(file)

                                pw.getCallingFunction()("statuscolor")(file, "olivedrab")
                                pw.getCallingFunction()("statustext")(file,
                                                                      text="File: %s\nStatus: %s" % (file, "imported"))
                            else:

                                pw.getCallingFunction()("statuscolor")(file, "firebrick")
                                pw.getCallingFunction()("statustext")(file, text="File: %s\nStatus: %s" % (file, "failed"))
                                failed["General error"].append(file)

                        except ExpatError as ex:
                            pw.getCallingFunction()("statuscolor")(file, "firebrick")
                            pw.getCallingFunction()("statustext")(file, text="File: %s\nStatus: %s" % (file, "failed"))
                            failed["Parsing error"].append(file)
                        except Exception as ex:
                            pw.getCallingFunction()("statuscolor")(file, "firebrick")
                            pw.getCallingFunction()("statustext")(file, text="File: %s\nStatus: %s" % (file, "failed"))
                            failed["General error"].append(file)

                        i += 1
                        pw.getCallingFunction()("value")(i)
                color = not color

            # update the TIC graph
            pw.hide()
            self.drawCanvas(self.ui.pl_tic)

            if len(failed) > 0:
                t=[]
                for q in failed.keys():
                    t.append(q)
                    t.append("\n")
                    for w in failed[q]:
                        t.append("  - ")
                        t.append(w)
                        t.append("\n")
                t="".join(t)
                QtGui.QMessageBox.warning(self, "MetExtract",
                                          "%d failed to import correctly. Please resolve the following problems before continuing.\n\nFailed files:\n%s" % (len(failed), t),
                                          QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
                return False

            if len(commonFilterLines) == 0:
                QtGui.QMessageBox.information(self, "MetExtract",
                                              "No common scan events among files. Please use files with at least one identical measurement method",
                                              QtGui.QMessageBox.Ok)
                return False

            # Create list with positive mode scan events
            qslpos = QtCore.QStringList()
            if "+" in commonFilterLines.keys():
                for fl in commonFilterLines["+"]:
                    filterLines.add(fl)
                    qslpos.append(fl)

            # Create list with negative mode scan events
            qslneg = QtCore.QStringList()
            if "-" in commonFilterLines.keys():
                for fl in commonFilterLines["-"]:
                    filterLines.add(fl)
                    qslneg.append(fl)

            # Add empty scan event to both ionisation modes if positive and negative ionisation scans are available
            if (len(qslneg) == 0 and len(qslpos) == 0) or (len(qslneg) >= 1 and len(qslpos) >= 1):
                qslneg.insert(0, "None")
                qslpos.insert(0, "None")

            # Add empty scan event if only one ionisation mode was used to the opposite ionisation mode
            if len(qslneg) == 0 and len(qslpos) >= 1:
                qslneg.insert(0, "None")
            if len(qslpos) == 0 and len(qslneg) >= 1:
                qslpos.insert(0, "None")

            # Set scan event in user GUI
            self.ui.positiveScanEvent.clear()
            self.ui.positiveScanEvent.addItems(qslpos)
            self.ui.negativeScanEvent.clear()
            self.ui.negativeScanEvent.addItems(qslneg)

            # Set selected scan event. If both ionisation modes have scan events, select the first of both
            # If only one ionisation mode was used, set the first there and set the other ionisation mode to empty
            self.ui.positiveScanEvent.setCurrentIndex(0)
            self.ui.negativeScanEvent.setCurrentIndex(0)
            if len(qslpos) > 1 and len(qslneg) > 1:
                self.ui.positiveScanEvent.setCurrentIndex(1)
                self.ui.negativeScanEvent.setCurrentIndex(1)

            # Allow user to process files and view the results
            self.ui.pickingTab.setEnabled(len(filterLines) > 0)
            self.ui.resultsTab.setEnabled(self.ui.processedFilesComboBox.count() > 0)

            self.checkIndFilesSettings()

            self.loadGroupsResultsFile(str(self.ui.groupsSave.text()))

            return True

    def smoothingWindowChanged(self):
        e=self.ui.eicSmoothingWindow.currentText()!="None"
        self.ui.eicSmoothingWindowSize.setEnabled(e)
        self.ui.smoothingWindowSizeLabel.setEnabled(e)

    def isoSearchChanged(self):
        self.ui.label_73.setVisible(self.ui.isoAbundance.checkState()==QtCore.Qt.Checked)
        self.ui.intensityThresholdIsotopologs.setVisible(self.ui.isoAbundance.checkState()==QtCore.Qt.Checked)
        #if self.ui.isoAbundance.checkState()==QtCore.Qt.Checked:
        #    self.ui.intensityThresholdIsotopologs.setValue(self.ui.intensityThreshold.value())


    # check if all imported LC-HRMS data files were processed with the same MetExtract settings
    def checkIndFilesSettings(self):
        done = []
        inValidfiles = ""
        unprocessed = 0
        i = 0
        settings = {}

        # check each imported LC-HRMS file individually
        for group in [t.data(QListWidgetItem.UserType).toPyObject() for t in natSort(self.ui.groupsList.findItems('*', QtCore.Qt.MatchWildcard), key=lambda x: str(x.data(QListWidgetItem.UserType).toPyObject().name))]:
            for file in natSort(group.files):
                if file not in done:
                    done.append(file)
                    if os.path.exists(file + ".identified.sqlite") and os.path.isfile(file + ".identified.sqlite"):
                        conn = None
                        curs = None
                        try:
                            conn = connect(file + ".identified.sqlite")
                            curs = conn.cursor()
                            hasConfigTable = False

                            # fetch processing configuration (table config)
                            for name, type in curs.execute("SELECT name, type FROM sqlite_master WHERE type='table' AND name='config'"):
                                hasConfigTable = True
                            if hasConfigTable:

                                # create an entry for each processed setting
                                for key, value in curs.execute("SELECT key, value FROM config"):
                                    if not (settings.has_key(key)):
                                        settings[key] = {}
                                    if not (settings[key].has_key(value)):
                                        settings[key][value] = 0
                                    settings[key][value] += 1
                            curs.close()
                            curs = None
                            conn.close()
                            conn = None
                        except:
                            try:
                                if curs is not None:
                                    curs.close()
                                if conn is not None:
                                    conn.close()
                            except:
                                logging.warning("Warning: Could not close intermediate (sql) file")
                    else:
                        if len(inValidfiles) > 0:
                            inValidfiles += "\n"
                        inValidfiles = inValidfiles + file
                        unprocessed += 1
                    i += 1
        # if only one setting key-value pair was used for each LC-HRMS file, the processing settings of all LC-HRMS
        # data files are identical
        corruptParameters = []
        for sett in settings.keys():
            if len(settings[sett]) != 1:
                if sett not in ["heteroAtoms", "adducts", "elements", "tracerConfiguration", "heteroElements", "adducts", "elementsForNL"]:
                    corruptParameters.append("  * "+sett+ ": "+", ".join(settings[sett]))

        if not self.silent:
            if len(inValidfiles) > 0 and unprocessed != i:
                QtGui.QMessageBox.warning(self, "MetExtract","Not all files were processed\nUnprocessed files:\n%s" % inValidfiles,QtGui.QMessageBox.Ok)
        if not self.silent:
            if len(corruptParameters)>0:
                QtGui.QMessageBox.warning(self, "MetExtract","Not all files were processed using the same parameters.\n"
                                                             "Individual files should be reprocessed\n\n%s"%("\n".join(sorted(corruptParameters))),
                                          QtGui.QMessageBox.Ok)
            #self.updateIndividualFileProcessing = False
            #self.ui.processIndividualFiles.setChecked(True)
            #self.updateIndividualFileProcessing = True
    #</editor-fold>

    #<editor-fold desc="### add/modify/remove group functions">
    def addGroup(self, name, files, minGrpFound, omitFeatures, useForMetaboliteGrouping, color, atPos=None):

        failed=defaultdict(list)

        pw = ProgressWrapper(1, parent=self, showIndProgress=True, indGroups={"files":files})
        pw.show()
        pw.getCallingFunction()("max")(len(files))
        pw.getCallingFunction()("text")("Checking group %s"%(name))
        i=0

        for f in files:
            pw.getCallingFunction()("text")("Checking group %s\nFile: %s"%(name, f))
            x = self.checkFileImport(f)
            if x==True:
                pw.getCallingFunction()("statuscolor")(f, "olivedrab")
                pw.getCallingFunction()("statustext")(f,
                                                      text="File: %s\nStatus: %s" % (f, "imported"))
            else:
                failed[x].append(f)
                pw.getCallingFunction()("statuscolor")(f, "Firebrick")
                pw.getCallingFunction()("statustext")(f,
                                                      text="File: %s\nStatus: %s" % (f, x))
            i+=1
            pw.getCallingFunction()("value")(i)
        pw.hide()

        if len(failed)==0:
            if atPos is None:
                atPos=self.ui.groupsList.count()
            qlwi=QListWidgetItem("%s (%d%s%s)"%(str(name), len(files), ", Annotate" if useForMetaboliteGrouping else "", ", Omit/%d"%minGrpFound if omitFeatures else ""))
            qlwi.setBackgroundColor(QColor(color))
            qlwi.setData(QListWidgetItem.UserType, SampleGroup(name, files, minGrpFound, omitFeatures, useForMetaboliteGrouping, color))

            self.ui.groupsList.insertItem(atPos, qlwi)
        else:
            t=[]
            for q in failed.keys():
                t.append(q)
                t.append("\n")
                for w in failed[q]:
                    t.append("  - ")
                    t.append(w)
                    t.append("\n")
            t="".join(t)
            QtGui.QMessageBox.warning(self, "MetExtract",
                                      "%d files failed to import correctly. Group '%s' will not be created.\n"
                                      "Please resolve the following issues (see documentation for further help)\n\n%s" % (len(failed), name, t), QtGui.QMessageBox.Ok)


    def showAddGroupDialogClicked(self):
        self.showAddGroupDialog()

    # show an import group dialog to the user
    def showAddGroupDialog(self, initWithFiles=[]):
        tdiag = groupEdit(self, self.lastOpenDir, colors=predefinedColors)
        if tdiag.executeDialog(groupfiles=initWithFiles, activeColor=predefinedColors[self.ui.groupsList.count()%len(predefinedColors)]) == QtGui.QDialog.Accepted:
            self.lastOpenDir = tdiag.getOpenDir()

            failed = defaultdict(list)
            pw = ProgressWrapper(1, parent=self, showIndProgress=True,
                                 indGroups={tdiag.getGroupName(): tdiag.getGroupFiles()})
            pw.show()
            pw.getCallingFunction()("max")(len(tdiag.getGroupFiles()))
            pw.getCallingFunction()("value")(0)

            i = 0
            for f in tdiag.getGroupFiles():
                pw.getCallingFunction()("value")(i)
                pw.getCallingFunction()("text")(f)
                pw.getCallingFunction()("statuscolor")(f, "orange")
                pw.getCallingFunction()("statustext")(f, text="File: %s\nStatus: %s" % (f, "importing"))

                i += 1

                x = self.checkFileImport(f)

                if x == True:
                    pw.getCallingFunction()("statuscolor")(f, "olivedrab")
                    pw.getCallingFunction()("statustext")(f, text="File: %s\nStatus: %s" % (f, "imported"))
                else:
                    failed[x].append(f)
                    pw.getCallingFunction()("statuscolor")(f, "firebrick")
                    pw.getCallingFunction()("statustext")(f, text="File: %s\nStatus: %s" % (f, "failed to import"))

            pw.hide()

            if len(failed) > 0:
                t=[]
                for q in failed.keys():
                    t.append(q)
                    t.append("\n")
                    for w in failed[q]:
                        t.append("  - ")
                        t.append(w)
                        t.append("\n")
                t="".join(t)
                QtGui.QMessageBox.warning(self, "MetExtract",
                                          "%d files failed to import correctly. Group will not be created.\n"
                                          "Please resolve the following issues (see documentation for further help)\n\n%s" % (len(failed), t), QtGui.QMessageBox.Ok)
                return None

            self.addGroup(name=tdiag.getGroupName(), files=tdiag.getGroupFiles(), minGrpFound=tdiag.getMinimumGroupFound(),
                          omitFeatures=tdiag.getOmitFeatures(), useForMetaboliteGrouping=tdiag.getUseForMetaboliteGrouping(),
                          color=str(tdiag.getGroupColor()))

            if self.updateLCMSSampleSettings():
                self.grpFileEdited = True

    # remove one group from the input
    def remGrp(self):
        todel = []
        for selectedIndex in self.ui.groupsList.selectedIndexes():
            todel.append(selectedIndex.row())

        todel.sort()
        todel.reverse()

        for ind in todel:
            self.ui.groupsList.takeItem(ind)

        self.updateLCMSSampleSettings()
        self.grpFileEdited = True

    # show group edit dialog to the user
    def editGroup(self):
        selRow = self.ui.groupsList.selectedIndexes()[0].row()

        grp = self.ui.groupsList.item(selRow).data(QListWidgetItem.UserType).toPyObject()
        t = groupEdit(colors=predefinedColors)
        if t.executeDialog(groupName=grp.name, groupfiles=grp.files, minimumGroupFound=grp.minFound,
                           omitFeatures=grp.omitFeatures, useForMetaboliteGrouping=grp.useForMetaboliteGrouping, activeColor=grp.color) == QtGui.QDialog.Accepted:

            self.addGroup(name=t.getGroupName(), files=t.getGroupFiles(), minGrpFound=t.getMinimumGroupFound(),
                          omitFeatures=t.getOmitFeatures(), useForMetaboliteGrouping=t.getUseForMetaboliteGrouping(),
                          color=str(t.getGroupColor()), atPos=selRow)

            self.ui.groupsList.takeItem(selRow+1)

            self.updateLCMSSampleSettings()
            self.grpFileEdited = True

    #</editor-fold>

    #<editor-fold desc="### load/save group compilation">
    def saveGroupsClicked(self):

        groupFile = QtGui.QFileDialog.getSaveFileName(caption="Select group file", directory=self.lastOpenDir,
                                                      filter="Group file (*.grp);;All files(*.*)")

        if groupFile is not None and len(groupFile) > 0:
            self.lastOpenDir = str(groupFile).replace("\\", "/")
            self.lastOpenDir = self.lastOpenDir[:self.lastOpenDir.rfind("/")]

            self.saveGroups(groupFile)

    # save group compilation and settings
    def saveGroups(self, groupFile=None):

        if groupFile is not None and len(groupFile) > 0:
            grps = QtCore.QSettings(groupFile, QtCore.QSettings.IniFormat)

            grps.beginGroup("ExperimentDescription")
            grps.setValue("ExperimentName", self.ui.exExperimentName_LineEdit.text())
            grps.setValue("Operator", self.ui.exOperator_LineEdit.text())
            grps.setValue("ExperimentID", self.ui.exExperimentID_LineEdit.text())
            grps.setValue("ExperimentComments", self.ui.exComments_TextEdit.toPlainText())
            grps.endGroup()

            for group in [t.data(QListWidgetItem.UserType).toPyObject() for t in natSort(self.ui.groupsList.findItems('*', QtCore.Qt.MatchWildcard), key=lambda x: str(x.data(QListWidgetItem.UserType).toPyObject().name))]:
                grps.beginGroup(group.name)
                for i in range(len(group.files)):
                    relFilePath = "./" + str(
                        os.path.relpath(group.files[i], os.path.split(str(groupFile))[0]).replace("\\", "/"))
                    grps.setValue(group.name + "__" + str(i), relFilePath)
                grps.setValue("Min_Peaks_Found", group.minFound)
                grps.setValue("OmitFeatures", group.omitFeatures)
                grps.setValue("Color", group.color)
                grps.setValue("UseForMetaboliteGrouping", group.useForMetaboliteGrouping)
                grps.endGroup()

            grps.beginGroup("ExperimentResults")

            resFile = str(self.ui.groupsSave.text().replace("\\", "/"))
            if resFile.startswith("./"):
                pat=os.path.split(str(groupFile))[0] + "/"+resFile
                grps.setValue("GroupSaveFile", resFile)
                self.ui.groupsSave.setText(pat)
            else:
                grps.setValue("GroupSaveFile", resFile)
                self.ui.groupsSave.setText(resFile)

            grps.endGroup()

            # ask user, if currently loaded settings shall also be saved to this compilation
            if QtGui.QMessageBox.question(self, "MetExtract", "Do you want to save the currently loaded settings with this group?",
                                          QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.Yes:
                self.saveSettingsFile(groupFile, clear=False)


    def loadGroupsClicked(self):

        groupFile = QtGui.QFileDialog.getOpenFileName(caption="Select group file", directory=self.lastOpenDir,
                                                      filter="Group file (*.grp);;All files(*.*)")

        if groupFile is not None and len(groupFile) > 0:
            self.lastOpenDir = str(groupFile).replace("\\", "/")
            self.lastOpenDir = self.lastOpenDir[:self.lastOpenDir.rfind("/")]

            self.loadGroups(groupFile)

    def loadGroups(self, groupFile=None, forceLoadSettings=False, askLoadSettings=True):

        if groupFile is not None and len(groupFile) > 0:

            if len(groupFile) > 0:
                self.lastOpenDir = str(groupFile).replace("\\", "/")
                self.lastOpenDir = self.lastOpenDir[:self.lastOpenDir.rfind("/")]

            grps = QtCore.QSettings(groupFile, QtCore.QSettings.IniFormat)

            grps.beginGroup("ExperimentDescription")

            if grps.contains("ExperimentName"):
                self.ui.exExperimentName_LineEdit.setText(str(grps.value("ExperimentName").toString()))
            if grps.contains("Operator"):
                self.ui.exOperator_LineEdit.setText(str(grps.value("Operator").toString()))
            if grps.contains("ExperimentID"):
                self.ui.exExperimentID_LineEdit.setText(str(grps.value("ExperimentID").toString()))
            if grps.contains("ExperimentComments"):
                self.ui.exComments_TextEdit.setPlainText(str(grps.value("ExperimentComments").toString()))

            grps.endGroup()

            procGrps=[]
            for grp in grps.childGroups():
                if str(grp) != "Settings" and str(grp) != "ExperimentDescription" and str(grp) != "ExperimentResults" and str(grp) != "MetExtract":
                    procGrps.append(grp)

            grpi=0
            groupsToAdd=[]
            for grp in natSort(procGrps):
                grps.beginGroup(grp)
                kids = []
                minFound = 1
                color = predefinedColors[grpi%len(predefinedColors)]
                omitFeatures = True
                useForMetaboliteGrouping = True
                for kid in grps.childKeys():
                    if str(kid) == "Min_Peaks_Found":
                        minFound = grps.value(kid).toInt()[0]
                    elif str(kid) == "OmitFeatures":
                        omitFeatures = grps.value(kid).toBool()
                    elif str(kid) == "Color":
                        color= str(grps.value(kid).toString())
                    elif str(kid) == "UseForMetaboliteGrouping":
                        useForMetaboliteGrouping = grps.value(kid).toBool()
                    elif str(kid).startswith(grp):
                        if os.path.isabs(str(grps.value(kid).toString()).replace("\\", "/")):
                            kids.append(str(grps.value(kid).toString()).replace("\\", "/"))
                        else:
                            kids.append(os.path.split(str(groupFile))[0].replace("\\", "/") + "/" + str(
                                grps.value(kid).toString()).replace("\\", "/"))

                groupsToAdd.append(Bunch(name=grp, files=kids, minGrpFound=minFound, omitFeatures=omitFeatures, useForMetaboliteGrouping=useForMetaboliteGrouping, color=color))

                grps.endGroup()
                grpi += 1

            for grpToAdd in natSort(groupsToAdd, key=lambda x:x.name):
                self.addGroup(name=grpToAdd.name, files=grpToAdd.files, minGrpFound=grpToAdd.minGrpFound, omitFeatures=grpToAdd.omitFeatures, useForMetaboliteGrouping=grpToAdd.useForMetaboliteGrouping, color=grpToAdd.color)

            grps.beginGroup("ExperimentResults")
            if grps.contains("GroupSaveFile"):
                resFile = str(grps.value("GroupSaveFile").toString())

                if resFile.startswith("./"):
                    resFile = os.path.split(str(groupFile))[0] + "/" + resFile

                self.ui.groupsSave.setText(resFile)
                self.loadGroupsResultsFile(resFile)

            grps.endGroup()

            self.updateLCMSSampleSettings()

            if forceLoadSettings or (askLoadSettings and QtGui.QMessageBox.question(self, "MetExtract","Do you want to load the associated settings with this group?",
                                                                                       QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.Yes):
                self.loadSettingsFile(str(groupFile))
            self.grpFile = str(groupFile)
            self.grpFileEdited = False
    #</editor-fold>

    #<editor-fold desc="### group results visualisation functions">
    # EXPERIMENTAL
    def loadGroupsResultsFile(self, groupsResFile):
        try:
            if os.path.exists(groupsResFile+".identified.sqlite") and os.path.isfile(groupsResFile+".identified.sqlite"):

                self.ui.resultsExperiment_TreeWidget.clear()
                self.ui.resultsExperiment_TreeWidget.setHeaderLabels(["OGroup", "MZ", "RT", "Xn", "Z", "IonMode"])
                widths=[80,90,90,90,90,90]
                for i in range(len(widths)):
                    self.ui.resultsExperiment_TreeWidget.setColumnWidth(i, widths[i])

                self.experimentResults=Bunch(conn=None, curs=None, file=None)

                self.experimentResults.conn=connect(groupsResFile+".identified.sqlite")
                self.experimentResults.curs=self.experimentResults.conn.cursor()

                metaboliteGroupTreeItems={}
                for metaboliteGroup in SQLSelectAsObject(self.experimentResults.curs, "SELECT DISTINCT 'metaboliteGroup' AS type, OGroup AS metaboliteGroupID FROM GroupResults ORDER BY rt"):
                    metaboliteGroupTreeItem=QtGui.QTreeWidgetItem(["%d"%metaboliteGroup.metaboliteGroupID])
                    metaboliteGroupTreeItem.bunchData=metaboliteGroup
                    self.ui.resultsExperiment_TreeWidget.addTopLevelItem(metaboliteGroupTreeItem)
                    metaboliteGroupTreeItems[metaboliteGroup.metaboliteGroupID]=metaboliteGroupTreeItem

                for fp in SQLSelectAsObject(self.experimentResults.curs, "SELECT 'featurePair' AS type, id, OGroup AS metaboliteGroupID, mz, lmz, dmz, rt, xn, charge, scanEvent, ionisationMode, tracer FROM GroupResults ORDER BY mz"):
                    featurePair=QtGui.QTreeWidgetItem([str(fp.id), "%.4f"%fp.mz, "%.2f"%(fp.rt/60.), str(fp.xn), str(fp.charge), str(fp.ionisationMode)])
                    featurePair.bunchData=fp
                    metaboliteGroupTreeItems[fp.metaboliteGroupID].addChild(featurePair)

                for grpID in metaboliteGroupTreeItems.keys():
                    kids=[]
                    for i in range(metaboliteGroupTreeItems[grpID].childCount()):
                        kids.append(metaboliteGroupTreeItems[grpID].child(i))
                    meanRT=mean([float(kid.text(2)) for kid in kids])
                    metaboliteGroupTreeItems[grpID].setText(1, "%d"%len(kids))
                    metaboliteGroupTreeItems[grpID].setText(2, "%.2f"%meanRT)

                ##TODO finish that
        except Exception as e:
            logging.error("Multiple file results could not be fetched correctly")
            pass




    def closeLoadedGroupsResultsFile(self):
        if hasattr(self, "experimentResults"):
            self.ui.resultsExperiment_TreeWidget.clear()
            self.experimentResults.curs.close()
            self.experimentResults.conn.close()
            delattr(self, "experimentResults")


    def resultsExperimentChanged(self):
        self.clearPlot(self.ui.resultsExperiment_plot)
        self.clearPlot(self.ui.resultsExperimentSeparatedPeaks_plot)
        self.clearPlot(self.ui.resultsExperimentMSScanPeaks_plot)

        plotItems=[]

        for item in self.ui.resultsExperiment_TreeWidget.selectedItems():
            if item.bunchData.type=="metaboliteGroup":
                for i in range(item.childCount()):
                    plotItems.append(item.child(i))
            if item.bunchData.type=="featurePair":
                plotItems.append(item)

        for item in plotItems:
            assert item.bunchData.type=="featurePair"
            self.ui.expRes_ID_LineEdit.setText(str(item.bunchData.id))
            self.ui.expRes_MZ_Spinner.setValue(item.bunchData.mz)
            rt=item.bunchData.rt/60.
            self.ui.expRes_RT_Spinner.setValue(rt)
            self.ui.expRes_Z_Spinner.setValue(item.bunchData.charge)
            self.ui.expRes_IonMode_LineEdit.setText(item.bunchData.ionisationMode)

            groups={}
            for group in SQLSelectAsObject(self.experimentResults.curs, "SELECT groupName, id FROM FileGroups"):
                groups[group.id]=group.groupName

            filesToGroup={}
            filesInGroups={}

            for fileMapping in SQLSelectAsObject(self.experimentResults.curs, "SELECT fileName, filePath, groupID FROM FileMapping"):
                filesToGroup[fileMapping.fileName]=fileMapping.groupID
                if not filesInGroups.has_key(fileMapping.groupID):
                    filesInGroups[fileMapping.groupID]=[]
                filesInGroups[fileMapping.groupID].append(fileMapping)

            foundIn={}
            for foundFP in SQLSelectAsObject(self.experimentResults.curs, "SELECT file AS fil, featurePairID, featureGroupID, areaN, areaL, featureType FROM FoundFeaturePairs WHERE resID=%d"%item.bunchData.id):
                foundIn[foundFP.fil]=foundFP

            offsetCount=0
            for groupID, groupName in groups.items():
                for file in filesInGroups[groupID]:
                    if file.fileName in foundIn.keys():

                        conn=connect(file.filePath+".identified.sqlite")
                        curs=conn.cursor()

                        groupID=filesToGroup[file.fileName]

                        msSpectrumID=None

                        for XICObj in SQLSelectAsObject(curs, "SELECT x.xic, x.xicL, x.times, c.massSpectrumID FROM XICs x INNER JOIN chromPeaks c ON x.id=c.eicID WHERE c.id=%d"%foundIn[file.fileName].featurePairID):
                            XICObj.xic=[float(f) for f in XICObj.xic.split(";")]
                            XICObj.xicL=[float(f) for f in XICObj.xicL.split(";")]
                            XICObj.times=[float(f) for f in XICObj.times.split(";")]

                            msSpectrumID=XICObj.massSpectrumID

                            useInds=[]
                            bestCenter=0
                            bestCenterDiff=1000000
                            timeWindow=4*60./2
                            for i, t in enumerate(XICObj.times):
                                if (item.bunchData.rt-timeWindow)<t<(item.bunchData.rt+timeWindow):
                                    useInds.append(i)
                                if abs(t-item.bunchData.rt)<bestCenterDiff:
                                    bestCenterDiff=abs(t-item.bunchData.rt)
                                    bestCenter=i
                            centerInt=1
                            if self.ui.resultsExperimentNormaliseXICs_checkBox.checkState()==QtCore.Qt.Checked:
                                centerInt=XICObj.xic[bestCenter]

                            self.ui.resultsExperiment_plot.axes.plot([t/60. for t in XICObj.times], [f/centerInt for f in XICObj.xic],   color=predefinedColors[groupID%len(predefinedColors)])
                            self.ui.resultsExperiment_plot.axes.plot([t/60. for t in XICObj.times], [-f/centerInt for f in XICObj.xicL], color=predefinedColors[groupID%len(predefinedColors)])
                            self.ui.resultsExperiment_plot.axes.set_xlim([rt-.5, rt+.5])


                            minInd=min(useInds)
                            maxInd=max(useInds)

                            self.ui.resultsExperimentSeparatedPeaks_plot.axes.plot([t/60. +offsetCount*0.33 for t in XICObj.times[minInd:maxInd]], [f/centerInt for f in XICObj.xic[minInd:maxInd]],   color=predefinedColors[groupID%len(predefinedColors)])
                            self.ui.resultsExperimentSeparatedPeaks_plot.axes.plot([t/60. +offsetCount*0.33 for t in XICObj.times[minInd:maxInd]], [-f/centerInt for f in XICObj.xicL[minInd:maxInd]], color=predefinedColors[groupID%len(predefinedColors)])

                        if msSpectrumID is not None:
                            for msSpectrum in SQLSelectAsObject(curs, "SELECT mzs, intensities, ionMode FROM massspectrum WHERE mID=%d"%msSpectrumID):
                                mzs=[float(f) for f in msSpectrum.mzs.split(";")]
                                intensities=[float(f) for f in msSpectrum.intensities.split(";")]

                                self.ui.resultsExperimentMSScanPeaks_plot.axes.stem(mzs, intensities, lineEdgeColor=predefinedColors[groupID%len(predefinedColors)],color=predefinedColors[groupID%len(predefinedColors)])


                        curs.close()
                        conn.close()

                offsetCount = offsetCount+1

        self.drawCanvas(self.ui.resultsExperiment_plot)
        self.drawCanvas(self.ui.resultsExperimentSeparatedPeaks_plot)
        self.drawCanvas(self.ui.resultsExperimentMSScanPeaks_plot)
    #</editor-fold>

    #<editor-fold desc="### load/save settings functions">
    def loadSettings(self):
        settingsFile = QtGui.QFileDialog.getOpenFileName(caption="Select settings file", directory=self.lastOpenDir, filter="Setting file (*.ini);;Group file with settings (*.grp);;All files(*.*)")
        if len(settingsFile) > 0:
            self.lastOpenDir = str(settingsFile).replace("\\", "/")
            self.lastOpenDir = self.lastOpenDir[:self.lastOpenDir.rfind("/")]
        self.loadSettingsFile(settingsFile)

    # load settings from key-value pair file
    def loadSettingsFile(self, settingsFile, checkExperimentType=True, silent=False):
        if len(settingsFile) > 0:
            sett = QtCore.QSettings(settingsFile, QtCore.QSettings.IniFormat)

            sett.beginGroup("MetExtract")

            if sett.contains("Version"):
                version=str(sett.value("Version").toString())
                if not(silent) and version != MetExtractVersion and QtGui.QMessageBox.question(self, "MetExtract", "Settings were created from a different MetExtract II version. Some settings may not be imported correctly."
                                                                      "\nDo you want to continue?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.No:
                    return

            if sett.contains("R-Version"):
                rversion=str(sett.value("R-Version").toString())
                if not(silent) and rversion != getRVersion() and QtGui.QMessageBox.question(self, "MetExtract", "Settings were created from a different R Version. Processed results may slightly be different."
                                                                      "\nDo you want to continue?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.No:
                    return

            sett.endGroup()

            sett.beginGroup("Settings")

            if sett.contains("processIndividualFiles"):
                self.ui.processIndividualFiles.setChecked(sett.value("processIndividualFiles").toBool())

            if checkExperimentType:
                if sett.contains("TracerExperiment") :
                    if sett.value("TracerExperiment").toBool() and self.labellingExperiment==METABOLOME:
                        QtGui.QMessageBox.warning(self, "MetExtract", "Warning: You are trying to open a TracExtract experiment "
                                                                      "with AllExtract. Labelling parameters will not be loaded. "
                                                                      "Please switch to TracExtract", QtGui.QMessageBox.Ok)
                    elif not(sett.value("TracerExperiment").toBool()) and self.labellingExperiment==TRACER:
                        QtGui.QMessageBox.warning(self, "MetExtract", "Warning: You are trying to open an AllExtract experiment "
                                                                      "with TracExtract. Labelling parameters will not be loaded. "
                                                                      "Please switch to AllExtract", QtGui.QMessageBox.Ok)
                else:
                    QtGui.QMessageBox.warning(self, "MetExtract", "Warning: Experiment type is not stored in the settings. Please ensure "
                                                                  "that the correct module has been loaded for your experiment",
                                              QtGui.QMessageBox.Ok)

            if sett.contains("tracerConfiguration") and self.labellingExperiment==TRACER:
                self.configuredTracers = loads(str(base64.b64decode(sett.value("tracerConfiguration").toString())))

            if sett.contains("LabellingElementA") and self.labellingExperiment==METABOLOME:
                self.ui.isotopeAText.setText(str(sett.value("LabellingElementA").toString()))
            if sett.contains("IsotopicAbundanceA") and self.labellingExperiment==METABOLOME:
                self.ui.isotopicAbundanceA.setValue(sett.value("IsotopicAbundanceA").toDouble()[0])
            if sett.contains("LabellingElementB") and self.labellingExperiment==METABOLOME:
                self.ui.isotopeBText.setText(str(sett.value("LabellingElementB").toString()))
            if sett.contains("IsotopicAbundanceB") and self.labellingExperiment==METABOLOME:
                self.ui.isotopicAbundanceB.setValue(sett.value("IsotopicAbundanceB").toDouble()[0])
            if sett.contains("useCValidation") and self.labellingExperiment==METABOLOME:
                self.ui.useCValidation.setCheckState(sett.value("useCValidation").toInt()[0])
            if sett.contains("useRatio") and self.labellingExperiment==METABOLOME:
                self.ui.useRatio.setChecked(sett.value("useRatio").toBool())
            if sett.contains("minRatio") and self.labellingExperiment==METABOLOME:
                self.ui.minRatio.setValue(sett.value("minRatio").toDouble()[0])
            if sett.contains("maxRatio") and self.labellingExperiment==METABOLOME:
                self.ui.maxRatio.setValue(sett.value("maxRatio").toDouble()[0])

            if sett.contains("IntensityThreshold"):
                self.ui.intensityThreshold.setValue(sett.value("IntensityThreshold").toInt()[0])
            if sett.contains("IntensityCutoff"):
                self.ui.intensityCutoff.setValue(sett.value("IntensityCutoff").toInt()[0])
            if sett.contains("intensityThresholdIsotopologs"):
                self.ui.intensityThresholdIsotopologs.setValue(sett.value("intensityThresholdIsotopologs").toInt()[0])

            if sett.contains("ScanStart"):
                self.ui.scanStartTime.setValue(sett.value("ScanStart").toDouble()[0])
            if sett.contains("ScanEnd"):
                self.ui.scanEndTime.setValue(sett.value("ScanEnd").toDouble()[0])
            if sett.contains("MinXCount"):
                self.ui.minXCount.setValue(sett.value("MinXCount").toInt()[0])
            if sett.contains("MaxXCount"):
                self.ui.maxXCount.setValue(sett.value("MaxXCount").toInt()[0])
            if sett.contains("MaxLoading"):
                self.ui.maxLoading.setValue(sett.value("MaxLoading").toInt()[0])
            if sett.contains("MaxMassDeviation"):
                self.ui.ppmRangeIdentification.setValue(sett.value("MaxMassDeviation").toDouble()[0])
            if sett.contains("IsotopicPatternCountA"):
                self.ui.isotopePatternCountA.setValue(sett.value("IsotopicPatternCountA").toInt()[0])
            if sett.contains("IsotopicPatternCountB"):
                self.ui.isotopePatternCountB.setValue(sett.value("IsotopicPatternCountB").toInt()[0])
            if sett.contains("lowAbundanceIsotopeCutoff"):
                self.ui.isoAbundance.setChecked(sett.value("lowAbundanceIsotopeCutoff").toBool())
            if sett.contains("IntensityAbundanceErrorA"):
                self.ui.baseRange.setValue(sett.value("IntensityAbundanceErrorA").toDouble()[0])
            if sett.contains("IntensityAbundanceErrorB"):
                self.ui.isotopeRange.setValue(sett.value("IntensityAbundanceErrorB").toDouble()[0])

            if sett.contains("ClustPPM"):
                self.ui.clustPPM.setValue(sett.value("ClustPPM").toDouble()[0])
            if sett.contains("minSpectraCount"):
                self.ui.minSpectraCount.setValue(sett.value("minSpectraCount").toInt()[0])
            if sett.contains("correctCCount"):
                self.ui.correctcCount.setChecked(sett.value("correctcCount").toBool())

            if sett.contains("EICppm"):
                self.ui.wavelet_EICppm.setValue(sett.value("EICppm").toDouble()[0])
            if sett.contains("EICSmoothingWindow"):
                eicSmoothingWindow = str(sett.value("EICSmoothingWindow").toString())
                pos = self.ui.eicSmoothingWindow.findText(eicSmoothingWindow)
                if pos >= 0:
                    self.ui.eicSmoothingWindow.setCurrentIndex(pos)
                else:
                    QtGui.QMessageBox.information(self, "MetExtract", "Invalid EIC smoothing option",
                                                  QtGui.QMessageBox.Ok);
                    self.ui.eicSmoothingWindow.setCurrentIndex(0)
            if sett.contains("EICSmoothingWindowSize"):
                self.ui.eicSmoothingWindowSize.setValue(sett.value("EICSmoothingWindowSize").toInt()[0])
            if sett.contains("Wavelet_MinScale"):
                self.ui.wavelet_minScale.setValue(sett.value("Wavelet_MinScale").toInt()[0])
            if sett.contains("Wavelet_MaxScale"):
                self.ui.wavelet_maxScale.setValue(sett.value("Wavelet_MaxScale").toInt()[0])
            if sett.contains("Wavelet_SNRTh"):
                self.ui.wavelet_SNRThreshold.setValue(sett.value("Wavelet_SNRTh").toInt()[0])
            if sett.contains("Peak_CenterError"):
                self.ui.peak_centerError.setValue(sett.value("Peak_CenterError").toInt()[0])
            if sett.contains("Peak_WidthError"):
                self.ui.peak_scaleError.setValue(sett.value("Peak_WidthError").toInt()[0])
            if sett.contains("Peak_minPeakCorr"):
                self.ui.minPeakCorr.setValue(sett.value("Peak_minPeakCorr").toDouble()[0])
            if sett.contains("checkBox_checkPeakRatio"):
                self.ui.checkBox_checkPeakRatio.setChecked(sett.value("checkBox_checkPeakRatio").toBool())
            if sett.contains("doubleSpinBox_minPeakRatio"):
                self.ui.doubleSpinBox_minPeakRatio.setValue(sett.value("doubleSpinBox_minPeakRatio").toDouble()[0])
            if sett.contains("doubleSpinBox_maxPeakRatio"):
                self.ui.doubleSpinBox_maxPeakRatio.setValue(sett.value("doubleSpinBox_maxPeakRatio").toDouble()[0])

            if sett.contains("calcIsoRatioNative"):
                self.ui.calcIsoRatioNative_spinBox.setValue(sett.value("calcIsoRatioNative").toInt()[0])
            if sett.contains("calcIsoRatioLabelled"):
                self.ui.calcIsoRatioLabelled_spinBox.setValue(sett.value("calcIsoRatioLabelled").toInt()[0])
            if sett.contains("calcIsoRatioMoiety"):
                self.ui.calcIsoRatioMoiety_spinBox.setValue(sett.value("calcIsoRatioMoiety").toInt()[0])

            if sett.contains("hAIntensityError"):
                self.ui.hAIntensityError.setValue(sett.value("hAIntensityError").toDouble()[0])
            if sett.contains("hAMinScans"):
                self.ui.hAMinScans.setValue(sett.value("hAMinScans").toInt()[0])
            if sett.contains("heteroAtoms") or sett.contains("heteroElements"):
                if sett.contains("heteroAtoms"):
                    self.heteroElements = loads(str(base64.b64decode(sett.value("heteroAtoms").toString())))
                elif sett.contains("heteroElements"):
                    self.heteroElements = loads(str(base64.b64decode(sett.value("heteroElements").toString())))
                if any([not(isinstance(d, ConfiguredHeteroAtom)) for d in self.heteroElements]):
                    QtGui.QMessageBox.warning(self, "MetExtract", "Could not parse hetero elements. Predefined elements will be used instead. "
                                                                  "Please adapt the list accordingly and save the compilation again.")
                    self.heteroElements=self.preConfigured_heteroElements

            if sett.contains("saveTSV"):
                self.ui.saveCSV.setChecked(sett.value("saveTSV").toBool())
            if sett.contains("saveMZXML"):
                self.ui.saveMZXML.setChecked(sett.value("saveMZXML").toBool())
            if sett.contains("writeMZXMLOptions"):
                writeMZXMLOptions = sett.value("writeMZXMLOptions").toInt()[0]
                self.ui.wm_ia.setChecked(writeMZXMLOptions & 1)
                self.ui.wm_iap.setChecked(writeMZXMLOptions & 2)
                self.ui.wm_imb.setChecked(writeMZXMLOptions & 4)
                self.ui.wm_ib.setChecked(writeMZXMLOptions & 8)
            if sett.contains("savePDF"):
                self.ui.savePDF.setChecked(sett.value("savePDF").toBool())

            if sett.contains("minCorrelation"):
                self.ui.minCorrelation.setValue(sett.value("minCorrelation").toDouble()[0])
            if sett.contains("minCorrelationConnections"):
                self.ui.minCorrelationConnections.setValue(sett.value("minCorrelationConnections").toDouble()[0])
            if sett.contains("adducts"):
                self.adducts = loads(str(base64.b64decode(sett.value("adducts").toString())))
                if any([not(isinstance(d, ConfiguredAdduct)) for d in self.adducts]):
                    QtGui.QMessageBox.warning(self, "MetExtract", "Could not parse adducts for ion annotation. Predefined adducts will be used instead. "
                                                                  "Please adapt the list accordingly and save the compilation again.")
                    self.adducts=self.preConfigured_adducts

            if sett.contains("elements") or sett.contains("elementsForNL"):
                if sett.contains("elements"):
                    self.elementsForNL = loads(str(base64.b64decode(sett.value("elements").toString())))
                elif sett.contains("elementsForNL"):
                    self.elementsForNL = loads(str(base64.b64decode(sett.value("elementsForNL").toString())))
                if any([not(isinstance(d, ConfiguredElement)) for d in self.elementsForNL]):
                    QtGui.QMessageBox.warning(self, "MetExtract", "Could not parse elements for neutral loss calculation. Predefined elements will be used instead. "
                                                                  "Please adapt the list accordingly and save the compilation again.")
                    self.elementsForNL=self.preConfigured_elementsForNL

            if sett.contains("processMultipleFiles"):
                self.ui.processMultipleFiles.setChecked(sett.value("processMultipleFiles").toBool())

            if sett.contains("GroupFeatures"):
                self.ui.groupResults.setChecked(sett.value("GroupFeatures").toBool())
            if sett.contains("GroupPPM"):
                self.ui.groupPpm.setValue(sett.value("GroupPPM").toDouble()[0])
            if sett.contains("GroupAlign"):
                self.ui.alignChromatograms.setChecked(sett.value("GroupAlign").toBool())
            if sett.contains("GroupNPolynom"):
                self.ui.polynomValue.setValue(sett.value("GroupNPolynom").toInt()[0])
            if sett.contains("GroupTimeWindow"):
                self.ui.groupingRT.setValue(sett.value("GroupTimeWindow").toDouble()[0])

            if sett.contains("MetaboliteClusterMinConnections"):
                self.ui.metaboliteClusterMinConnections.setValue(sett.value("MetaboliteClusterMinConnections").toInt()[0])
            if sett.contains("minConnectionRate"):
                self.ui.minConnectionRate.setValue(sett.value("minConnectionRate").toDouble()[0])

            if sett.contains("GroupIntegrateMissedPeaks"):
                self.ui.integratedMissedPeaks.setChecked(sett.value("GroupIntegrateMissedPeaks").toBool())
            if sett.contains("integrationMaxTimeDifference"):
                self.ui.integrationMaxTimeDifference.setValue(sett.value("integrationMaxTimeDifference").toDouble()[0])
            if sett.contains("reintegrateIntensityCutoff"):
                self.ui.reintegrateIntensityCutoff.setValue(sett.value("reintegrateIntensityCutoff").toDouble()[0])

            sett.endGroup()

            self.updateTracerInfo()

    def saveSettings(self):
        settingsFile = QtGui.QFileDialog.getSaveFileName(caption="Select settings file", directory=self.lastOpenDir,
                                                         filter="Setting file (*.ini);;All files (*.*)")
        if len(settingsFile) > 0:
            self.lastOpenDir = str(settingsFile).replace("\\", "/")
            self.lastOpenDir = self.lastOpenDir[:self.lastOpenDir.rfind("/")]
        if len(settingsFile) > 0:
            self.lastOpenDir = str(settingsFile).replace("\\", "/")
            self.lastOpenDir = self.lastOpenDir[:self.lastOpenDir.rfind("/")]
        self.saveSettingsFile(settingsFile)

    # save currently loaded settings to key-value pair file
    def saveSettingsFile(self, settingsFile, clear=True):
        if len(settingsFile) > 0:
            sett = QtCore.QSettings(settingsFile, QtCore.QSettings.IniFormat)
            if clear:
                sett.clear()
            sett.beginGroup("Settings")

            sett.setValue("processIndividualFiles", self.ui.processIndividualFiles.checkState() == QtCore.Qt.Checked)

            sett.setValue("TracerExperiment", self.labellingExperiment==TRACER)

            if self.labellingExperiment==TRACER:
                sett.setValue("tracerConfiguration", base64.b64encode(dumps(self.configuredTracers)))
            else:
                sett.setValue("LabellingElementA", self.ui.isotopeAText.text())
                sett.setValue("IsotopicAbundanceA", self.ui.isotopicAbundanceA.value())
                sett.setValue("LabellingElementB", self.ui.isotopeBText.text())
                sett.setValue("IsotopicAbundanceB", self.ui.isotopicAbundanceB.value())
                sett.setValue("useCValidation", str(self.ui.useCValidation.checkState()))
                sett.setValue("useRatio", self.ui.useRatio.checkState() == QtCore.Qt.Checked)
                sett.setValue("minRatio", self.ui.minRatio.value())
                sett.setValue("maxRatio", self.ui.maxRatio.value())

            sett.setValue("IntensityThreshold", self.ui.intensityThreshold.value())
            sett.setValue("IntensityCutoff", self.ui.intensityCutoff.value())
            sett.setValue("ScanStart", self.ui.scanStartTime.value())
            sett.setValue("ScanEnd", self.ui.scanEndTime.value())
            sett.setValue("MinXCount", self.ui.minXCount.value())
            sett.setValue("MaxXCount", self.ui.maxXCount.value())
            sett.setValue("MaxLoading", self.ui.maxLoading.value())
            sett.setValue("MaxMassDeviation", self.ui.ppmRangeIdentification.value())
            sett.setValue("IsotopicPatternCountA", self.ui.isotopePatternCountA.value())
            sett.setValue("IsotopicPatternCountB", self.ui.isotopePatternCountB.value())
            sett.setValue("lowAbundanceIsotopeCutoff", self.ui.isoAbundance.checkState() == QtCore.Qt.Checked)
            sett.setValue("intensityThresholdIsotopologs", self.ui.intensityThresholdIsotopologs.value())
            sett.setValue("IntensityAbundanceErrorA", self.ui.baseRange.value())
            sett.setValue("IntensityAbundanceErrorB", self.ui.isotopeRange.value())

            sett.setValue("ClustPPM", self.ui.clustPPM.value())
            sett.setValue("minSpectraCount", self.ui.minSpectraCount.value())
            sett.setValue("correctCCount", self.ui.correctcCount.checkState() == QtCore.Qt.Checked)
            sett.setValue("EICppm", self.ui.wavelet_EICppm.value())
            sett.setValue("EICSmoothingWindow", str(self.ui.eicSmoothingWindow.currentText())),
            sett.setValue("EICSmoothingWindowSize", self.ui.eicSmoothingWindowSize.value()),
            sett.setValue("Wavelet_MinScale", self.ui.wavelet_minScale.value())
            sett.setValue("Wavelet_MaxScale", self.ui.wavelet_maxScale.value())
            sett.setValue("Wavelet_SNRTh", self.ui.wavelet_SNRThreshold.value())
            sett.setValue("Peak_CenterError", self.ui.peak_centerError.value())
            sett.setValue("Peak_WidthError", self.ui.peak_scaleError.value())
            sett.setValue("Peak_minPeakCorr", self.ui.minPeakCorr.value())
            sett.setValue("checkBox_checkPeakRatio", self.ui.checkBox_checkPeakRatio.isChecked())
            sett.setValue("doubleSpinBox_minPeakRatio", self.ui.doubleSpinBox_minPeakRatio.value())
            sett.setValue("doubleSpinBox_maxPeakRatio", self.ui.doubleSpinBox_maxPeakRatio.value())

            sett.setValue("calcIsoRatioNative", self.ui.calcIsoRatioNative_spinBox.value())
            sett.setValue("calcIsoRatioLabelled", self.ui.calcIsoRatioLabelled_spinBox.value())
            sett.setValue("calcIsoRatioMoiety", self.ui.calcIsoRatioMoiety_spinBox.value())

            sett.setValue("hAIntensityError", self.ui.hAIntensityError.value())
            sett.setValue("hAMinScans", self.ui.hAMinScans.value())
            sett.setValue("heteroElements", base64.b64encode(dumps(self.heteroElements)))

            sett.setValue("saveTSV", self.ui.saveCSV.checkState() == QtCore.Qt.Checked)
            sett.setValue("saveMZXML", self.ui.saveMZXML.isChecked())
            writeMZXMLOptions = 0
            if self.ui.wm_ia.checkState() == QtCore.Qt.Checked:
                writeMZXMLOptions |= 1
            if self.ui.wm_iap.checkState() == QtCore.Qt.Checked:
                writeMZXMLOptions |= 2
            if self.ui.wm_imb.checkState() == QtCore.Qt.Checked:
                writeMZXMLOptions |= 4
            if self.ui.wm_ib.checkState() == QtCore.Qt.Checked:
                writeMZXMLOptions |= 8
            sett.setValue("writeMZXMLOptions", writeMZXMLOptions)

            sett.setValue("savePDF", self.ui.savePDF.checkState() == QtCore.Qt.Checked)

            sett.setValue("minCorrelation", self.ui.minCorrelation.value())
            sett.setValue("minCorrelationConnections", self.ui.minCorrelationConnections.value())
            sett.setValue("adducts", base64.b64encode(dumps(self.adducts)))
            sett.setValue("elementsForNL", base64.b64encode(dumps(self.elementsForNL)))

            sett.setValue("processMultipleFiles", self.ui.processMultipleFiles.checkState() == QtCore.Qt.Checked)

            sett.setValue("GroupFeatures", self.ui.groupResults.isChecked())
            sett.setValue("GroupPPM", self.ui.groupPpm.value())
            sett.setValue("GroupAlign", self.ui.alignChromatograms.checkState() == QtCore.Qt.Checked)
            sett.setValue("GroupNPolynom", self.ui.polynomValue.value())
            sett.setValue("GroupTimeWindow", self.ui.groupingRT.value())

            sett.setValue("MetaboliteClusterMinConnections", self.ui.metaboliteClusterMinConnections.value())
            sett.setValue("minConnectionRate", self.ui.minConnectionRate.value())

            sett.setValue("GroupIntegrateMissedPeaks", self.ui.integratedMissedPeaks.isChecked())
            sett.setValue("integrationMaxTimeDifference", self.ui.integrationMaxTimeDifference.value())
            sett.setValue("reintegrateIntensityCutoff", self.ui.reintegrateIntensityCutoff.value())

            sett.endGroup()

            sett.beginGroup("MetExtract")

            sett.setValue("Version",MetExtractVersion)
            sett.setValue("RVersion", getRVersion())

            sett.endGroup()
    #</editor-fold>

    def updateGroupPPM(self):
        self.ui.groupPpm.setValue(self.ui.ppmRangeIdentification.value() * 2)

    def selectGroupsFile(self):
        grpFile = QtGui.QFileDialog.getSaveFileName(caption="Select groups file", directory=self.lastOpenDir,
                                                    filter="TSV file (*.tsv);;CSV file (*.csv);;All files (*.*)")
        if len(grpFile) > 0:
            self.lastOpenDir = str(grpFile).replace("\\", "/")
            self.lastOpenDir = self.lastOpenDir[:self.lastOpenDir.rfind("/")]
        if len(grpFile) > 0:
            self.lastOpenDir = str(grpFile).replace("\\", "/")
            self.lastOpenDir = self.lastOpenDir[:self.lastOpenDir.rfind("/")]

        if len(grpFile) > 0:
            self.ui.groupsSave.setText(grpFile)
        else:
            self.ui.groupsSave.setText("./results.tsv")

    def exitMe(self):
        QtGui.QApplication.exit()

    def aboutMe(self):
        lic="THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR "\
            "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, "\
            "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE "\
            "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER "\
            "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, "\
            "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN "\
            "THE SOFTWARE."
        QtGui.QMessageBox.information(self, "MetExtract",
                                      "MetExtract %s\n\n(c) Centre for Analytical Chemistry, IFA Tulln\nUniversity of Natural Resources and Life Sciences, Vienna\n\n%s" % (MetExtractVersion, lic),
                                      QtGui.QMessageBox.Ok)

    # open local MetExtract documentation
    # taken from http://stackoverflow.com/questions/4216985/call-to-operating-system-to-open-url
    def helpMe(self):
        import subprocess
        import webbrowser
        import sys

        url=get_main_dir()+"/documentation/index.html"
        if sys.platform == "darwin":    # in case of OS X
            subprocess.Popen(['open', url])
        else:
            webbrowser.open_new_tab(url)



    # helper method for multiprocessing module of LC-HRMS file processing
    def runProcess(self, dontSave=False, askStarting=True):

        self.terminateJobs=False

        self.closeCurrentOpenResultsFile()
        self.closeLoadedGroupsResultsFile()

        # check certain requirements for LC-HRMS file processing
        if str(self.ui.negativeScanEvent.currentText()) == "None" and str(
                self.ui.positiveScanEvent.currentText()) == "None":
            QtGui.QMessageBox.information(self, "MetExtract",
                                          "Please select at least one scan event prior to calculations",
                                          QtGui.QMessageBox.Ok)
            return 1

        if str(self.ui.negativeScanEvent.currentText()) == str(self.ui.positiveScanEvent.currentText()):
            QtGui.QMessageBox.information(self, "MetExtract",
                                          "It seems two identical scan events are selected for the positive and negative mode. "
                                          "Select only one mode if just a single ionisation mode is present in your files.")
            return 1

        if askStarting:
            if self.ui.isotopePatternCountA.value() == 1 and self.ui.isotopePatternCountB.value() == 1:
                if QtGui.QMessageBox.question(self, "MetExtract",
                                              "Are you sure you want to search for metabolites without using isotopic peaks?",
                                              QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.No:
                    return 1

        if self.ui.intensityCutoff.value() >= self.ui.intensityThreshold.value() and self.ui.intensityCutoff > 0:
            QtGui.QMessageBox.information(self, "MetExtract",
                                          "The intensity cutoff needs to be less than the intensity threshold. Adjust the two parameters accordingly")
            return 1

        ma, ea = getIsotopeMass(self.ui.isotopeAText.text())
        mb, eb = getIsotopeMass(self.ui.isotopeBText.text())

        if not (self.labellingExperiment==TRACER) and (len(ea) < 0 or ea != eb) and not (
                        ea == "Hydrogen" and eb == "Deuterium"):
            QtGui.QMessageBox.question(self, "MetExtract",
                                       "You cannot use two different elements for the labelling process. Please enter isotopes of the same element.\nCalculation aborted",
                                       QtGui.QMessageBox.Ok)
            return 1

        if not dontSave:
            z = QtGui.QMessageBox.question(self, "MetExtract",
                                           "Do you want to save the loaded files and configuration?",
                                           QtGui.QMessageBox.Yes | QtGui.QMessageBox.No | QtGui.QMessageBox.Abort)
            if z == QtGui.QMessageBox.Yes:
                self.saveGroupsClicked()
            elif z == QtGui.QMessageBox.Abort:
                return 1

        if askStarting:
            z = QtGui.QMessageBox.question(self, "MetExtract",
                                           "Do you want to start the processing?",
                                           QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
            if z == QtGui.QMessageBox.No:
                return 1

        definedGroups=[t.data(QListWidgetItem.UserType).toPyObject() for t in natSort(self.ui.groupsList.findItems('*', QtCore.Qt.MatchWildcard), key=lambda x: str(x.data(QListWidgetItem.UserType).toPyObject().name))]
        for group in definedGroups:
            for i in range(len(group.files)):
                group.files[i]=str(group.files[i]).replace("\\", "/")

        # get all individual LC-HRMS files for processing
        files = []

        indGroups = {}
        for group in definedGroups:
            grName=str(group.name)
            indGroups[grName] = []
            for file in natSort(group.files):
                indGroups[grName].append(str(file))
                if not (file in files):
                    files.append(file)

        errorCount = 0

        cpus = min(cpu_count(), self.ui.cpuCores.value())

        if self.terminateJobs:
            return

        overallStart = time.time()

        writeMZXMLOptions = 0
        if self.ui.saveMZXML.isChecked():
            if self.ui.wm_ia.checkState() == QtCore.Qt.Checked:
                writeMZXMLOptions |= 1
            if self.ui.wm_iap.checkState() == QtCore.Qt.Checked:
                writeMZXMLOptions |= 2
            if self.ui.wm_imb.checkState() == QtCore.Qt.Checked:
                writeMZXMLOptions |= 4
            if self.ui.wm_ib.checkState() == QtCore.Qt.Checked:
                writeMZXMLOptions |= 8

        # process individual files
        if self.ui.processIndividualFiles.isChecked():
            logging.info("")
            logging.info("Processing individual LC-HRMS data on %d CPU core(s).."%min(len(files), cpus))



            # Initialize multiprocessing pool
            p = Pool(processes=min(len(files), cpus))  #, maxtasksperchild=1) # only in python >=2.7
            manager = Manager()
            lock = manager.Lock()
            queue = manager.Queue()

            pIds = {}
            for i in range(len(files)):
                pIds[i + 1] = files[i]

            # start the multiprocessing
            res = p.imap_unordered(runFile, [
                RunIdentification(files[i],
                                  exOperator=str(self.ui.exOperator_LineEdit.text()),
                                  exExperimentID=str(self.ui.exExperimentID_LineEdit.text()),
                                  exComments=str(self.ui.exComments_TextEdit.toPlainText()),
                                  exExperimentName=str(self.ui.exExperimentName_LineEdit.text()),
                                  writePDF=self.ui.savePDF.checkState() == QtCore.Qt.Checked,
                                  writeTSV=self.ui.saveCSV.checkState() == QtCore.Qt.Checked,
                                  writeMZXML=writeMZXMLOptions,
                                  metabolisationExperiment=self.labellingExperiment==TRACER,
                                  intensityThreshold=self.ui.intensityThreshold.value(),
                                  intensityCutoff=self.ui.intensityCutoff.value(),
                                  labellingisotopeA=str(self.ui.isotopeAText.text()),
                                  labellingisotopeB=str(self.ui.isotopeBText.text()),
                                  xOffset=self.isotopeBmass - self.isotopeAmass,
                                  useRatio=self.ui.useRatio.checkState()==QtCore.Qt.Checked,
                                  minRatio=self.ui.minRatio.value(),
                                  maxRatio=self.ui.maxRatio.value(),
                                  useCValidation=int(str(self.ui.useCValidation.checkState())),
                                  configuredTracers=self.configuredTracers, startTime=self.ui.scanStartTime.value(),
                                  stopTime=self.ui.scanEndTime.value(), maxLoading=self.ui.maxLoading.value(),
                                  xMin=self.ui.minXCount.value(), xMax=self.ui.maxXCount.value(),
                                  ppm=self.ui.ppmRangeIdentification.value(),
                                  isotopicPatternCountLeft=self.ui.isotopePatternCountA.value(),
                                  isotopicPatternCountRight=self.ui.isotopePatternCountB.value(),
                                  lowAbundanceIsotopeCutoff=self.ui.isoAbundance.checkState() == QtCore.Qt.Checked,
                                  intensityThresholdIsotopologs=self.ui.intensityThresholdIsotopologs.value(),
                                  purityN=self.ui.isotopicAbundanceA.value(),
                                  purityL=self.ui.isotopicAbundanceB.value(), intensityErrorN=self.ui.baseRange.value(),
                                  intensityErrorL=self.ui.isotopeRange.value(),
                                  minSpectraCount=self.ui.minSpectraCount.value(), clustPPM=self.ui.clustPPM.value(),
                                  chromPeakPPM=self.ui.wavelet_EICppm.value(),
                                  eicSmoothingWindow=str(self.ui.eicSmoothingWindow.currentText()),
                                  eicSmoothingWindowSize=self.ui.eicSmoothingWindowSize.value(),
                                  scales=[self.ui.wavelet_minScale.value(), self.ui.wavelet_maxScale.value()],
                                  snrTh=self.ui.wavelet_SNRThreshold.value(),
                                  peakCenterError=self.ui.peak_centerError.value(),
                                  peakScaleError=self.ui.peak_scaleError.value(),
                                  minPeakCorr=self.ui.minPeakCorr.value(),
                                  checkPeaksRatio=self.ui.checkBox_checkPeakRatio.isChecked(),
                                  minPeaksRatio=self.ui.doubleSpinBox_minPeakRatio.value(),
                                  maxPeaksRatio=self.ui.doubleSpinBox_maxPeakRatio.value(),
                                  calcIsoRatioNative=self.ui.calcIsoRatioNative_spinBox.value(),
                                  calcIsoRatioLabelled=self.ui.calcIsoRatioLabelled_spinBox.value(),
                                  calcIsoRatioMoiety=self.ui.calcIsoRatioMoiety_spinBox.value(),
                                  minCorrelationConnections=self.ui.minCorrelationConnections.value(),
                                  positiveScanEvent=str(self.ui.positiveScanEvent.currentText()),
                                  negativeScanEvent=str(self.ui.negativeScanEvent.currentText()),
                                  correctCCount=self.ui.correctcCount.checkState() == QtCore.Qt.Checked,
                                  minCorrelation=self.ui.minCorrelation.value(),
                                  hAIntensityError=self.ui.hAIntensityError.value(),
                                  hAMinScans=self.ui.hAMinScans.value(), adducts=self.adducts, elements=self.elementsForNL,
                                  heteroAtoms=self.heteroElements, lock=lock, queue=queue, pID=i + 1,
                                  rVersion=getRVersion(), meVersion="MetExtract (%s)" % MetExtractVersion) for i in
                range(len(files))])

            pw = ProgressWrapper(pwCount=min(len(files), cpus) + 1, showIndProgress=True, indGroups=indGroups, parent=self,
                                 closeCallback=CallBackMethod(_target=interruptIndividualFilesProcessing, selfObj=self, pool=p).getRunMethod())

            for w in range(1, min(len(files), cpus) + 1):
                pw.setUntils(
                    {.25: "#E3E9E2", .35: "#D3D9D2", .65: "#52626D", .70: "#6C7B8B", .75: "#92A1AB", .8: "#CACDC9",
                     1.: "darkorange"}, w)
                pw.setMax(1., w)

            pw.show()
            pwMain = pw.getCallingFunction(0)
            pwMain("max")(len(files))
            pwMain("value")(0)


            start = time.time()

            # monitor processing of individual LC-HRMS files and report to the user
            loop = True
            freeSlots = range(min(len(files), cpus))
            assignedThreads = {}
            while loop and not self.terminateJobs:
                completed = res._index
                if completed == len(files):
                    loop = False
                else:
                    pwMain("value")(completed)

                    mess = {}
                    while not (queue.empty()):
                        mes = queue.get(block=False, timeout=1)
                        if mes.pid not in mess:
                            mess[mes.pid] = {}
                        mess[mes.pid][mes.mes] = mes

                    for v in mess.values():
                        if "end" in v.keys() or "failed" in v.keys():
                            mes = None
                            if "end" in v.keys():
                                mes = v["end"]
                            elif "failed" in v.keys():
                                mes = v["failed"]
                            freeS = assignedThreads[mes.pid]
                            pw.getCallingFunction(assignedThreads[mes.pid] + 1)("text")("")
                            pw.getCallingFunction(assignedThreads[mes.pid] + 1)("value")(0)

                            pw.getCallingFunction()("statuscolor")(pIds[mes.pid],"olivedrab" if mes.mes == "end" else "firebrick")
                            pw.getCallingFunction()("statustext")(pIds[mes.pid], text="File: %s\nStatus: %s" % (pIds[mes.pid], "finished" if mes.mes == "end" else "failed"))

                            if freeS == -1:
                                logging.error("Something went wrong..")
                                logging.error("Progress bars do not work correctly, but files will be processed and \"finished..\" will be printed..")
                            else:
                                assignedThreads[mes.pid] = -1
                                freeSlots.append(freeS)

                    for v in mess.values():
                        if "start" in v.keys():
                            mes = v["start"]
                            if len(freeSlots) > 0:
                                w = freeSlots.pop()
                                assignedThreads[mes.pid] = w

                                pw.getCallingFunction()("statuscolor")(pIds[mes.pid], "orange")
                                pw.getCallingFunction()("statustext")(pIds[mes.pid],text="File: %s\nStatus: %s\nProcess ID: %d" % (pIds[mes.pid], "processing", mes.pid))
                            else:
                                logging.error("Something went wrong..")
                                logging.error("Progress bars do not work correctly, but files will be processed and \"finished..\" will be printed..")

                    for v in mess.values():
                        for mes in v.values():
                            if mes.mes in ["log", "text", "max", "value"]:
                                if assignedThreads.has_key(mes.pid):
                                    pw.getCallingFunction(assignedThreads[mes.pid] + 1)(mes.mes)(mes.val)
                                else:
                                    logging.error("Error %d" % mes.pid)

                    elapsed = (time.time() - start) / 60.
                    hours = ""
                    if elapsed >= 60.:
                        if elapsed < 120.:
                            hours = "1 hour "
                        else:
                            hours = "%d hours " % (elapsed // 60)

                    pwMain("text")("%s %.2f min elapsed %d / %d files done (%d parallel)" % (
                        hours, elapsed % 60, completed, len(files), min(cpus, len(files))))
                    time.sleep(.5)

            pw.setSkipCallBack(True)
            pw.hide()

            # Log time used for processing of individual files
            elapsed = (time.time() - start) / 60.
            hours = ""
            if elapsed >= 60.:
                if elapsed < 120.:
                    hours = "1 hour "
                else:
                    hours = "%d hours " % (elapsed // 60)
            mins = "%.2f min(s)" % (elapsed % 60.)

            if not self.terminateJobs:
                logging.info("Individual files processed (%s%s).." % (hours, mins))

            p.close()
            p.terminate()
            p.join()



        if self.terminateJobs:
            return

        # bracket/group from individual LC-HRMS data / re-integrate missed peaks
        if self.ui.processMultipleFiles.checkState() == QtCore.Qt.Checked:
            start = time.time()

            # bracket/group results from individual LC-HRMS data
            if self.ui.groupResults.isChecked():
                logging.info("Bracketing of individual LC-HRMS results..")
                start = time.time()

                pw = ProgressWrapper(1, parent=self)
                pw.show()
                pw.getCallingFunction()("text")("Bracketing results")
                pw.getCallingFunction()("header")("Bracketing..")

                try:
                    #Group results
                    generalProcessingParams=Bunch(exOperator=str(self.ui.exOperator_LineEdit.text()),
                                                  exExperimentID=str(self.ui.exExperimentID_LineEdit.text()),
                                                  exComments=str(self.ui.exComments_TextEdit.toPlainText()),
                                                  exExperimentName=str(self.ui.exExperimentName_LineEdit.text()),
                                                  writePDF=self.ui.savePDF.checkState() == QtCore.Qt.Checked,
                                                  writeTSV=self.ui.saveCSV.checkState() == QtCore.Qt.Checked,
                                                  writeMZXML=writeMZXMLOptions,
                                                  metabolisationExperiment=self.labellingExperiment==TRACER,
                                                  intensityThreshold=self.ui.intensityThreshold.value(),
                                                  intensityCutoff=self.ui.intensityCutoff.value(),
                                                  labellingisotopeA=str(self.ui.isotopeAText.text()),
                                                  labellingisotopeB=str(self.ui.isotopeBText.text()),
                                                  xOffset=self.isotopeBmass - self.isotopeAmass,
                                                  useRatio=self.ui.useRatio.checkState()==QtCore.Qt.Checked,
                                                  minRatio=self.ui.minRatio.value(),
                                                  maxRatio=self.ui.maxRatio.value(),
                                                  useCValidation=int(str(self.ui.useCValidation.checkState())),
                                                  configuredTracers="[%s]"%",".join([str(t) for t in self.configuredTracers]), startTime=self.ui.scanStartTime.value(),
                                                  stopTime=self.ui.scanEndTime.value(), maxLoading=self.ui.maxLoading.value(),
                                                  xMin=self.ui.minXCount.value(), xMax=self.ui.maxXCount.value(),
                                                  ppm=self.ui.ppmRangeIdentification.value(),
                                                  isotopicPatternCountLeft=self.ui.isotopePatternCountA.value(),
                                                  isotopicPatternCountRight=self.ui.isotopePatternCountB.value(),
                                                  lowAbundanceIsotopeCutoff=self.ui.isoAbundance.checkState() == QtCore.Qt.Checked,
                                                  intensityThresholdIsotopologs=self.ui.intensityThresholdIsotopologs.value(),
                                                  purityN=self.ui.isotopicAbundanceA.value(),
                                                  purityL=self.ui.isotopicAbundanceB.value(), intensityErrorN=self.ui.baseRange.value(),
                                                  intensityErrorL=self.ui.isotopeRange.value(),
                                                  minSpectraCount=self.ui.minSpectraCount.value(), clustPPM=self.ui.clustPPM.value(),
                                                  chromPeakPPM=self.ui.wavelet_EICppm.value(),
                                                  eicSmoothingWindow=str(self.ui.eicSmoothingWindow.currentText()),
                                                  eicSmoothingWindowSize=self.ui.eicSmoothingWindowSize.value(),
                                                  scales=[self.ui.wavelet_minScale.value(), self.ui.wavelet_maxScale.value()],
                                                  snrTh=self.ui.wavelet_SNRThreshold.value(),
                                                  peakCenterError=self.ui.peak_centerError.value(),
                                                  peakScaleError=self.ui.peak_scaleError.value(),
                                                  minPeakCorr=self.ui.minPeakCorr.value(),
                                                  checkPeaksRatio=self.ui.checkBox_checkPeakRatio.isChecked(),
                                                  minPeaksRatio=self.ui.doubleSpinBox_minPeakRatio.value(),
                                                  maxPeaksRatio=self.ui.doubleSpinBox_maxPeakRatio.value(),
                                                  calcIsoRatioNative=self.ui.calcIsoRatioNative_spinBox.value(),
                                                  calcIsoRatioLabelled=self.ui.calcIsoRatioLabelled_spinBox.value(),
                                                  calcIsoRatioMoiety=self.ui.calcIsoRatioMoiety_spinBox.value(),
                                                  minCorrelationConnections=self.ui.minCorrelationConnections.value(),
                                                  positiveScanEvent=str(self.ui.positiveScanEvent.currentText()),
                                                  negativeScanEvent=str(self.ui.negativeScanEvent.currentText()),
                                                  correctCCount=self.ui.correctcCount.checkState() == QtCore.Qt.Checked,
                                                  minCorrelation=self.ui.minCorrelation.value(),
                                                  hAIntensityError=self.ui.hAIntensityError.value(),
                                                  hAMinScans=self.ui.hAMinScans.value(),
                                                  adducts="[%s]"%",".join([str(a) for a in self.adducts]),
                                                  elements="[%s]"%",".join([str(e) for e in self.elementsForNL]),
                                                  heteroAtoms="[%s]"%",".join([str(h) for h in self.heteroElements]),
                                                  rVersion=getRVersion(), meVersion="MetExtract (%s)" % MetExtractVersion)
                    procProc = FuncProcess(_target=bracketResults,
                                            indGroups=indGroups, minX=self.ui.minXCount.value(), maxX=self.ui.maxXCount.value(),
                                            groupSizePPM=self.ui.groupPpm.value(),
                                            maxTimeDeviation=self.ui.groupingRT.value() * 60.,
                                            maxLoading=self.ui.maxLoading.value(),
                                            positiveScanEvent=str(self.ui.positiveScanEvent.currentText()),
                                            negativeScanEvent=str(self.ui.negativeScanEvent.currentText()),
                                            file=str(self.ui.groupsSave.text()),
                                            align=(self.ui.alignChromatograms.checkState() == QtCore.Qt.Checked),
                                            nPolynom=self.ui.polynomValue.value(),
                                            rVersion=getRVersion(), meVersion="MetExtract (%s)" % MetExtractVersion,
                                            generalProcessingParams=generalProcessingParams)
                    procProc.addKwd("pwMaxSet", procProc.getQueue())
                    procProc.addKwd("pwValSet", procProc.getQueue())
                    procProc.start()

                    pw.setCloseCallback(closeCallBack=CallBackMethod(_target=interruptBracketingOfFeaturePairs, selfObj=self, funcProc=procProc).getRunMethod())

                    # check for status updates
                    while procProc.isAlive():
                        QtGui.QApplication.processEvents();

                        while not (procProc.getQueue().empty()):
                            mes = procProc.getQueue().get(block=False, timeout=1)

                            # No idea why / where there are sometimes other objects than Bunch(mes, val), but they occur
                            if isinstance(mes, Bunch) and hasattr(mes, "mes") and hasattr(mes, "val"):
                                pw.getCallingFunction()(mes.mes)(mes.val)
                            else:
                                logging.critical("UNKNONW OBJECT IN PROCESSING QUEUE:", mes)

                        time.sleep(.5)

                    # Log time used for bracketing
                    elapsed = (time.time() - start) / 60.
                    hours = ""
                    if elapsed >= 60.:
                        if elapsed < 120.:
                            hours = "1 hour "
                        else:
                            hours = "%d hours " % (elapsed // 60)
                    mins = "%.2f min(s)" % (elapsed % 60.)

                    if self.terminateJobs:
                        return
                    else:
                        logging.info("Bracketing finished (%s%s).." % (hours, mins))

                    #Arrange grouped results and add statistics columns
                    groups = {}
                    outputOrder = []

                    pw.getCallingFunction()("text")("Adding statistics columns")
                    if False:
                        for group in definedGroups:
                            preFix="_Stat_N"
                            grpName=group.name+preFix
                            grpAdd(groups, group.name+preFix, group.minFound,
                                   [grp[(grp.rfind("/") + 1):max(grp.lower().rfind(".mzxml"), grp.lower().rfind(".mzml"))] + "_Area_N" for grp in
                                    natSort(group.files)])
                            outputOrder.append(grpName)
                        for group in definedGroups:
                            preFix="_Stat_L"
                            grpName=group.name+preFix
                            grpAdd(groups, group.name+preFix, group.minFound,
                                   [grp[(grp.rfind("/") + 1):max(grp.lower().rfind(".mzxml"), grp.lower().rfind(".mzml"))] + "_Area_L" for grp in
                                    natSort(group.files)])
                            outputOrder.append(grpName)
                    grpStats=[]
                    for group in definedGroups:
                        preFix="_Stat_fold"
                        grpName=group.name+preFix
                        grpAdd(groups, group.name+preFix, group.minFound,
                               [grp[(grp.rfind("/") + 1):max(grp.lower().rfind(".mzxml"), grp.lower().rfind(".mzml"))] + "_fold" for grp in
                                natSort(group.files)])
                        outputOrder.append(grpName)
                        grpStats.append((str(group.name+"_Stat_fold"), group.minFound, group.omitFeatures))

                    addStatsColumnToResults(str(self.ui.groupsSave.text()), groups, str(self.ui.groupsSave.text()), outputOrder)

                    #remove feature pairs not found more than n times (according to user specified omit value)
                    grpOmit(str(self.ui.groupsSave.text()), grpStats, str(self.ui.groupsSave.text()))
                    logging.info("Statistic columns added (and feature pairs omitted)..")


                except Exception as ex:
                    import traceback
                    traceback.print_exc()
                    logging.error(str(traceback))

                    QtGui.QMessageBox.warning(self, "MetExtract", "Error during bracketing of files: '%s'" % str(ex), QtGui.QMessageBox.Ok)
                    errorCount += 1


                if self.terminateJobs:
                    return

                try:
                    # Calculate metabolite groups
                    pw.getCallingFunction()("text")("Convoluting feature pairs")

                    runIdentificationInstance=RunIdentification(files[0],
                                  exOperator=str(self.ui.exOperator_LineEdit.text()),
                                  exExperimentID=str(self.ui.exExperimentID_LineEdit.text()),
                                  exComments=str(self.ui.exComments_TextEdit.toPlainText()),
                                  exExperimentName=str(self.ui.exExperimentName_LineEdit.text()),
                                  writePDF=False,
                                  writeTSV=False,
                                  writeMZXML=0,
                                  metabolisationExperiment=self.labellingExperiment==TRACER,
                                  intensityThreshold=self.ui.intensityThreshold.value(),
                                  intensityCutoff=self.ui.intensityCutoff.value(),
                                  labellingisotopeA=str(self.ui.isotopeAText.text()),
                                  labellingisotopeB=str(self.ui.isotopeBText.text()),
                                  xOffset=self.isotopeBmass - self.isotopeAmass,
                                  useRatio=self.ui.useRatio.checkState()==QtCore.Qt.Checked,
                                  minRatio=self.ui.minRatio.value(),
                                  maxRatio=self.ui.maxRatio.value(),
                                  useCValidation=int(str(self.ui.useCValidation.checkState())),
                                  configuredTracers=self.configuredTracers, startTime=self.ui.scanStartTime.value(),
                                  stopTime=self.ui.scanEndTime.value(), maxLoading=self.ui.maxLoading.value(),
                                  xMin=self.ui.minXCount.value(), xMax=self.ui.maxXCount.value(),
                                  ppm=self.ui.ppmRangeIdentification.value(),
                                  isotopicPatternCountLeft=self.ui.isotopePatternCountA.value(),
                                  isotopicPatternCountRight=self.ui.isotopePatternCountB.value(),
                                  lowAbundanceIsotopeCutoff=self.ui.isoAbundance.checkState() == QtCore.Qt.Checked,
                                  purityN=self.ui.isotopicAbundanceA.value(),
                                  purityL=self.ui.isotopicAbundanceB.value(), intensityErrorN=self.ui.baseRange.value(),
                                  intensityErrorL=self.ui.isotopeRange.value(),
                                  minSpectraCount=self.ui.minSpectraCount.value(), clustPPM=self.ui.clustPPM.value(),
                                  chromPeakPPM=self.ui.wavelet_EICppm.value(),
                                  eicSmoothingWindow=str(self.ui.eicSmoothingWindow.currentText()),
                                  eicSmoothingWindowSize=self.ui.eicSmoothingWindowSize.value(),
                                  scales=[self.ui.wavelet_minScale.value(), self.ui.wavelet_maxScale.value()],
                                  snrTh=self.ui.wavelet_SNRThreshold.value(),
                                  peakCenterError=self.ui.peak_centerError.value(),
                                  peakScaleError=self.ui.peak_scaleError.value(),
                                  minPeakCorr=self.ui.minPeakCorr.value(),
                                  minCorrelationConnections=self.ui.minCorrelationConnections.value(),
                                  positiveScanEvent=str(self.ui.positiveScanEvent.currentText()),
                                  negativeScanEvent=str(self.ui.negativeScanEvent.currentText()),
                                  correctCCount=self.ui.correctcCount.checkState() == QtCore.Qt.Checked,
                                  minCorrelation=self.ui.minCorrelation.value(),
                                  hAIntensityError=self.ui.hAIntensityError.value() / 100.,
                                  hAMinScans=self.ui.hAMinScans.value(), adducts=self.adducts, elements=self.elementsForNL,
                                  heteroAtoms=self.heteroElements, lock=None, queue=None, pID=1,
                                  rVersion=getRVersion(), meVersion="MetExtract (%s)" % MetExtractVersion)

                    calculateMetaboliteGroups(str(self.ui.groupsSave.text()), definedGroups,
                                              minConnectionsInFiles=self.ui.metaboliteClusterMinConnections.value(), minConnectionRate=self.ui.minConnectionRate.value(),
                                              runIdentificationInstance=runIdentificationInstance)
                    logging.info("Metabolite groups calculated..")

                except Exception as ex:
                    import traceback
                    traceback.print_exc()
                    logging.error(str(traceback))

                    QtGui.QMessageBox.warning(self, "MetExtract", "Error during convolution of feature pairs: '%s'" % str(ex), QtGui.QMessageBox.Ok)
                    errorCount += 1
                finally:
                    pw.setSkipCallBack(True)
                    pw.hide()

            if self.terminateJobs:
                return

            # re-integrate missed peaks
            if self.ui.integratedMissedPeaks.isChecked():
                logging.info("Re-integrating of individual LC-HRMS results..")
                start = time.time()

                pw = ProgressWrapper(min(len(files), cpus) + 1, showLog=False, parent=self)
                pw.show()
                pw.getCallingFunction()("text")("Integrating..")
                pw.getCallingFunction()("header")("Integrating..")

                try:
                    #Reintegrate missed peaks in files
                    fDict = {}
                    for group in definedGroups:
                        for grp in natSort(group.files):
                            f = grp
                            f = f.replace("\\", "/")
                            fDict[f] = f[(f.rfind("/") + 1):max(f.lower().rfind(".mzxml"), f.lower().rfind(".mzml"))]

                    integrateResultsFile(str(self.ui.groupsSave.text()), str(self.ui.groupsSave.text()), fDict, "MZ",
                                         "RT", "Xn", "Charge", "L_MZ", "Ionisation_Mode", "Num",
                                         ppm=self.ui.groupPpm.value(),
                                         maxRTShift=self.ui.integrationMaxTimeDifference.value(),
                                         scales=[self.ui.wavelet_minScale.value(), self.ui.wavelet_maxScale.value()],
                                         reintegrateIntensityCutoff=self.ui.reintegrateIntensityCutoff.value(),
                                         snrTH=self.ui.wavelet_SNRThreshold.value(),
                                         smoothingWindow=str(self.ui.eicSmoothingWindow.currentText()),
                                         smoothingWindowSize=self.ui.eicSmoothingWindowSize.value(),
                                         positiveScanEvent=str(self.ui.positiveScanEvent.currentText()),
                                         negativeScanEvent=str(self.ui.negativeScanEvent.currentText()),
                                         pw=pw, selfObj=self, cpus=min(len(files), cpus))
                    # Log time used for bracketing
                    elapsed = (time.time() - start) / 60.
                    hours = ""
                    if elapsed >= 60.:
                        if elapsed < 120.:
                            hours = "1 hour "
                        else:
                            hours = "%d hours " % (elapsed // 60)
                    mins = "%.2f min(s)" % (elapsed % 60.)
                    logging.info("Re-integrating finished (%s%s).." % (hours, mins))

                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    logging.error(str(traceback))

                    QtGui.QMessageBox.warning(self, "MetExtract", "Error during reintegrating files: '%s'" % str(e),
                                              QtGui.QMessageBox.Ok)
                    errorCount += 1
                finally:
                    pw.setSkipCallBack(True)
                    pw.hide()

            if self.terminateJobs:
                return




        self.updateLCMSSampleSettings(force=True)

        # Log time used for bracketing
        elapsed = (time.time() - overallStart) / 60.
        hours = ""
        if elapsed >= 60.:
            if elapsed < 120.:
                hours = "1 hour "
            else:
                hours = "%d hours " % (elapsed // 60)
        mins = "%.2f min(s)" % (elapsed % 60.)

        if errorCount == 0:
            logging.info("Processing successfully finished (%s%s)..\n"%(hours, mins))
        else:
            logging.warning("Processing finished with %d errors (%s%s)..\n" % (errorCount, hours, mins))


    def groupFilesChanges(self, sta):
        self.ui.label_26.setEnabled(sta)
        self.ui.groupPpm.setEnabled(sta)
        self.ui.alignChromatograms.setEnabled(sta)
        self.ui.label_18.setEnabled(sta)
        self.ui.groupingRT.setEnabled(sta)

    def processMultipleFilesChanged(self, sta):
        self.ui.groupResults.setEnabled(sta)
        self.ui.label_26.setEnabled(sta)
        self.ui.groupPpm.setEnabled(sta)
        self.ui.alignChromatograms.setEnabled(sta)
        self.ui.label_18.setEnabled(sta)
        self.ui.groupingRT.setEnabled(sta)
        self.ui.integratedMissedPeaks.setEnabled(sta)
        self.ui.label_20.setEnabled(sta)
        self.ui.groupsSave.setEnabled(sta)
        self.ui.groupsSelectFile.setEnabled(sta)

    def saveMZXMLChanged(self, sta):
        self.ui.wm_ia.setEnabled(sta)
        self.ui.wm_iap.setEnabled(sta)
        self.ui.wm_imb.setEnabled(sta)
        self.ui.wm_ib.setEnabled(sta)

    def procIndFilesChanges(self, sta):
        self.ui.procIndFiles.setEnabled(sta)

    #<editor-fold desc="### visualisation of results of single sample">
    def closeCurrentOpenResultsFile(self):
        if hasattr(self, "currentOpenResultsFile") and self.currentOpenResultsFile is not None:
            self.currentOpenResultsFile.curs.close()
            self.currentOpenResultsFile.conn.close()
            self.currentOpenResultsFile.file = None
            self.currentOpenResultsFile=None

    def openFileAsCurrentOpenResultsFile(self, file):
        if os.path.exists(file + ".identified.sqlite") and os.path.isfile(file + ".identified.sqlite"):
            self.currentOpenResultsFile = Bunch()
            self.currentOpenResultsFile.file = file
            self.currentOpenResultsFile.conn = connect(file + ".identified.sqlite")
            self.currentOpenResultsFile.curs = self.currentOpenResultsFile.conn.cursor()
            return True
        else:
            return False

    def selectedResultChanged(self, ind):

        self.ui.res_ExtractedData.clear()
        self.ui.chromPeakName.setText("")

        cInd = self.ui.processedFilesComboBox.currentIndex()
        b = self.ui.processedFilesComboBox.itemData(cInd).toPyObject()

        if not hasattr(b, "file") or b.file is None:
            return -1

        self.closeCurrentOpenResultsFile()
        if self.openFileAsCurrentOpenResultsFile(b.file):

            it = QtGui.QTreeWidgetItem(["MZs"])
            self.ui.res_ExtractedData.addTopLevelItem(it)
            it.myType = "MZs"
            count = 0
            children=[]
            pw=ProgressWrapper(pwCount=4)


            ## Load mz pairs
            pw.setText("Fetching mzs", i=0)
            pw.setText("", i=1)
            pw.setText("", i=2)
            pw.setTextu("", i=3)
            pw.show()

            numberOfMZs=0
            for row in SQLSelectAsObject(self.currentOpenResultsFile.curs, "SELECT count(id) AS co FROM MZs"):
                numberOfMZs=row.co
            pw.setMax(numberOfMZs)

            maxMZsFetch=50000

            if numberOfMZs<maxMZsFetch:

                pw.setText("Fetching mzs (%d)"%numberOfMZs, i=0)

                for mzRes in SQLSelectAsObject(self.currentOpenResultsFile.curs, "SELECT id, mz, xcount, scanid, loading, scantime, intensity FROM MZs ORDER BY scanid"):
                    d = QtGui.QTreeWidgetItem(it, [str(s) for s in [mzRes.mz, mzRes.xcount, mzRes.scanid, "%.2f min / %.2f sec"%(mzRes.scantime/60., mzRes.scantime), mzRes.loading, "%.1f"%mzRes.intensity]])
                    d.myType = "mz"
                    d.myData=mzRes
                    d.myID=int(mzRes.id)
                    children.append(d)
                    count += 1

                    pw.setValueu(count, i=0)
                pw.setText("%d MZs fetched"%numberOfMZs, i=0)
            else:
                pw.setTextu("Mzs not fetched (too many; %d)"%numberOfMZs)
            it.addChildren(children)
            it.setText(1, "%d"%numberOfMZs)


            it = QtGui.QTreeWidgetItem(["MZ bins"])
            it.myType = "MZBins"
            self.ui.res_ExtractedData.addTopLevelItem(it)
            mzbins = []

            for mzbin in SQLSelectAsObject(self.currentOpenResultsFile.curs, "SELECT id, mz FROM MZBins ORDER BY mz"):
                mzbins.append(mzbin)
            count = 0
            children=[]


            ## Load mz bins
            pw.setText("Fetching MZBins (%d)"%len(mzbins), i=1)
            pw.setMax(len(mzbins), i=1)
            if len(mzbins)<2500:
                for mzbin in mzbins:
                    d = QtGui.QTreeWidgetItem([str(mzbin.mz)])
                    d.myType = "mzbin"
                    d.myID = int(mzbin.id)
                    children.append(d)
                    countinner = 0
                    minInner = 1000000.
                    maxInner = 0.
                    xcount = 0



                    for mzRes in SQLSelectAsObject(self.currentOpenResultsFile.curs, "SELECT id, mz, "
                                                                                     "xcount, "
                                                                                     "scanid, "
                                                                                     "loading, "
                                                                                     "scantime, "
                                                                                     "intensity "
                                                                                     "FROM MZs m, MZBinsKids k "
                                                                                     "WHERE m.id==k.mzID AND k.mzbinID=%d ORDER BY m.scanid" % mzbin.id):
                        minInner = min(float(mzRes.mz), minInner)
                        maxInner = max(maxInner, mzRes.mz)
                        xcount = int(mzRes.xcount)
                        if numberOfMZs<maxMZsFetch:
                            dd = QtGui.QTreeWidgetItem([str(s) for s in [mzRes.mz, mzRes.xcount, mzRes.scanid, "%.2f min / %.2f sec"%(mzRes.scantime/60., mzRes.scantime), mzRes.loading, "%.1f"%mzRes.intensity]])
                            dd.myType = "mz"
                            dd.myData=mzRes
                            dd.myID=int(mzRes.id)
                            d.addChild(dd)
                        countinner += 1

                    d.setText(0, "%.5f (%d)" % (mzbin.mz, countinner))
                    d.setText(1, "%.4f" % ((maxInner - minInner) * 1000000. / minInner))
                    d.setText(2, "%d" % xcount)
                    count += 1
                    pw.setValueu(count, i=1)
                pw.setText("%d MZBins fetched"%count, i=1)
            else:
                pw.setTextu("MZBins not fetched (too many; %d)"%len(mzbins), i=1)


            it.addChildren(children)
            it.setText(1, "%d" % len(mzbins))



            ## Load feature pairs
            it = QtGui.QTreeWidgetItem(["Feature pairs"])
            self.ui.res_ExtractedData.addTopLevelItem(it)
            it.myType = "Features"

            numberOfFeaturePairs=0
            for row in SQLSelectAsObject(self.currentOpenResultsFile.curs, "SELECT count(chromPeaks.id) AS co "
                                                                           "FROM chromPeaks LEFT JOIN tracerConfiguration ON tracerConfiguration.id=chromPeaks.tracer"):
                numberOfFeaturePairs=row.co

            pw.setTextu("Fetching feature pairs (%d)"%numberOfFeaturePairs, i=2)
            pw.setMaxu(numberOfFeaturePairs, i=2)

            count = 0
            children=[]
            for row in SQLSelectAsObject(self.currentOpenResultsFile.curs, "SELECT chromPeaks.id AS cpID, "
                                                                           "mz, "
                                                                           "xcount, "
                                                                           "NPeakCenterMin, "
                                                                           "LPeakCenterMin, "
                                                                           "NPeakCenter, "
                                                                           "NPeakScale, "
                                                                           "eicID, "
                                                                           "LPeakScale, "
                                                                           "peaksCorr, "
                                                                           "peaksRatio, "
                                                                           "LPeakCenter, "
                                                                           "LPeakScale, "
                                                                           "Loading, "
                                                                           "heteroAtoms, "
                                                                           "name AS tracerName, "
                                                                           "adducts, "
                                                                           "heteroAtoms, "
                                                                           "chromPeaks.ionMode AS ionMode, "
                                                                           "assignedName, "
                                                                           "NBorderLeft, "
                                                                           "NBorderRight, "
                                                                           "LBorderLeft, "
                                                                           "LBorderRight, "
                                                                           "NPeakArea, "
                                                                           "LPeakArea, "
                                                                           "chromPeaks.massSpectrumID AS massSpectrumID, "
                                                                           "assignedMZs "
                                                                           "FROM chromPeaks LEFT JOIN tracerConfiguration ON tracerConfiguration.id=chromPeaks.tracer "
                                                                           "ORDER BY tracerConfiguration.id, NPeakCenter, mz, xcount"):
                adducts = ""
                lk = loads(base64.b64decode(row.adducts))
                if len(lk) > 0:
                    adducts = ", ".join(lk)

                heteroAtoms = []
                lk = loads(base64.b64decode(row.heteroAtoms))
                for hetAtom in lk:
                    pIso = lk[hetAtom]
                    for hetAtomCount in pIso:
                        heteroAtoms.append("%s%d"%(hetAtom, hetAtomCount))
                heteroAtoms=", ".join(heteroAtoms)

                d = QtGui.QTreeWidgetItem([str(row.mz) + " (/" + str(row.ionMode) + str(row.Loading) + ") ",
                                           "%.2f / %.2f" % (float(row.NPeakCenterMin) / 60., float(row.LPeakCenterMin) / 60.),
                                           str(row.xcount),
                                           "%s / %s "%(adducts, heteroAtoms),
                                           "%.0f / %.0f" % (float(row.NPeakScale), float(row.LPeakScale)),
                                           "%.3f"%(row.peaksCorr),
                                           "%.3f / %.3f"%(row.peaksRatio, row.NPeakArea/row.LPeakArea),
                                           "%.1f / %.1f"%(row.NPeakArea, row.LPeakArea),
                                           "%d"%row.assignedMZs,
                                           str(row.tracerName)])

                xp = ChromPeakPair(NPeakCenter=int(row.NPeakCenter), LPeakScale=float(row.LPeakScale), LPeakCenter=int(row.LPeakCenter),
                               NPeakScale=float(row.NPeakScale), NSNR=0, NPeakArea=-1, mz=float(row.mz), xCount=int(row.xcount),
                               NBorderLeft=float(row.NBorderLeft), NBorderRight=float(row.NBorderRight),
                               LBorderLeft=float(row.LBorderLeft), LBorderRight=float(row.LBorderRight),
                               NPeakCenterMin=float(row.NPeakCenterMin), LPeakCenterMin=float(row.LPeakCenterMin), eicID=int(row.eicID), massSpectrumID=int(row.massSpectrumID),
                               assignedName=str(row.assignedName), id=int(row.cpID), loading=int(row.Loading), peaksCorr=float(row.peaksCorr), peaksRatio=float(row.peaksRatio),
                               tracer=str(row.tracerName), ionMode=str(row.ionMode), heteroAtoms=heteroAtoms, adducts=adducts)

                d.myType = "feature"
                d.myID = int(row.cpID)
                d.myData = xp
                children.append(d)
                count += 1

                pw.setValueu(count, i=2)
            pw.setTextu("%d feature pairs fetched", i=2)

            it.addChildren(children)
            it.setExpanded(False)
            it.setText(1, "%d" % count)

            it = QtGui.QTreeWidgetItem(["Feature groups"])
            self.ui.res_ExtractedData.addTopLevelItem(it)
            it.myType = "Feature Groups"
            it.setExpanded(True)
            count = 0
            fGs = []

            children=[]
            for row in SQLSelectAsObject(self.currentOpenResultsFile.curs, "SELECT featureGroups.id AS fgID, "
                                                                           "featureName AS featureName, "
                                                                           "name AS tracerName "
                                                                           "FROM featureGroups INNER JOIN tracerConfiguration ON featureGroups.tracer=tracerConfiguration.id "
                                                                           "ORDER BY tracerConfiguration.id, featureGroups.id"):
                fGs.append(row)


            ## Load Feature groups
            pw.setText("Fetching feature groups (%d)"%len(fGs), i=3)
            pw.setMax(len(fGs), i=3)
            for fG in fGs:

                d = QtGui.QTreeWidgetItem([str(fG.featureName), "", str(fG.fgID), "", "", "", "", "", "", str(fG.tracerName), ""])
                d.myType = "featureGroup"
                d.myID = fG.fgID
                d.myData = fG
                d.setExpanded(True)
                cpCount = 0
                sumRt = 0.
                for row in SQLSelectAsObject(self.currentOpenResultsFile.curs, "SELECT c.id AS cpID, "
                                                                               "c.mz AS mz, "
                                                                               "c.xcount AS xcount, "
                                                                               "c.NPeakCenterMin AS NPeakCenterMin, "
                                                                               "c.LPeakCenterMin AS LPeakCenterMin, "
                                                                               "c.NPeakCenter AS NPeakCenter, "
                                                                               "c.NPeakScale AS NPeakScale, "
                                                                               "c.NBorderLeft AS NBorderLeft, "
                                                                               "c.NBorderRight AS NBorderRight, "
                                                                               "c.LBorderLeft AS LBorderLeft, "
                                                                               "c.LBorderRight AS LBorderRight, "
                                                                               "c.eicID AS eicID, "
                                                                               "f.fDesc AS fDesc, "
                                                                               "c.LPeakScale AS LPeakScale, "
                                                                               "c.LPeakCenter AS LPeakCenter, "
                                                                               "c.peaksCorr AS peaksCorr, "
                                                                               "c.peaksRatio AS peaksRatio, "
                                                                               "c.Loading AS Loading, "
                                                                               "c.NPeakArea AS NPeakArea, "
                                                                               "c.LPeakArea AS LPeakArea, "
                                                                               "t.name AS tracerName, "
                                                                               "adducts AS adducts, "
                                                                               "heteroAtoms AS heteroAtoms, "
                                                                               "c.assignedName AS assignedName, "
                                                                               "c.ionMode AS ionMode, "
                                                                               "c.massSpectrumID AS massSpectrumID, "
                                                                               "c.assignedMZs AS assignedMZs "
                                                                               "FROM chromPeaks c JOIN featureGroupFeatures f ON c.id==f.fID INNER JOIN tracerConfiguration t ON t.id=c.tracer "
                                                                               "WHERE f.fGroupID=%d ORDER BY c.mz, c.xcount" %
                                fG.fgID):

                    adducts = ""
                    lk = loads(base64.b64decode(row.adducts))
                    if len(lk) > 0:
                        adducts = ", ".join(lk)

                    heteroAtoms = []
                    lk = loads(base64.b64decode(row.heteroAtoms))
                    for hetAtom in lk:
                        pIso = lk[hetAtom]
                        for hetAtomCount in pIso:
                            heteroAtoms.append("%s%d"%(hetAtom, hetAtomCount))
                    heteroAtoms=", ".join(heteroAtoms)

                    xp = ChromPeakPair(NPeakCenter=int(row.NPeakCenter), loading=int(row.Loading), LPeakScale=float(row.LPeakScale),
                                   LPeakCenter=int(row.LPeakCenter), NPeakScale=float(row.NPeakScale), NSNR=0, NPeakArea=-1,
                                   mz=float(row.mz), xCount=int(row.xcount), NPeakCenterMin=float(row.NPeakCenterMin),
                                   NBorderLeft=float(row.NBorderLeft), NBorderRight=float(row.NBorderRight),
                                   LBorderLeft=float(row.LBorderLeft), LBorderRight=float(row.LBorderRight),
                                   LPeakCenterMin=float(row.LPeakCenterMin), eicID=int(row.eicID), massSpectrumID=int(row.massSpectrumID),
                                   assignedName=str(row.assignedName), id=int(row.cpID),
                                   tracer=str(row.tracerName), ionMode=str(row.ionMode), adducts=adducts, heteroAtoms=heteroAtoms)
                    xp.fDesc = str(row.fDesc)
                    xp.peaksCorr = float(row.peaksCorr)
                    xp.peaksRatio = float(row.peaksRatio)

                    g = QtGui.QTreeWidgetItem([str(row.mz) + " (/" + str(row.ionMode) + str(row.Loading) + ") ",
                                               "%.2f / %.2f" % (float(row.NPeakCenterMin) / 60., float(row.LPeakCenterMin) / 60.),
                                               str(row.xcount),
                                               "%s / %s "%(adducts, heteroAtoms),
                                               "%.0f / %.0f" % (float(row.NPeakScale), float(row.LPeakScale)),
                                               "%.3f"%(row.peaksCorr),
                                               "%.3f / %.3f"%(row.peaksRatio, row.NPeakArea/row.LPeakArea),
                                               "%.1f / %.1f"%(row.NPeakArea, row.LPeakArea),
                                               "%d"%row.assignedMZs,
                                               str(row.tracerName)])

                    g.myType = "feature"
                    g.myID = int(row.cpID)
                    g.myData = xp
                    d.addChild(g)

                    sumRt = sumRt + xp.NPeakCenterMin
                    cpCount += 1
                d.setText(2, str(cpCount))
                d.setText(1, "%.2f" % (sumRt / cpCount / 60.))
                children.append(d)
                count += 1
                pw.setValueu(count, i=3)
            it.addChildren(children)
            it.setText(1, "%d" % count)

            pw.hide()

            ## Load parameters
            it = QtGui.QTreeWidgetItem(["Parameters"]);
            it.myType = "parameter"
            self.ui.res_ExtractedData.addTopLevelItem(it)
            config = {}
            for row in self.currentOpenResultsFile.curs.execute("SELECT key, value FROM config"):
                config[str(row[0])] = str(row[1])

            if config["metabolisationExperiment"] == "True":
                itl = QtGui.QTreeWidgetItem(["Tracers"]);
                it.addChild(itl);
                itl.myType = "parameter";
                itl.myType = "parameter"
                for tracer in loads(base64.b64decode(config["configuredTracers"])):
                    itle = QtGui.QTreeWidgetItem([tracer.name]);
                    itl.addChild(itle);
                    itle.myType = "parameter"
                    itlee = QtGui.QTreeWidgetItem(["Elements", "%d" % tracer.elementCount]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
                    itlee = QtGui.QTreeWidgetItem(["Labelling", "%s/%s" % (tracer.isotopeA, tracer.isotopeB)]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
                    itlee = QtGui.QTreeWidgetItem(["Delta mz", "%.5f" % (
                        getIsotopeMass(tracer.isotopeB)[0] - getIsotopeMass(tracer.isotopeA)[0])]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
                    itlee = QtGui.QTreeWidgetItem(
                        ["%s purity" % tracer.isotopeA, "%.2f%%" % (float(tracer.enrichmentA) * 100.)]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
                    itlee = QtGui.QTreeWidgetItem(
                        ["%s purity" % tracer.isotopeB, "%.2f%%" % (float(tracer.enrichmentB) * 100.)]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
                    itlee = QtGui.QTreeWidgetItem(["%s amount" % tracer.isotopeA, "%.2f" % tracer.amountA]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
                    itlee = QtGui.QTreeWidgetItem(["%s amount" % tracer.isotopeB, "%.2f" % tracer.amountB]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
                    itlee = QtGui.QTreeWidgetItem(
                        ["%s/%s ratio" % (tracer.isotopeA, tracer.isotopeB), "%.4f" % tracer.monoisotopicRatio]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
                    itlee = QtGui.QTreeWidgetItem(["Upper error", "%.2f%%" % (float(tracer.maxRelNegBias) * 100.)]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
                    itlee = QtGui.QTreeWidgetItem(["Lower error", "%.2f%%" % (float(tracer.maxRelPosBias) * 100.)]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
            else:
                itl = QtGui.QTreeWidgetItem(["Full metabolome labelling experiment"]);
                it.addChild(itl);
                itl.myType = "parameter"
                itle = QtGui.QTreeWidgetItem(["%d%s/%d%s" % (
                    float(config["isotopeA"]), config["labellingElement"], float(config["isotopeB"]),
                    config["labellingElement"])]);
                itl.addChild(itle);
                itle.myType = "parameter"
                itle = QtGui.QTreeWidgetItem(["%d%s purity" % (float(config["isotopeA"]), config["labellingElement"]),
                                              "%.2f%%" % (float(config["purityN"]) * 100.)]);
                itl.addChild(itle);
                itle.myType = "parameter"
                itle = QtGui.QTreeWidgetItem(["%d%s purity" % (float(config["isotopeB"]), config["labellingElement"]),
                                              "%.2f%%" % (float(config["purityL"]) * 100.)]);
                itl.addChild(itle);
                itle.myType = "parameter"

                if "useRatio" in config.keys() and config["useRatio"]:
                    if "minRatio" in config.keys():
                        itle = QtGui.QTreeWidgetItem(["Min. ratio" , config["minRatio"]]);
                        itl.addChild(itle);
                        itle.myType = "parameter"
                    if "maxRatio" in config.keys():
                        itle = QtGui.QTreeWidgetItem(["Max. ratio" , config["maxRatio"]]);
                        itl.addChild(itle);
                        itle.myType = "parameter"

            itl = QtGui.QTreeWidgetItem(["MZ picking"]);
            it.addChild(itl);
            itl.myType = "parameter"
            if "positiveScanEvent" in config.keys():
                itle = QtGui.QTreeWidgetItem(["ScanEvent (positive)", config["positiveScanEvent"]]);
                itl.addChild(itle)
                itle.myType = "parameter"
            if "negativeScanEvent" in config.keys():
                itle = QtGui.QTreeWidgetItem(["ScanEvent (negative)", config["negativeScanEvent"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "xMin" in config.keys() and "xMax" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Atoms range", "%d - %d" % (int(config["xMin"]), int(config["xMax"]))]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "startTime" in config.keys() and "stopTime" in config.keys():
                itle = QtGui.QTreeWidgetItem(
                    ["Scan range", "%.2f - %.2f min" % (float(config["startTime"]), float(config["stopTime"]))]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "intensityThreshold" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Intensity threshold", "%.0f" % float(config["intensityThreshold"])]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "intensityCutoff" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Intensity cutoff", "%.0f" % float(config["intensityCutoff"])]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "maxLoading" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Max. charge", "%d" % int(config["maxLoading"])]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "ppm" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Mass deviation (+/- ppm)", "%.1f" % (float(config["ppm"]))]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "isotopicPatternCountLeft" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Isotopic pattern count (A)", config["isotopicPatternCountLeft"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "isotopicPatternCountRight" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Isotopic pattern count (B)", config["isotopicPatternCountRight"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "lowAbundanceIsotopeCutoff" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Consider isotopologue abundance", config["lowAbundanceIsotopeCutoff"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "Intensity abundance error" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Intensity abundance error"]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "intensityErrorN" in config.keys():
                itlee = QtGui.QTreeWidgetItem(["Non-labelled ion", config["intensityErrorN"]]);
                itle.addChild(itlee);
                itlee.myType = "parameter"
            if "intensityErrorL" in config.keys():
                itlee = QtGui.QTreeWidgetItem(["Labelled ion", config["intensityErrorL"]]);
                itle.addChild(itlee);
                itlee.myType = "parameter"

            itl = QtGui.QTreeWidgetItem(["MZ clustering"]);
            it.addChild(itl);
            itl.myType = "parameter"
            if "clustPPM" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Clustering ppm", config["clustPPM"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "minSpectraCount" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Min. spectra", config["minSpectraCount"]]);
                itl.addChild(itle);
                itle.myType = "parameter"

            itl = QtGui.QTreeWidgetItem(["Chromatographic peak picking"]);
            it.addChild(itl);
            itl.myType = "parameter"
            if "chromPeakPPM" in config.keys():
                itle = QtGui.QTreeWidgetItem(["EIC ppm", config["chromPeakPPM"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "eicSmoothing" in config.keys():
                itle = QtGui.QTreeWidgetItem(["EIC smoothing", config["eicSmoothing"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "eicSmoothing" in config.keys() and bool(config["eicSmoothing"]) and "eicSmoothingWindowSize" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Smoothing window size", config["eicSmoothingWindowSize"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "scales" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Min. scale", "%.0f" % (loads(base64.b64decode(config["scales"]))[0])]);
                itl.addChild(itle);
                itle.myType = "parameter"
                itle = QtGui.QTreeWidgetItem(["Max. scale", "%.0f" % (loads(base64.b64decode(config["scales"]))[1])]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "snrTh" in config.keys():
                itle = QtGui.QTreeWidgetItem(["SNR threshold", config["snrTh"]]);
                itl.addChild(itle);
                itle.myType = "parameter"

            itl = QtGui.QTreeWidgetItem(["Peak matching"]);
            it.addChild(itl);
            itl.myType = "parameter"
            if "peakCenterError" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Center error", config["peakCenterError"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "minPeakCorr" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Min. corr", config["minPeakCorr"]]);
                itl.addChild(itle);
                itle.myType = "parameter"

            itl = QtGui.QTreeWidgetItem(["Group features"]);
            it.addChild(itl);
            itl.myType = "parameter"
            if "minCorrelation" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Min. corr", config["minCorrelation"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            if "minCorrelationConnections" in config.keys():
                itle = QtGui.QTreeWidgetItem(["Min. correlation connections", config["minCorrelationConnections"]]);
                itl.addChild(itle);
                itle.myType = "parameter"
            itle = QtGui.QTreeWidgetItem(["Adducts"]);
            itl.addChild(itle);
            itle.myType = "parameter"
            if "adducts" in config.keys():
                for adduct in loads(base64.b64decode(config["adducts"])):
                    itlee = QtGui.QTreeWidgetItem([str(adduct.name), str(adduct.mzoffset), str(adduct.polarity)]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"
            itle = QtGui.QTreeWidgetItem(["Neutral loss (elements)"]);
            itl.addChild(itle);
            itle.myType = "parameter"
            if "elements" in config.keys():
                for elem, elemDetails in loads(base64.b64decode(config["elements"])).items():
                    itlee = QtGui.QTreeWidgetItem(
                        [str(elem), "%.4f" % elemDetails.weight, str(elemDetails.numberValenzElectrons),
                        "%d-%d" % (elemDetails.minCount, elemDetails.maxCount)]);
                    itle.addChild(itlee);
                    itlee.myType = "parameter"

            self.sortTreeChildren(self.ui.res_ExtractedData.topLevelItem(3), 1)

            self.filterEdited(str(self.ui.dataFilter.text()))

            sectionSizes=[235,65,50,105,110,70,85,70,50,60,75]
            for i in range(11):
                self.ui.res_ExtractedData.header().resizeSection(i, sectionSizes[i])




            ## Setup diagnostics
            it = QtGui.QTreeWidgetItem(["Diagnostics"]);
            it.myType = "diagnostic"
            self.ui.res_ExtractedData.addTopLevelItem(it)

            itl = QtGui.QTreeWidgetItem(["Observed intensities"]);
            it.addChild(itl);
            itl.myType = "diagnostic - observed intensities"

            itl = QtGui.QTreeWidgetItem(["Relative mz error"]);
            it.addChild(itl);
            itl.myType = "diagnostic - relative mz error"

            itl = QtGui.QTreeWidgetItem(["Relative mzbin deviation"]);
            it.addChild(itl);
            itl.myType = "diagnostic - relative mzbin deviation"

            itl = QtGui.QTreeWidgetItem(["Feature pair correlations"]);
            it.addChild(itl);
            itl.myType = "diagnostic - feature pair correlations"

            itl = QtGui.QTreeWidgetItem(["Feature pair M+1/M ratio"]);
            it.addChild(itl);
            itl.myType = "diagnostic - feature pair mp1 to m ratio"

            itl = QtGui.QTreeWidgetItem(["Feature pair M'-1/M' ratio"]);
            it.addChild(itl);
            itl.myType = "diagnostic - feature pair mPp1 to mP ratio"

            itl = QtGui.QTreeWidgetItem(["Feature pair assigned MZs"]);
            it.addChild(itl);
            itl.myType = "diagnostic - feature pair assigned mzs"

            itl = QtGui.QTreeWidgetItem(["Feature pair MZ deviation"]);
            it.addChild(itl);
            itl.myType = "diagnostic - feature pair mz deviation"




        else:
            QtGui.QMessageBox.warning(self, "MetExtract", "No results are available for file %s", QtGui.QMessageBox.Ok)

    def deColorQTreeWidgetItem(self, item):
        item.setBackgroundColor(0, QColor("white"))
        item.setBackgroundColor(1, QColor("white"))
        item.setBackgroundColor(2, QColor("white"))
        item.setBackgroundColor(3, QColor("white"))
        item.setBackgroundColor(4, QColor("white"))
        item.setBackgroundColor(5, QColor("white"))
        item.setBackgroundColor(6, QColor("white"))
        item.setBackgroundColor(7, QColor("white"))
        item.setBackgroundColor(8, QColor("white"))
        item.setBackgroundColor(9, QColor("white"))

        for i in range(item.childCount()):
            self.deColorQTreeWidgetItem(item.child(i))

    def getParametersFromCurrentRes(self, parameter):
        p = self.ui.res_ExtractedData.topLevelItem(4)
        childs = [p]
        while len(childs) > 0:
            child = childs.pop(0)
            for i in range(child.childCount()):
                childs.append(child.child(i))
            if child.text(0) == parameter:
                return child.text(1)

    def isTracerMetabolisationExperiment(self):
        rows = []
        for row in self.currentOpenResultsFile.curs.execute(
                "SELECT key, value FROM config WHERE key='metabolisationExperiment'"):
            rows.append(copy(row))
        assert len(rows) == 1
        return str(rows[0][1]).lower() == "true"

    def getAllowedIsotopeRatioErrorsForResult(self):
        rows = []
        for row in self.currentOpenResultsFile.curs.execute(
                "SELECT key, value FROM config WHERE key='intensityErrorN'"):
            rows.append(copy(row))
        assert len(rows) == 1
        intErrN = float(rows[0][1])

        rows = []
        for row in self.currentOpenResultsFile.curs.execute(
                "SELECT key, value FROM config WHERE key='intensityErrorL'"):
            rows.append(copy(row))
        assert len(rows) == 1
        intErrL = float(rows[0][1])

        return intErrN, intErrL

    def getTICsForResult(self):

        if hasattr(self.currentOpenResultsFile, "tics"):
            return self.currentOpenResultsFile.tics

        tics={}
        for row in self.currentOpenResultsFile.curs.execute("SELECT id, loading, scanevent, times, intensities FROM tics"):
            id, loading, scanevent, times, intensities=row
            id=int(id)
            loading=str(loading)
            scanevent=str(scanevent)
            times=[float(t)/60. for t in times.split(";")]
            intensities=[float(t) for t in intensities.split(";")]

            tics[loading]=Bunch(id=id, loading=loading, scanEvent=scanevent, times=times, intensities=intensities)

        setattr(self.currentOpenResultsFile, "tics", tics)

        return tics


    def getLabellingParametersForResult(self, featureID):
        if not (self.isTracerMetabolisationExperiment()):
            rows = []
            for row in self.currentOpenResultsFile.curs.execute("SELECT key, value FROM config WHERE key='isotopeA'"):
                rows.append(copy(row))
            assert len(rows) == 1
            isoA = float(rows[0][1])

            rows = []
            for row in self.currentOpenResultsFile.curs.execute("SELECT key, value FROM config WHERE key='isotopeB'"):
                rows.append(copy(row))
            assert len(rows) == 1
            isoB = float(rows[0][1])

            rows = []
            for row in self.currentOpenResultsFile.curs.execute("SELECT key, value FROM config WHERE key='purityN'"):
                rows.append(copy(row))
            assert len(rows) == 1
            purN = float(rows[0][1])

            rows = []
            for row in self.currentOpenResultsFile.curs.execute("SELECT key, value FROM config WHERE key='purityL'"):
                rows.append(copy(row))
            assert len(rows) == 1
            purL = float(rows[0][1])

            return isoB - isoA, purN, purL
        else:
            rows = []
            for row in self.currentOpenResultsFile.curs.execute("SELECT id, tracer FROM chromPeaks WHERE id=%d" % featureID):
                rows.append(copy(row))
            assert len(rows) == 1
            trcid = int(rows[0][1])

            rows = []
            for row in self.currentOpenResultsFile.curs.execute("SELECT deltaMZ, purityN, purityL FROM tracerConfiguration WHERE id=%d" % trcid):
                rows.append(copy(row))
            assert len(rows) == 1
            dmz = float(rows[0][0])
            purN = float(rows[0][1])
            purL = float(rows[0][2])

            return dmz, purN, purL

    def selectedResChanged(self):

        for i in range(self.ui.res_ExtractedData.topLevelItemCount()):
            self.deColorQTreeWidgetItem(self.ui.res_ExtractedData.topLevelItem(i))

        selectedItems = self.ui.res_ExtractedData.selectedItems()

        changePlots = False
        for item in selectedItems:
            if hasattr(item, "myType"):
                if item.myType == "parameter":
                    continue
                else:
                    changePlots = True
        if not changePlots:
            return

        self.clearPlot(self.ui.pl1)
        self.clearPlot(self.ui.pl2)
        self.clearPlot(self.ui.pl3)

        useColi = 0
        maxIntX = 0
        maxIntY = 0
        minIntX = 0
        minIntY = 0
        minTime = 1000
        maxTime = 1
        x_vals = []
        y_vals = []
        mzs = []
        peaks = []
        plotTypes = set()
        selIndex = 0
        selFeatureGroups = []

        featuresPosSelected = False
        for item in selectedItems:
            if not (hasattr(item, "myType")):
                continue

            if item.myType == "MZs" or item.myType == "mz":
                plotTypes.add("MZs")
            if item.myType == "MZBins" or item.myType == "mzbin":
                plotTypes.add("MZBins")
            if item.myType == "Features" or item.myType == "feature":
                plotTypes.add("Features")
            if item.myType == "Feature Groups" or item.myType == "feature group":
                plotTypes.add("Feature Groups")
            if item.myType.lower().startswith("diagnostic"):
                plotTypes.add("diagnostic")


            if item.myType == "Features" or item.myType == "feature":
                if hasattr(item, "myData"):
                    cp = item.myData
                    if cp.ionMode == "-":
                        pass
                    else:
                        featuresPosSelected = True

        if len(plotTypes) > 1:
            QtGui.QMessageBox.warning(self, "MetExtract","Selecting different result types in not supported. Please select only one or multipel MZs, MZBins, FeaturePairs or FeatureGroups at a time",QtGui.QMessageBox.Ok)
            self.clearPlot(self.ui.pl1)
            return

        for item in selectedItems:
            if not (hasattr(item, "myType")):
                continue

            #<editor-fold desc="#mz result">
            if item.myType == "MZs" or item.myType == "mz":
                plotTypes.add("MZs")
                self.ui.res_ExtractedData.setHeaderLabels(QtCore.QStringList(["MZ", "Xn", "Scan id", "Rt", "Charge", "Intensity", "", "", "", "", ""]))

                t = item
                if len(x_vals) == 0:
                    x = []
                    y = []
                    if t.myType == "mz":
                        t = t.parent()
                    for j in range(t.childCount()):
                        child = t.child(j)
                        assert child.myType == "mz"
                        x.append(child.myData.scantime / 60.)
                        y.append(child.myData.mz)

                    maxIntY = max(y)
                    maxIntX = max(x)

                    x_vals = [x]
                    y_vals = [y]

                if item.myType == "mz":
                    if len(x_vals) == 1:
                        x_vals.append([])
                        y_vals.append([])
                    x_vals[1].append(item.myData.scantime / 60.)
                    y_vals[1].append(item.myData.mz)
            #</editor-fold>

            #<editor-fold desc="#mzbin results">
            elif item.myType == "MZBins" or item.myType == "mzbin":
                self.ui.res_ExtractedData.setHeaderLabels(QtCore.QStringList(["MZ", "Delta ppm", "Xn", "", "", "", "", "", "", "", ""]))
                if item.myType == "MZBins":
                    plotTypes.add("MZBins")
                    for i in range(item.childCount()):
                        kid = item.child(i)
                        assert kid.myType == "mzbin"
                        x = []
                        y = []
                        for j in range(kid.childCount()):
                            child = kid.child(j)
                            assert child.myType == "mz"
                            x.append(child.myData.scantime / 60.)
                            y.append(child.myData.mz)
                        x_vals.append(x)
                        y_vals.append(y)
                    if(len(x_vals)>0):
                        maxIntX = max(max(x_vals), maxIntX)
                        maxIntY = max(max(y_vals), maxIntY)
                    else:
                        maxIntX=1
                        maxIntY=1

                elif item.myType == "mzbin":
                    plotTypes.add("mzbin")
                    x = []
                    y = []
                    if item.myType == "mzbin" and len(x_vals) == 0:
                        t = item.parent()
                        for o in range(t.childCount()):
                            childo = t.child(o)
                            for j in range(childo.childCount()):
                                childj = childo.child(j)

                                assert childj.myType == "mz"
                                x.append(childj.myData.scantime / 60.)
                                y.append(childj.myData.mz)
                        x_vals.append(x)
                        y_vals.append(y)

                    x = []
                    y = []
                    for j in range(item.childCount()):
                        child = item.child(j)
                        assert child.myType == "mz"
                        x.append(child.myData.scantime / 60.)
                        y.append(child.myData.mz)
                    x_vals.append(x)
                    y_vals.append(y)
                if len(x_vals)>0:
                    maxIntX = max(max(x_vals), maxIntX)
                    maxIntY = max(max(y_vals), maxIntY)
                else:
                    maxIntX=1
                    maxIntY=1
            #</editor-fold>

            #<editor-fold desc="#feature results">
            elif item.myType == "Features" or item.myType == "feature":
                self.ui.chromPeakName.setText("")

                self.ui.res_ExtractedData.setHeaderLabels(QtCore.QStringList(
                    ["MZ (/Ionmode Z)", "Rt min", "Xn", "Adducts, hetero atoms", "Scale M / M'", "Peak cor", "M:M' peaks ratio / area ratio", "Area M / M'", "Scans", "Tracer"]))

                if item.myType == "Features":
                    mzs = []
                    rts = []
                    plotTypes.add("Features")
                    for i in range(item.childCount()):
                        child = item.child(i)
                        assert child.myType == "feature"
                        mzs.append(child.myData.mz);
                        rts.append(child.myData.NPeakCenterMin / 60.)
                    self.drawPlot(self.ui.pl1, plotIndex=0, x=rts, y=mzs, ylab="m/z", useCol=0, scatter=True,
                                  plot=False)

                if item.myType == "feature":
                    cp = item.myData
                    plotTypes.add("feature")
                    mzs.append(cp.mz)
                    peaks.append(cp.NPeakCenterMin / 60.)
                    xic = []
                    xicL = []
                    times = []
                    maxE=1

                    invert=1
                    if self.ui.negEIC.isChecked():
                        invert=-1

                    for row in self.currentOpenResultsFile.curs.execute(
                                    "SELECT xic, xicL, xicfirstiso, xicLfirstiso, xicLfirstisoconjugate, xic_smoothed, xicL_smoothed, times, allPeaks FROM XICs WHERE id==%d" % cp.eicID):
                        xic                   = [float(t) for t in row[0].split(";")]
                        xicL                  = [float(t) for t in row[1].split(";")]
                        xicfirstiso           = [float(t) for t in row[2].split(";")]
                        xicLfirstiso          = [float(t) for t in row[3].split(";")]
                        xicLfirstisoconjugate = [float(t) for t in row[4].split(";")]
                        xic_smoothed                   = [float(t) for t in row[5].split(";")]
                        xicL_smoothed                  = [float(t) for t in row[6].split(";")]
                        times = [float(t) / 60. for t in row[7].split(";")]
                        allPeaks = loads(base64.b64decode(str(row[8])))


                    if self.ui.scaleFeatures.isChecked():
                        s = int(cp.NPeakCenter - cp.NBorderLeft * 1)
                        e = int(cp.NPeakCenter + cp.NBorderRight * 1)
                        maxP = s + max(range(e - s), key=lambda x: xic[s + x])
                        maxE  = mean([xic[maxP - 3], xic[maxP - 2], xic[maxP - 1], xic[maxP], xic[maxP + 1], xic[maxP + 2],xic[maxP + 3]])
                        maxEL = mean([xicL[maxP - 3], xicL[maxP - 2], xicL[maxP - 1], xicL[maxP], xicL[maxP + 1], xicL[maxP + 2], xicL[maxP + 3]])
                        if maxE == 0:
                            maxE=1
                        if maxEL == 0:
                            maxEL=1
                        if not self.ui.scaleLabelledFeatures.isChecked():
                            maxEL=maxE

                        xic                   = [u / maxE for u in xic]
                        xicfirstiso           = [u / maxE for u in xicfirstiso]
                        xicL                  = [u / maxEL for u in xicL]
                        xicLfirstiso          = [u / maxEL for u in xicLfirstiso]
                        xicLfirstisoconjugate = [u / maxEL for u in xicLfirstisoconjugate]
                        xic_smoothed                   = [u / maxE for u in xic_smoothed]
                        xicL_smoothed                  = [u / maxEL for u in xicL_smoothed]

                    xicL                  = [invert*u for u in xicL]
                    xicLfirstiso          = [invert*u for u in xicLfirstiso]
                    xicLfirstisoconjugate = [invert*u for u in xicLfirstisoconjugate]
                    xicL_smoothed                  = [invert*u for u in xicL_smoothed]

                    if self.ui.flattenXIC.isChecked():
                        ps = int(cp.NPeakCenter - cp.NPeakScale * 2)
                        pe = int(cp.NPeakCenter + cp.NPeakScale * 2)
                        for u in range(len(xic)):
                            if u < ps or u > pe:
                                xic[u] = 0
                                xicfirstiso[u] = 0
                        ps = int(cp.LPeakCenter - cp.LPeakScale * 2)
                        pe = int(cp.LPeakCenter + cp.LPeakScale * 2)
                        for u in range(len(xicL)):
                            if u < ps or u > pe:
                                xicL[u] = 0
                                xicLfirstiso[u] = 0
                                xicL_smoothed[u] = 0
                                xicLfirstiso[u] = 0
                    try:
                        minTime = min(minTime, min(
                            times[int(cp.NPeakCenter - cp.NPeakScale * 1): int(cp.NPeakCenter + cp.NPeakScale * 1)]))
                        maxTime = max(maxTime, max(
                            times[int(cp.NPeakCenter - cp.NPeakScale * 1): int(cp.NPeakCenter + cp.NPeakScale * 1)]))
                        maxIntY = max(maxIntY, max(
                            xic[int(cp.NPeakCenter - cp.NPeakScale * 1): int(cp.NPeakCenter + cp.NPeakScale * 1)]))
                        minIntY = min(minIntY, min(
                            xicL[int(cp.LPeakCenter - cp.LPeakScale * 1): int(cp.LPeakCenter + cp.LPeakScale * 1)]))
                    except:
                        pass

                    item.setBackgroundColor(0, QColor(predefinedColors[useColi % len(predefinedColors)]))
                    item.setBackgroundColor(1, QColor(predefinedColors[useColi % len(predefinedColors)]))
                    item.setBackgroundColor(2, QColor(predefinedColors[useColi % len(predefinedColors)]))
                    item.setBackgroundColor(3, QColor(predefinedColors[useColi % len(predefinedColors)]))
                    item.setBackgroundColor(4, QColor(predefinedColors[useColi % len(predefinedColors)]))
                    item.setBackgroundColor(5, QColor(predefinedColors[useColi % len(predefinedColors)]))
                    item.setBackgroundColor(6, QColor(predefinedColors[useColi % len(predefinedColors)]))
                    item.setBackgroundColor(7, QColor(predefinedColors[useColi % len(predefinedColors)]))

                    self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xic,
                                  fill=[int(cp.NPeakCenter - cp.NBorderLeft),
                                        int(cp.NPeakCenter + cp.NBorderRight)], rearrange=len(selectedItems) == 1,
                                  label="%.4f (%d)"%(cp.mz, cp.xCount), useCol=useColi)

                    if self.ui.showSmoothedEIC_checkBox.isChecked():
                        self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xic_smoothed,
                                      fill=[int(cp.NPeakCenter - cp.NBorderLeft),
                                            int(cp.NPeakCenter + cp.NBorderRight)], rearrange=len(selectedItems) == 1,
                                      label=None, useCol=useColi, linestyle="--")


                    if self.ui.showIsotopologues.isChecked():
                        self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xicfirstiso,
                                      fill=[int(cp.NPeakCenter - cp.NBorderLeft),
                                            int(cp.NPeakCenter + cp.NBorderRight)], rearrange=len(selectedItems) == 1,
                                      label=None, useCol=useColi)

                    self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xicL,
                                  fill=[int(cp.LPeakCenter - cp.LBorderLeft),
                                        int(cp.LPeakCenter + cp.LBorderRight)], rearrange=len(selectedItems) == 1,
                                  label=None, useCol=useColi)

                    if self.ui.showSmoothedEIC_checkBox.isChecked():
                        self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xicL_smoothed,
                                      fill=[int(cp.LPeakCenter - cp.LBorderLeft),
                                            int(cp.LPeakCenter + cp.LBorderRight)], rearrange=len(selectedItems) == 1,
                                      label=None, useCol=useColi, linestyle="--")

                    if self.ui.showIsotopologues.isChecked():
                        self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xicLfirstiso,
                                      fill=[int(cp.LPeakCenter - cp.LBorderLeft),
                                            int(cp.LPeakCenter + cp.LBorderRight)], rearrange=len(selectedItems) == 1,
                                      label=None, useCol=useColi)

                    maxindex, maxvalue = max(
                        enumerate(xic[int(cp.NPeakCenter - 1):int(cp.NPeakCenter + 1)], start=int(cp.NPeakCenter - 1)),
                        key=itemgetter(1))
                    if self.ui.plotAddLabels.checkState() == QtCore.Qt.Checked:
                        self.addAnnotation(self.ui.pl1, "%.5f\n%.2f min" % (cp.mz, cp.NPeakCenterMin / 60.),
                                           (times[maxindex], xic[maxindex]), (times[maxindex], xic[maxindex]), 0,
                                           fcColor=predefinedColors[useColi % len(predefinedColors)],
                                           ecColor=predefinedColors[useColi % len(predefinedColors)],
                                           arrowColor=predefinedColors[useColi % len(predefinedColors)], alpha=0.25)



                    if self.ui.showDiagnostics.isChecked():
                        for pa in allPeaks["peaksN"]:
                            if pa.peakIndex!=cp.NPeakCenter:
                                self.addAnnotation(self.ui.pl1, "" ,
                                                   (times[pa.peakIndex], xic[pa.peakIndex]), (times[pa.peakIndex], xic[pa.peakIndex]), 0,
                                                   fcColor="slategrey",
                                                   ecColor="slategrey",
                                                   arrowColor="slategrey", alpha=0.15, add=30)
                        for pa in allPeaks["peaksL"]:
                            self.addAnnotation(self.ui.pl1, "" ,
                                               (times[pa.peakIndex], xicL[pa.peakIndex]), (times[pa.peakIndex], xicL[pa.peakIndex]), 0,
                                               fcColor="slategrey",
                                               ecColor="slategrey",
                                               arrowColor="slategrey", alpha=0.15, up=False, add=30)



                        for row in self.currentOpenResultsFile.curs.execute(
                                    "SELECT c.id, c.eicID, c.NPeakCenterMin, c.NPeakCenter, c.mz, c.xcount, c.loading, (SELECT fg.featureName FROM FeatureGroupFeatures fgf INNER JOIN FeatureGroups fg ON fgf.fGroupID=fg.id WHERE fgf.fID=c.id) AS FGroupID FROM chromPeaks c WHERE c.eicID IN (%d)" % cp.eicID):

                            if cp.NPeakCenter!=row[3]:
                                self.addAnnotation(self.ui.pl1, "%s\n%.5f\n%.2f min"%(str(row[7]), row[4], row[2]/60.) ,
                                                   (times[row[3]], xic[row[3]]), (times[row[3]], xic[row[3]]), 0,
                                                   fcColor="slategrey",
                                                   ecColor="slategrey",
                                                   arrowColor="slategrey", alpha=0.15)


                    isMetabolisationExperiment = self.isTracerMetabolisationExperiment()
                    mzD, purN, purL = self.getLabellingParametersForResult(cp.id)

                    toDrawMzs = []
                    toDrawInts = []
                    for row in self.currentOpenResultsFile.curs.execute("SELECT mzs, intensities, ionmode FROM massspectrum WHERE mID=%d"%cp.massSpectrumID):

                        mzs = [float(t) for t in row[0].split(";")]
                        intensities = [float(t) for t in row[1].split(";")]

                        for mz in mzs:
                            toDrawMzs.append(mz)
                            toDrawMzs.append(mz)
                            toDrawMzs.append(mz)

                        for intensity in intensities:
                            toDrawInts.append(0)
                            if cp.ionMode == "-" and featuresPosSelected:
                                toDrawInts.append(-intensity)
                            else:
                                toDrawInts.append(intensity)
                            toDrawInts.append(0)

                    self.drawPlot(self.ui.pl3, plotIndex=0, x=toDrawMzs, y=toDrawInts, useCol=useColi, multipleLocator=None,
                                  alpha=0.1, title="", xlab="MZ")

                    bm = min(range(len(toDrawMzs)), key=lambda i: abs(toDrawMzs[i] - cp.mz)) + 1
                    bml = min(range(len(toDrawMzs)),
                              key=lambda i: abs(toDrawMzs[i] - (cp.mz + mzD * cp.xCount / cp.loading))) + 1

                    intLeft = toDrawInts[bm]
                    intRight = toDrawInts[bml]

                    h = 0
                    if cp.ionMode == "-" and featuresPosSelected:
                        h = min(intLeft, intRight)
                    else:
                        h = max(intLeft, intRight)

                    if self.ui.MSLabels.checkState() == QtCore.Qt.Checked:
                        self.addAnnotation(self.ui.pl3, "mz: %.5f\nl-mz: %.5f\nd-mz: %.5f\nXn: %d Z: %s%d" % (
                            cp.mz, cp.mz + mzD * cp.xCount / cp.loading, mzD * cp.xCount, cp.xCount, cp.ionMode,
                            cp.loading), (cp.mz + mzD * cp.xCount / cp.loading / 2., h), (10, 120), rotation=0,
                                           up=not (cp.ionMode == "-" and featuresPosSelected))

                    self.addArrow(self.ui.pl3, (cp.mz, toDrawInts[bm]), (cp.mz, h), drawArrowHead=True)
                    self.addArrow(self.ui.pl3, (cp.mz, h), (cp.mz + mzD * cp.xCount / cp.loading, h),
                                  ecColor="slategrey")
                    self.addArrow(self.ui.pl3, (cp.mz + mzD * cp.xCount / cp.loading, toDrawInts[bml]),(cp.mz + mzD * cp.xCount / cp.loading, h), drawArrowHead=True)

                    annotationHeight=h

                    if self.ui.MSIsos.checkState() == QtCore.Qt.Checked:
                        bml = min(range(len(toDrawMzs)),
                                  key=lambda w: abs(toDrawMzs[w] - (cp.mz + 1.00335 / cp.loading))) + 1

                        if cp.ionMode == "-" and featuresPosSelected:
                            h = min(toDrawInts[bm], toDrawInts[bml])
                        else:
                            h = max(toDrawInts[bm], toDrawInts[bml])

                        self.addArrow(self.ui.pl3, (cp.mz + 1.00335 / cp.loading, toDrawInts[bml]),
                                      (cp.mz + 1.00335 / cp.loading, annotationHeight), drawArrowHead=True)

                        bm = min(range(len(toDrawMzs)),
                                 key=lambda w: abs(toDrawMzs[w] - (cp.mz + 1.00335 * (cp.xCount - 1) / cp.loading))) + 1
                        bml = min(range(len(toDrawMzs)),
                                  key=lambda w: abs(toDrawMzs[w] - (cp.mz + 1.00335 * cp.xCount / cp.loading))) + 1

                        if cp.ionMode == "-" and featuresPosSelected:
                            h = min(toDrawInts[bm], toDrawInts[bml])
                        else:
                            h = max(toDrawInts[bm], toDrawInts[bml])

                        self.addArrow(self.ui.pl3,
                                      (cp.mz + 1.00335 * (cp.xCount - 1) / cp.loading, toDrawInts[bm]),
                                      (cp.mz + 1.00335 * (cp.xCount - 1) / cp.loading, annotationHeight), drawArrowHead=True)


                        self.addArrow(self.ui.pl3, (cp.mz, 0), (cp.mz, intLeft), linewidth=5, ecColor="orange")
                        self.addArrow(self.ui.pl3, (cp.mz + 1.00335 * cp.xCount / cp.loading, 0),
                                      (cp.mz + 1.00335 * cp.xCount / cp.loading, intRight), linewidth=5,
                                      ecColor="orange")

                        intErrN, intErrL = self.getAllowedIsotopeRatioErrorsForResult()


                        for iso in [1, 2, 3]:
                            ratioN = getNormRatio(purN, cp.xCount, iso)
                            ratioL = getNormRatio(purL, cp.xCount, iso)
                            self.addArrow(self.ui.pl3, (cp.mz + (1.00335 * iso) / cp.loading, 0), (
                            cp.mz + (1.00335 * iso) / cp.loading, intLeft * max(0, (ratioN - intErrN))),
                                          linewidth=5, alpha=.1, ecColor="orange")
                            self.addArrow(self.ui.pl3, (
                            cp.mz + (1.00335 * iso) / cp.loading, intLeft * max(0, (ratioN - intErrN))),
                                          (cp.mz + (1.00335 * iso) / cp.loading,
                                          intLeft * (ratioN + intErrN)), linewidth=5, alpha=.1,
                                          ecColor="DarkSeaGreen")
                            self.addArrow(self.ui.pl3, (
                            cp.mz + (1.00335 * iso) / cp.loading, intLeft * max(0, (ratioN - .005))), (
                                          cp.mz + (1.00335 * iso) / cp.loading,
                                          intLeft * max(0, (ratioN + .005))), linewidth=5, alpha=.1,
                                          ecColor="yellow")

                            self.addArrow(self.ui.pl3, (cp.mz + 1.00335 * (cp.xCount - iso) / cp.loading, 0), (
                                cp.mz + 1.00335 * (cp.xCount - iso) / cp.loading,
                                intRight * max(0, (ratioL - intErrL))), linewidth=5, alpha=.1,
                                          ecColor="orange")
                            self.addArrow(self.ui.pl3, (cp.mz + 1.00335 * (cp.xCount - iso) / cp.loading,
                                                        intRight * max(0, (ratioL - intErrL))), (
                                              cp.mz + 1.00335 * (cp.xCount - iso) / cp.loading,
                                              intRight * (ratioL + intErrL)), linewidth=5, alpha=.1,
                                          ecColor="DarkSeaGreen")
                            self.addArrow(self.ui.pl3, (cp.mz + 1.00335 * (cp.xCount - iso) / cp.loading,
                                                        intRight * max(0, (ratioL - .005))), (
                                              cp.mz + 1.00335 * (cp.xCount - iso) / cp.loading,
                                              intRight * max(0, (ratioL + .005))), linewidth=5,
                                          alpha=.1, ecColor="yellow")

                        if self.ui.drawFPIsotopologues.checkState() == QtCore.Qt.Checked:
                            for iso in [1, 2, 3]:
                                self.addArrow(self.ui.pl3, (cp.mz - (1.00335 * iso) / cp.loading, 0),
                                              (cp.mz - (1.00335 * iso) / cp.loading, intLeft * .1), linewidth=5,
                                              ecColor="yellow")
                                self.addArrow(self.ui.pl3, (cp.mz + 1.00335 * (cp.xCount + iso) / cp.loading, 0),
                                              (cp.mz + 1.00335 * (cp.xCount + iso) / cp.loading, intRight * .1),
                                              linewidth=5, ecColor="yellow")

                            self.addArrow(self.ui.pl3, (cp.mz - 1.00335 / (cp.loading * 2), 0),
                                          (cp.mz - 1.00335 / (cp.loading * 2), intLeft * .05), linewidth=5,
                                          ecColor="yellow")
                            self.addArrow(self.ui.pl3, (cp.mz + 1.00335 / (cp.loading * 2), 0),
                                          (cp.mz + 1.00335 / (cp.loading * 2), intLeft * .05), linewidth=5,
                                          ecColor="yellow")
                            self.addArrow(self.ui.pl3, (
                                cp.mz + 1.00335 * cp.xCount / cp.loading + 1.00335 / (cp.loading * 2), 0), (
                                              cp.mz + 1.00335 * cp.xCount / cp.loading + 1.00335 / (cp.loading * 2),
                                              intRight * .05), linewidth=5, ecColor="yellow")
                            self.addArrow(self.ui.pl3, (
                                cp.mz + 1.00335 * cp.xCount / cp.loading - 1.00335 / (cp.loading * 2), 0), (
                                              cp.mz + 1.00335 * cp.xCount / cp.loading - 1.00335 / (cp.loading * 2),
                                              intRight * .05), linewidth=5, ecColor="yellow")
                useColi += 1
            #</editor-fold>

            #<editor-fold desc="#featureGroup results">
            elif item.myType == "featureGroup" or item.myType == "Feature Groups":
                self.ui.res_ExtractedData.setHeaderLabels(
                    QtCore.QStringList(["Name", "# Peaks", "Rt", "Adducts", "Heteroatoms", "M Scale", "M' Rt (min)", "M' Scale", "Corr", "Tracer", "Ratio"]))


                item.setBackgroundColor(0, QColor(predefinedColors[(useColi) % len(predefinedColors)]))
                item.setBackgroundColor(1, QColor(predefinedColors[(useColi) % len(predefinedColors)]))
                item.setBackgroundColor(2, QColor(predefinedColors[(useColi) % len(predefinedColors)]))
                item.setBackgroundColor(3, QColor(predefinedColors[(useColi) % len(predefinedColors)]))
                item.setBackgroundColor(4, QColor(predefinedColors[(useColi) % len(predefinedColors)]))
                item.setBackgroundColor(5, QColor(predefinedColors[(useColi) % len(predefinedColors)]))
                item.setBackgroundColor(6, QColor(predefinedColors[(useColi) % len(predefinedColors)]))
                item.setBackgroundColor(7, QColor(predefinedColors[(useColi) % len(predefinedColors)]))

                if item.myType == "Feature Groups":
                    plotTypes.add("Feature Groups")
                    for i in range(item.childCount()):
                        child = item.child(i)
                        mzs = [];
                        rts = []
                        for j in range(child.childCount()):
                            feature = child.child(j)
                            assert feature.myType == "feature"

                            mzs.append(feature.myData.mz);
                            rts.append(feature.myData.NPeakCenterMin / 60.)
                        self.drawPlot(self.ui.pl1, plotIndex=0, x=rts, y=mzs, ylab="m/z", useCol=i, scatter=True,plot=True)

                if item.myType == "featureGroup":
                    selFeatureGroups.append(item)
                    self.ui.chromPeakName.setText(item.myData.featureName)
                    plotTypes.add("feature")
                    meanRT = 0
                    countK = 0

                    self.clearPlot(self.ui.pl3)

                    toDrawIntsPos = []
                    toDrawIntsNeg = []
                    toDrawMZsPos = []
                    toDrawMZsNeg = []
                    hasPos = False
                    hasNeg = False
                    massSpectraAvailable = False


                    msi = 0
                    mzs = {}
                    intensities = {}
                    mstime = {}

                    try:
                        for msspectrum in SQLSelectAsObject(self.currentOpenResultsFile.curs, "SELECT mzs, intensities, time AS tim, ionMode FROM massspectrum WHERE fgID=%d" % item.myID):
                            msi += 1
                            ionMode = str(msspectrum.ionMode)
                            mzs[ionMode] = [float(u) for u in str(msspectrum.mzs).split(";")]
                            intensities[ionMode] = [float(u) for u in str(msspectrum.intensities).split(";")]
                            mstime[ionMode] = float(msspectrum.tim)
                    except:
                        pass


                    hasPos = "+" in mzs.keys()
                    hasNeg = "-" in mzs.keys()

                    massSpectraAvailable = hasPos or hasNeg

                    if hasPos:
                        for mz in mzs["+"]:
                            toDrawMZsPos.append(mz)
                            toDrawMZsPos.append(mz)
                            toDrawMZsPos.append(mz)
                        for intensity in intensities["+"]:
                            toDrawIntsPos.append(0)
                            toDrawIntsPos.append(intensity)
                            toDrawIntsPos.append(0)

                        self.drawPlot(self.ui.pl3, plotIndex=0, x=toDrawMZsPos, y=toDrawIntsPos, useCol=useColi,
                                      multipleLocator=None, alpha=0.1, title="", xlab="MZ")

                    if hasNeg:
                        negInt = 1.

                        if "+" in mzs.keys():
                            negInt = -1.

                        for mz in mzs["-"]:
                            toDrawMZsNeg.append(mz)
                            toDrawMZsNeg.append(mz)
                            toDrawMZsNeg.append(mz)
                        for intensity in intensities["-"]:
                            toDrawIntsNeg.append(0)
                            toDrawIntsNeg.append(intensity * negInt)
                            toDrawIntsNeg.append(0)

                        self.drawPlot(self.ui.pl3, plotIndex=0, x=toDrawMZsNeg, y=toDrawIntsNeg, useCol=useColi,
                                      multipleLocator=None, alpha=0.1, title="", xlab="MZ")
                    mzs=[]
                    childIDs = []
                    maxInt = 0
                    for childi in range(item.childCount()):
                        if not (item.child(childi).isHidden()):

                            child = item.child(childi).myData

                            if massSpectraAvailable:

                                intMul = 1.
                                toDrawInts = []
                                toDrawMzs = []

                                if child.ionMode == "-" and hasPos:
                                    mInt = -min(toDrawIntsNeg)
                                    toDrawInts = [-f for f in toDrawIntsNeg]
                                    toDrawMzs = toDrawMZsNeg
                                    intMul = -1.
                                elif child.ionMode == "-":
                                    mInt = max(toDrawIntsNeg)
                                    toDrawInts = toDrawIntsNeg
                                    toDrawMzs = toDrawMZsNeg
                                else:
                                    mInt = max(toDrawIntsPos)
                                    toDrawInts = toDrawIntsPos
                                    toDrawMzs = toDrawMZsPos


                                mzD, purN, purL = self.getLabellingParametersForResult(child.id)

                                bm = min(range(len(toDrawMzs)), key=lambda i: abs(toDrawMzs[i] - child.mz)) + 1
                                bml = min(range(len(toDrawMzs)), key=lambda i: abs(
                                    toDrawMzs[i] - (child.mz + mzD * child.xCount / child.loading))) + 1

                                intLeft = toDrawInts[bm]
                                intRight = toDrawInts[bml]

                                h = max(toDrawInts[bm], toDrawInts[bml])

                                h = h
                                if self.ui.MSLabels.checkState() == QtCore.Qt.Checked:
                                    self.addAnnotation(self.ui.pl3,
                                                       "mz: %.5f\nl-mz: %.5f\nd-mz: %.5f\nXn: %d Z: %s%d" % (
                                                           child.mz, child.mz + mzD * child.xCount / child.loading,
                                                           mzD * child.xCount, child.xCount, child.ionMode,
                                                           child.loading),
                                                       (child.mz + mzD * child.xCount / child.loading / 2., h * intMul),
                                                       (10, 120), rotation=0, offset=(-10, 20), up=intMul > 0)

                                self.addArrow(self.ui.pl3, (child.mz, (toDrawInts[bm]) * intMul),
                                              (child.mz, h * intMul), drawArrowHead=True)
                                self.addArrow(self.ui.pl3, (child.mz, h * intMul),
                                              (child.mz + mzD * child.xCount / child.loading, h * intMul),
                                              ecColor="slategrey")
                                self.addArrow(self.ui.pl3, (
                                child.mz + mzD * child.xCount / child.loading, (toDrawInts[bml]) * intMul),
                                              (child.mz + mzD * child.xCount / child.loading, h * intMul),
                                              drawArrowHead=True)

                                if self.ui.MSIsos.checkState() == QtCore.Qt.Checked:
                                    bml = min(range(len(toDrawMzs)), key=lambda i: abs(
                                        toDrawMzs[i] - (child.mz + 1.00335 / child.loading))) + 1
                                    h = max(toDrawInts[bm], toDrawInts[bml])
                                    self.addArrow(self.ui.pl3, (child.mz, h * intMul),
                                                  (child.mz + 1.00335 / child.loading, h * intMul), ecColor="slategrey")
                                    self.addArrow(self.ui.pl3,
                                                  (child.mz + 1.00335 / child.loading, (toDrawInts[bml]) * intMul),
                                                  (child.mz + 1.00335 / child.loading, h * intMul), drawArrowHead=True)

                                    bm = min(range(len(toDrawMzs)), key=lambda i: abs(
                                        toDrawMzs[i] - (child.mz + 1.00335 * (child.xCount - 1) / child.loading))) + 1
                                    bml = min(range(len(toDrawMzs)), key=lambda i: abs(
                                        toDrawMzs[i] - (child.mz + 1.00335 * child.xCount / child.loading))) + 1
                                    h = max(toDrawInts[bm], toDrawInts[bml])
                                    self.addArrow(self.ui.pl3, (
                                        child.mz + 1.00335 * (child.xCount - 1) / child.loading,
                                        (toDrawInts[bm]) * intMul),
                                                  (child.mz + 1.00335 * (child.xCount - 1) / child.loading, h * intMul),
                                                  drawArrowHead=True)
                                    self.addArrow(self.ui.pl3,
                                                  (child.mz + 1.00335 * (child.xCount - 1) / child.loading, h * intMul),
                                                  (child.mz + 1.00335 * child.xCount / child.loading, h * intMul),
                                                  ecColor="slategrey")

                                    self.addArrow(self.ui.pl3, (child.mz, 0), (child.mz, intLeft * intMul), linewidth=5,
                                                  ecColor="orange")
                                    self.addArrow(self.ui.pl3, (child.mz + 1.00335 * child.xCount / child.loading, 0), (
                                    child.mz + 1.00335 * child.xCount / child.loading, intRight * intMul), linewidth=5,
                                                  ecColor="orange")

                                    intErrN, intErrL = self.getAllowedIsotopeRatioErrorsForResult()

                                    for iso in [1, 2, 3]:
                                        ratioN = getNormRatio(purN, child.xCount, iso)
                                        ratioL = getNormRatio(purL, child.xCount, iso)

                                        self.addArrow(self.ui.pl3, (child.mz + (1.00335 * iso) / child.loading, 0), (
                                        child.mz + (1.00335 * iso) / child.loading,
                                        intMul * intLeft * max(0., (ratioN - intErrN))), linewidth=5, alpha=.1, ecColor="orange")
                                        self.addArrow(self.ui.pl3, (child.mz + (1.00335 * iso) / child.loading,
                                                                    intMul * intLeft * max(0., (ratioN - intErrN))), (
                                                      child.mz + (1.00335 * iso) / child.loading,
                                                      intMul * intLeft * (ratioN + intErrN)), linewidth=5, alpha=.1,
                                                      ecColor="DarkSeaGreen")

                                        self.addArrow(self.ui.pl3,
                                                      (child.mz + 1.00335 * (child.xCount - iso) / child.loading, 0), (
                                                child.mz + 1.00335 * (child.xCount - iso) / child.loading,
                                                intMul * intRight * max(0., (ratioL - intErrL))), linewidth=5, alpha=.1,
                                                      ecColor="orange")
                                        self.addArrow(self.ui.pl3, (
                                        child.mz + 1.00335 * (child.xCount - iso) / child.loading,
                                        intMul * intRight * max(0., (ratioL - intErrL))), (
                                                          child.mz + 1.00335 * (child.xCount - iso) / child.loading,
                                                          intMul * intRight * (ratioL + intErrL)), linewidth=5, alpha=.1,
                                                      ecColor="DarkSeaGreen")

                                    if self.ui.drawFPIsotopologues.checkState() == QtCore.Qt.Checked:
                                        for iso in [1, 2, 3]:
                                            self.addArrow(self.ui.pl3, (child.mz - (1.00335 * iso) / child.loading, 0),
                                                          (child.mz - (1.00335 * iso) / child.loading,
                                                           intLeft * .1 * intMul), linewidth=5, ecColor="yellow")
                                            self.addArrow(self.ui.pl3, (
                                            child.mz + 1.00335 * (child.xCount + iso) / child.loading, 0), (
                                                              child.mz + 1.00335 * (child.xCount + iso) / child.loading,
                                                              intRight * .1 * intMul), linewidth=5, ecColor="yellow")

                                        self.addArrow(self.ui.pl3, (child.mz - 1.00335 / (child.loading * 2), 0), (
                                        child.mz - 1.00335 / (child.loading * 2), intLeft * .05 * intMul), linewidth=5,
                                                      ecColor="yellow")
                                        self.addArrow(self.ui.pl3, (child.mz + 1.00335 / (child.loading * 2), 0), (
                                        child.mz + 1.00335 / (child.loading * 2), intLeft * .05 * intMul), linewidth=5,
                                                      ecColor="yellow")
                                        self.addArrow(self.ui.pl3, (
                                            child.mz + 1.00335 * child.xCount / child.loading + 1.00335 / (
                                                child.loading * 2), 0), (
                                                          child.mz + 1.00335 * child.xCount / child.loading + 1.00335 / (
                                                              child.loading * 2), intRight * .05 * intMul), linewidth=5,
                                                      ecColor="yellow")
                                        self.addArrow(self.ui.pl3, (
                                            child.mz + 1.00335 * child.xCount / child.loading - 1.00335 / (
                                                child.loading * 2), 0), (
                                                          child.mz + 1.00335 * child.xCount / child.loading - 1.00335 / (
                                                              child.loading * 2), intRight * .05 * intMul), linewidth=5,
                                                      ecColor="yellow")

                            childIDs.append(child.id)

                            mzs.append(child.mz)
                            peaks.append(child.NPeakCenterMin / 60.)

                            xic = []
                            xicL = []
                            times = []

                            invert=1
                            if self.ui.negEIC.isChecked():
                                invert=-1

                            for row in self.currentOpenResultsFile.curs.execute("SELECT xic, xicL, xicfirstiso, xicLfirstiso, xicLfirstisoconjugate, xic_smoothed, xicL_smoothed, times FROM XICs WHERE id==%d" % child.eicID):
                                xic                   = [float(t) for t in row[0].split(";")]
                                xicL                  = [float(t) for t in row[1].split(";")]
                                xicfirstiso           = [float(t) for t in row[2].split(";")]
                                xicLfirstiso          = [float(t) for t in row[3].split(";")]
                                xicLfirstisoconjugate = [float(t) for t in row[4].split(";")]
                                xic_smoothed          = [float(t) for t in row[5].split(";")]
                                xicL_smoothed         = [float(t) for t in row[6].split(";")]
                                times = [float(t) / 60. for t in row[7].split(";")]

                            minTime = min(minTime, min(times[int(child.NPeakCenter - child.NPeakScale * 1):int(
                                child.NPeakCenter + child.NPeakScale * 1)]))
                            maxTime = max(maxTime, max(times[int(child.NPeakCenter - child.NPeakScale * 1):int(
                                child.NPeakCenter + child.NPeakScale * 1)]))

                            if self.ui.scaleFeatures.isChecked():
                                s = int(child.NPeakCenter - child.NBorderLeft * 1)
                                e = int(child.NPeakCenter + child.NBorderRight * 1)
                                maxP = s + max(range(e - s), key=lambda x: xic[s + x])
                                maxE  = mean([xic[maxP - 3], xic[maxP - 2], xic[maxP - 1], xic[maxP], xic[maxP + 1],xic[maxP + 2], xic[maxP + 3]])
                                maxEL = mean([xicL[maxP - 3], xicL[maxP - 2], xicL[maxP - 1], xicL[maxP], xicL[maxP + 1],xicL[maxP + 2], xicL[maxP + 3]])

                                if maxE == 0:
                                    maxE=1
                                if maxEL == 0:
                                    maxEL=1
                                if not self.ui.scaleLabelledFeatures.isChecked():
                                    maxEL=maxE

                                xic = [u / maxE for u in xic]
                                xic_smoothed = [u / maxE for u in xic_smoothed]
                                xicfirstiso = [u / maxE for u in xicfirstiso]

                                xicL = [u / maxEL for u in xicL]
                                xicL_smoothed = [u / maxEL for u in xicL_smoothed]
                                xicLfirstiso = [u / maxEL for u in xicLfirstiso]
                                xicLfirstisoconjugate = [u / maxEL for u in xicLfirstisoconjugate]

                            xicL = [invert*u for u in xicL]
                            xicL_smoothed = [invert*u for u in xicL_smoothed]
                            xicLfirstiso = [invert*u for u in xicLfirstiso]
                            xicLfirstisoconjugate = [invert*u for u in xicLfirstisoconjugate]

                            if self.ui.flattenXIC.isChecked():
                                ps = int(child.NPeakCenter - child.NPeakScale * 2)
                                pe = int(child.NPeakCenter + child.NPeakScale * 2)
                                for u in range(len(xic)):
                                    if u < ps or u > pe:
                                        xic[u] = 0

                                ps = int(child.LPeakCenter - child.LPeakScale * 2)
                                pe = int(child.LPeakCenter + child.LPeakScale * 2)
                                for u in range(len(xicL)):
                                    if u < ps or u > pe:
                                        xicL[u] = 0

                            maxInt = max(maxInt, max(xic[int(child.NPeakCenter - child.NPeakScale * 1): int(
                                child.NPeakCenter + child.NPeakScale * 1)]))

                            maxIntY = max(maxIntY, max(xic[int(child.NPeakCenter - child.NPeakScale * 1): int(
                                child.NPeakCenter + child.NPeakScale * 1)]))
                            minIntY = min(minIntY, min(xicL[int(child.LPeakCenter - child.LPeakScale * 1): int(
                                child.LPeakCenter + child.LPeakScale * 1)]))

                            self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xic, fill=[
                                max(int(child.LPeakCenter - child.LBorderLeft),
                                    int(child.NPeakCenter - child.NBorderLeft)),
                                min(int(child.NPeakCenter + child.NBorderRight),
                                    int(child.LPeakCenter + child.LBorderRight))], rearrange=len(selectedItems) == 1,
                                          label="%.4f (%d)"%(child.mz, child.xCount), useCol=useColi)  #useCol=selIndex*2)

                            if self.ui.showSmoothedEIC_checkBox.isChecked():
                                self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xic_smoothed, fill=[
                                    max(int(child.LPeakCenter - child.LBorderLeft),
                                        int(child.NPeakCenter - child.NBorderLeft)),
                                    min(int(child.NPeakCenter + child.NBorderRight),
                                        int(child.LPeakCenter + child.LBorderRight))], rearrange=len(selectedItems) == 1,
                                              label="", useCol=useColi, linestyle="--")

                            if self.ui.showIsotopologues.isChecked():
                                self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xicfirstiso, fill=[
                                    max(int(child.LPeakCenter - child.LBorderLeft),
                                        int(child.NPeakCenter - child.NBorderLeft)),
                                    min(int(child.NPeakCenter + child.NBorderRight),
                                        int(child.LPeakCenter + child.LBorderRight))],
                                              rearrange=len(selectedItems) == 1, label=None, useCol=useColi)

                            self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xicL, fill=[
                                max(int(child.LPeakCenter - child.LBorderLeft),
                                    int(child.NPeakCenter - child.NBorderLeft)),
                                min(int(child.NPeakCenter + child.NBorderRight),
                                    int(child.LPeakCenter + child.LBorderRight))], rearrange=len(selectedItems) == 1,
                                          label=None, useCol=useColi)

                            if self.ui.showSmoothedEIC_checkBox.isChecked():
                                self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xicL_smoothed, fill=[
                                    max(int(child.LPeakCenter - child.LBorderLeft),
                                        int(child.NPeakCenter - child.NBorderLeft)),
                                    min(int(child.NPeakCenter + child.NBorderRight),
                                        int(child.LPeakCenter + child.LBorderRight))], rearrange=len(selectedItems) == 1,
                                              label=None, useCol=useColi, linestyle="--")


                            if self.ui.showIsotopologues.isChecked():
                                self.drawPlot(self.ui.pl1, plotIndex=0, x=times, y=xicLfirstiso, fill=[
                                    max(int(child.LPeakCenter - child.LBorderLeft),
                                        int(child.NPeakCenter - child.NBorderLeft)),
                                    min(int(child.NPeakCenter + child.NBorderRight),
                                        int(child.LPeakCenter + child.LBorderRight))],
                                              rearrange=len(selectedItems) == 1,
                                              label=None, useCol=useColi)

                            maxindex, maxvalue = max(
                                enumerate(xic[int(child.NPeakCenter - 1):int(child.NPeakCenter + 1)],
                                          start=int(child.NPeakCenter - 1)), key=itemgetter(1))
                            meanRT += child.NPeakCenterMin / 60.
                            countK += 1

                    meanRT /= countK

                    if self.ui.plotAddLabels.checkState() == QtCore.Qt.Checked:
                        self.addAnnotation(self.ui.pl1, "%s\n@ %.2f" % (item.myData.featureName, meanRT), (meanRT, maxInt),
                                           (times[maxindex], xic[maxindex]), 0)
                useColi += 1

            if len(selFeatureGroups)>=1:
                childIDs=[]
                for selFeatureGroup in selFeatureGroups:
                    for childi in range(selFeatureGroup.childCount()):
                        child = selFeatureGroup.child(childi).myData
                        childIDs.append(child.id)
                if len(childIDs)>=1:
                    try:
                        self.clearPlot(self.ui.pl2)
                        plt.cla()
                        self.ui.pl2.fig.subplots_adjust(left=0.15, bottom=0.05, right=0.99, top=0.85)

                        fRows = 0
                        minCorr = 1
                        for row in self.currentOpenResultsFile.curs.execute("SELECT key, value FROM config WHERE key='minCorrelation'"):
                            fRows += 1
                            minCorr = float(row[1])
                        assert 0 < fRows <= 1, "Min Correlation not found or found multiple times in settings"

                        data = []
                        featureIDToColsNum = {}
                        texts = []
                        for i in range(len(childIDs)):
                            arow = []
                            for j in range(len(childIDs)):
                                arow.append(0)
                            data.append(arow)
                            texts.append("")
                            featureIDToColsNum[childIDs[i]] = i

                        cIds = ",".join(["%d" % f for f in childIDs])

                        for row in self.currentOpenResultsFile.curs.execute(
                                        "SELECT c.id, c.mz, c.xcount, c.ionMode FROM chromPeaks c WHERE c.id IN (%s)" % (cIds)):
                            id, mz, xcount, ionMode=row

                            fI1 = featureIDToColsNum[id]
                            texts[fI1]="%s%.4f/%d"%(ionMode, mz, xcount)

                        for row in self.currentOpenResultsFile.curs.execute(
                                        "SELECT ff.fID1, ff.fID2, ff.corr FROM featurefeatures ff WHERE ff.fID1 IN (%s) AND ff.fID2 IN (%s)" % (
                                        cIds, cIds)):
                            correlation = row[2]

                            fI1 = featureIDToColsNum[row[0]]
                            fI2 = featureIDToColsNum[row[1]]

                            if fI1 > fI2:
                                a = fI2
                                fI2 = fI1
                                fI1 = a

                            data[fI1][fI2] = correlation
                            data[fI2][fI1] = correlation

                        dv=[]
                        for i in range(len(data)):
                            data[i][i] = 1
                            dv.extend(data[i])

                        hc=HCA_general.HCA_generic()
                        tree=hc.generateTree(data)

                        datOrd=range(len(data))

                        datOrd=hc.getObjsOrderInTree(tree)
                        dnew=[]
                        for i in datOrd:
                            dnew.append([data[i][j] for j in datOrd])
                        data=dnew

                        data[len(data) - 1][len(data) - 1] = 1

                        colorDict = {'red': ((0.0, 0, convertXaToX(47 / 255., 0.2)),
                                             (0.5, convertXaToX(47 / 255., 0.5), convertXaToX(178 / 255., 0.225)), (
                        minCorr + (1 - minCorr) / 2., convertXaToX(178 / 255., 0.525),
                        convertXaToX(154 / 255., .5)), (1.0, convertXaToX(154 / 255., .9), 0)),

                                     'green': ((0.0, 0, convertXaToX(79 / 255., 0.2)),
                                               (0.5, convertXaToX(79 / 255., 0.5), convertXaToX(34 / 255., 0.225)),
                                               (minCorr + (1 - minCorr) / 2., convertXaToX(34 / 255., 0.525),
                                                convertXaToX(205 / 255., .5)),
                                               (1.0, convertXaToX(205 / 255., .9), 0)),

                                     'blue': ((0.0, 0, convertXaToX(79 / 255., 0.2)),
                                              (0.5, convertXaToX(79 / 255., 0.5), convertXaToX(34 / 255., 0.225)), (
                                     minCorr + (1 - minCorr) / 2., convertXaToX(34 / 255., 0.525),
                                     convertXaToX(50 / 255., .5)), (1.0, convertXaToX(50 / 255., 0.9), .0))}
                        plt.register_cmap(name='custom_colormap', data=colorDict)
                        cax = self.ui.pl2.twinxs[0].matshow(data, cmap=get_cmap("custom_colormap"))

                        self.ui.pl2.axes.set_xticks([i for i in range(len(data))])
                        self.ui.pl2.axes.set_xticklabels([texts[i] for i in datOrd], rotation=90)
                        self.ui.pl2.axes.set_yticks([i for i in range(len(data))])
                        self.ui.pl2.axes.set_yticklabels([texts[i] for i in datOrd])

                        self.ui.pl2.axes.set_yticks([i - .5 for i in range(len(data) + 1)], minor=True)
                        self.ui.pl2.axes.set_xticks([i - .5 for i in range(len(data) + 1)], minor=True)

                        self.ui.pl2.axes.grid(ls='solid', which='minor', color='white', linewidth=2)

                        cax.set_clim(-1, 1)
                        cb = self.ui.pl2.fig.colorbar(cax, ticks=[-1, 0, minCorr])
                        cb.outline.set_color('white')
                        cb.outline.set_linewidth(0)

                        self.ui.pl2.fig.canvas.draw()
                    except Exception as ex:
                        import traceback
                        traceback.print_exc()
                        logging.error(str(traceback))

            #</editor-fold>



            selIndex += 1

        #<editor-fold desc="#featureGroup results">

        if len(selectedItems)==1 and selectedItems[0].myType == "diagnostic - observed intensities":
            intensities=[]
            for row in self.currentOpenResultsFile.curs.execute("SELECT mz, intensity, lmz, tmz, xcount, loading FROM mzs"):
                mz, intensity, lmz, tmz, xcount, loading=row
                intensities.append(intensity)

            self.ui.pl1.twinxs[0].hist(intensities, 50, normed=False, facecolor='green', alpha=0.5)
            self.drawCanvas(self.ui.pl1)
        elif len(selectedItems)==1 and selectedItems[0].myType == "diagnostic - relative mz error":
            mzdifferrors=[]
            for row in self.currentOpenResultsFile.curs.execute("SELECT mz, lmz, tmz, xcount, loading FROM mzs"):
                mz, lmz, tmz, xcount, loading=row
                mzdifferrors.append((lmz-mz-tmz)*1000000/mz)

            self.ui.pl1.twinxs[0].hist(mzdifferrors, 50, normed=False, facecolor='green', alpha=0.5)
            self.drawCanvas(self.ui.pl1)
        elif len(selectedItems)==1 and selectedItems[0].myType == "diagnostic - relative mzbin deviation":
            mzdeviations=[]
            for row in self.currentOpenResultsFile.curs.execute("SELECT mzbinid, min(mzs.mz), max(mzs.mz), count(mzs.mz) FROM mzbinskids INNER JOIN mzs ON mzs.id=mzbinskids.mzid GROUP BY mzbinid"):
                binid, mzmin, mzmax, n=row
                if n>1:
                    mzdeviations.append((mzmax-mzmin)*1000000/mzmin)

            self.ui.pl1.twinxs[0].hist(mzdeviations, 50, normed=False, facecolor='green', alpha=0.5)
            self.drawCanvas(self.ui.pl1)
        elif len(selectedItems)==1 and selectedItems[0].myType == "diagnostic - feature pair correlations":
            peaksCorr=[]
            for row in self.currentOpenResultsFile.curs.execute("SELECT peaksCorr FROM chromPeaks"):
                peakCorr,=row
                peaksCorr.append(peakCorr)

            self.ui.pl1.twinxs[0].hist(peaksCorr, 10, normed=False, facecolor='green', alpha=0.5)
            self.drawCanvas(self.ui.pl1)
        elif len(selectedItems)==1 and selectedItems[0].myType == "diagnostic - feature pair mp1 to m ratio":
            peaksRatio=[]
            for row in self.currentOpenResultsFile.curs.execute("SELECT peaksRatioMp1 FROM chromPeaks"):
                peakRatio,=row
                peaksRatio.append(peakRatio)

            self.ui.pl1.twinxs[0].hist(peaksRatio, 30, normed=False, facecolor='green', alpha=0.5)
            self.drawCanvas(self.ui.pl1)
        elif len(selectedItems)==1 and selectedItems[0].myType == "diagnostic - feature pair mPp1 to mP ratio":
            peaksRatio=[]
            for row in self.currentOpenResultsFile.curs.execute("SELECT peaksRatioMPm1 FROM chromPeaks"):
                peakRatio,=row
                peaksRatio.append(peakRatio)

            self.ui.pl1.twinxs[0].hist(peaksRatio, 30, normed=False, facecolor='green', alpha=0.5)
            self.drawCanvas(self.ui.pl1)
        elif len(selectedItems)==1 and selectedItems[0].myType == "diagnostic - feature pair assigned mzs":
            assignedmzs=[]
            for row in self.currentOpenResultsFile.curs.execute("SELECT assignedMZs FROM chromPeaks"):
                n,=row
                assignedmzs.append(n)

            self.ui.pl1.twinxs[0].hist(assignedmzs, 30, normed=False, facecolor='green', alpha=0.5)
            self.drawCanvas(self.ui.pl1)
        elif len(selectedItems)==1 and selectedItems[0].myType == "diagnostic - feature pair mz deviation":
            devMeans=[]
            devVals=[]
            for row in self.currentOpenResultsFile.curs.execute("SELECT mzdifferrors FROM chromPeaks"):
                mzdifferrors,=row
                mzdifferrors=loads(base64.b64decode(mzdifferrors))
                devMeans.append([mzdifferrors.mean if mzdifferrors.mean is not None else -1])
                devVals.extend([v for v in mzdifferrors.vals if v is not None])

            self.ui.pl1.twinxs[0].hist(devMeans, 30, normed=False, facecolor='green', alpha=0.5, label="Means")
            self.ui.pl1.twinxs[0].hist(devVals, 30, normed=False, facecolor='blue', alpha=0.5, label="Scan-level")
            self.ui.pl1.twinxs[0].legend(loc='upper right')
            self.drawCanvas(self.ui.pl1)


        #</editor-fold>

        #<editor-fold desc="#feature pair plotting">
        elif "Features" in plotTypes or "Feature Groups" in plotTypes or "feature" in plotTypes:
            if self.ui.scaleFeatures.checkState() == QtCore.Qt.Checked:
                self.ui.pl1.axes.set_ylabel("Intensity [counts; normalised]")
            else:
                self.ui.pl1.axes.set_ylabel("Intensity [counts]")
            if self.ui.autoZoomPlot.checkState() == QtCore.Qt.Checked:
                self.drawCanvas(self.ui.pl1, ylim=(minIntY * 1.6, maxIntY * 1.6), xlim=(minTime * 0.85, maxTime * 1.15))
            else:
                self.drawCanvas(self.ui.pl1)
        #</editor-fold>

        #<editor-fold desc="#mz and mzbin plotting">
        elif "MZs" in plotTypes or "mz" in plotTypes or "mzbin" in plotTypes or "MZBins" in plotTypes:
            colour_LUT = ['#0000FF', '#00FF00', '#FF00FF', '#00FFFF', '#FFFF00', '#FFFFFF', '#F0F0F0', '#0F0F0F']

            # the scatter plot:
            maxIntX = max([max(x) for x in x_vals])
            maxIntY = max([max(y) for y in y_vals])
            colors = []
            axScatter = self.ui.pl1.axes

            for i in range(self.ui.res_ExtractedData.topLevelItem(2).childCount()):
                child = self.ui.res_ExtractedData.topLevelItem(2).child(i)
                assert child.myType == "feature"

                tppm = float(self.getParametersFromCurrentRes("Mass deviation (+/- ppm)"))
                ppmupper = (1 + tppm / 1000000.)
                ppmlower = (1 - tppm / 1000000.)

                axScatter.plot([child.myData.NPeakCenterMin / 60. - child.myData.NPeakScale / 60.,
                                child.myData.NPeakCenterMin / 60. + child.myData.NPeakScale / 60.,
                                child.myData.NPeakCenterMin / 60. + child.myData.NPeakScale / 60.,
                                child.myData.NPeakCenterMin / 60. - child.myData.NPeakScale / 60.,
                                child.myData.NPeakCenterMin / 60. - child.myData.NPeakScale / 60.],
                               [child.myData.mz * ppmlower, child.myData.mz * ppmlower, child.myData.mz * ppmupper,
                                child.myData.mz * ppmupper, child.myData.mz * ppmlower], color="red")

                axScatter.plot([child.myData.NPeakCenterMin / 60. - child.myData.NPeakScale / 60.,
                                child.myData.NPeakCenterMin / 60. + child.myData.NPeakScale / 60.],
                               [child.myData.mz, child.myData.mz], color="red")

            for i in range(len(x_vals)):
                colour = colour_LUT[i % len(colour_LUT)]
                axScatter.scatter(x_vals[i], y_vals[i], color=colour, linewidth=0)
                colors.append(colour)

            self.ui.pl1.twinxs[0].set_ylabel("m/z")
            self.ui.pl1.twinxs[0].set_xlabel("Retention time [minutes]")
            self.drawCanvas(self.ui.pl1, ylim=[0, maxIntY], xlim=[0, maxIntX])
        #</editor-fold>

        else:
            self.drawCanvas(self.ui.pl1)

        self.ui.pl3.fig.canvas.draw()
    #</editor-fold>

    #<editor-fold desc="### general plotting functions">
    def addArrow(self, plt, point, at, plotIndex=0, fcColor="white", ecColor="white", arrowColor="slategrey", alpha=1,
                 arrowAlpha=1, drawArrowHead=False, linewidth=1):

        if point != at:
            if not drawArrowHead:
                plt.twinxs[plotIndex].arrow(at[0], at[1], point[0] - at[0], point[1] - at[1], color=ecColor,
                                            shape="right", head_width=(point[0] - at[0]) * 0.3,
                                            head_length=(point[1] - at[1]) * 0., linewidth=linewidth)
            else:
                plt.twinxs[plotIndex].annotate("", xy=point, xytext=at, xycoords='data', textcoords='data', va="center",
                                               ha="center",
                                               bbox=dict(boxstyle='round', fc=fcColor, ec=ecColor, alpha=alpha),
                                               rotation=0, arrowprops=dict(arrowstyle="->", connectionstyle='arc3',
                                                                           color=arrowColor, alpha=arrowAlpha))

    def addCircle(self, plt, at, plotIndex):
        pass

    def addAnnotation(self, plt, text, point, at, plotIndex=0, rotation=0, arrowAlpha=0.5, offset=(-10, 80),
                      fcColor="white", ecColor="white", arrowColor="firebrick", alpha=0.8, up=True, add=80):
        if not up:
            add = -add
        plt.twinxs[plotIndex].annotate(text, xy=point, xytext=(-10, add), xycoords='data', textcoords='offset points',
                                       va="center", ha="center",
                                       bbox=dict(boxstyle='round', fc=fcColor, ec=ecColor, alpha=alpha),
                                       rotation=rotation,
                                       arrowprops=dict(arrowstyle='wedge', color=arrowColor, alpha=arrowAlpha))

    def clearPlot(self, plt, setXtoZero=False):
        for ax in plt.twinxs:
            plt.fig.delaxes(ax)
        plt.axes = plt.fig.add_subplot(111)

        simpleaxis(plt.axes)
        if setXtoZero:
            plt.axes.spines['bottom'].set_position('zero')

        plt.axes.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.twinxs = [plt.axes]

    def drawPlot(self, plt, plotIndex=0, x=range(10), y=range(1, 11), fill=[], label="", b=None, rearrange=True,
                 useCol=None, multipleLocator=5, alpha=globAlpha, title="", xlab="Retention time [minutes]",
                 ylab="Intensity [counts]", plot=True, scatter=False, linestyle="-"):
        try:
            if b is None:
                b = predefinedColors

            if useCol is None:
                useCol = plotIndex

            ylim = None
            xlim = None
            if len(plt.twinxs) > 0 and not rearrange:
                ylim = plt.twinxs[0].get_ylim()
                xlim = plt.twinxs[0].get_xlim()

            if plotIndex == len(plt.twinxs):
                plt.twinxs.append(plt.axes.twinx())
            if plotIndex == 0:
                pass

            ax = plt.twinxs[plotIndex]

            if multipleLocator is not None:
                sf = ScalarFormatter(useOffset=True, useMathText=True)
                sf.set_scientific(True)
                sf.set_powerlimits((10, 0))
                ax.yaxis.set_major_formatter(sf)
                #ax.xaxis.set_major_locator(MultipleLocator(multipleLocator))
                #ax.xaxis.set_minor_locator(MultipleLocator(0.25))

            ax.set_title(title)
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab)

            plotCol=useCol
            if isinstance(useCol, int):
                plotCol=b[useCol % len(b)]

            if plot:
                ax.plot(x, y, color=plotCol, label=label, linestyle=linestyle)
            if scatter:
                ax.scatter(x, y, color=plotCol)
            if len(fill) > 0 and self.ui.plotMarkArea.checkState() == QtCore.Qt.Checked:
                ax.fill_between(x[fill[0]:fill[1]], y[fill[0]:fill[1]], color=plotCol, alpha=alpha)

            #if ylim is not None:
            #    for ax in plt.twinxs:
            #        ax.set_ylim(ylim)
            #        ax.set_xlim(xlim)
        except Exception as ex:
            logging.warning("Exception caught", ex)

    def setLimts(self, plt, ylim, xlim):
        for ax in plt.twinxs:
            ax.set_ylim(ylim[0], ylim[1])
            ax.set_xlim(xlim[0], xlim[1])

    def drawCanvas(self, plt, ylim=None, xlim=None):
        #if ylim is None:
        #    ylim = plt.twinxs[0].get_ylim()
        #if xlim is None:
        #    xlim = plt.twinxs[0].get_xlim()

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
            if self.ui.showLegend.isChecked():
                ax.legend()

        plt.canvas.draw()
    #</editor-fold>


    def multipleFilesPlotAdd(self):

        if self.ui.pl2.pictureShown:
            self.clearPlot(self.ui.pl2)
            self.ui.pl2.pictureShown = False

        maxInt = 0
        mzs = []
        peaks = []
        for item in self.ui.res_ExtractedData.selectedItems():
            if item.myType == "feature":
                if self.ui.pl2.type is None or self.ui.pl2.type == "feature":
                    self.ui.pl2.type = "feature"
                else:
                    QtGui.QMessageBox.warning(self, "MetExtract",
                                              "Selecting different result types in not supported. Please select only mz, mzbins or chromatographic peaks at a time",
                                              QtGui.QMessageBox.Ok)
                    continue

                mzs.append(item.myData.id)
                peaks.append(item.myData.NPeakCenterMin / 60.)
                xic = []
                times = []
                for row in self.currentOpenResultsFile.curs.execute(
                                "select xic, times, xicL from XICs where id==%d" % item.myData.eicID):
                    xic = [float(t) for t in row[0].split(";")]
                    times = [float(t) / 60. for t in row[1].split(";")]
                    xicL = [float(t) for t in row[2].split(";")]
                self.ui.pl2.xics.append(xic)
                self.ui.pl2.times.append(times)
                self.ui.pl2.peaks.append([item.myData])
            if item.myType == "MZs":
                if self.ui.pl2.type is None or self.ui.pl2.type == "MZs":
                    self.ui.pl2.type = "MZs"
                else:
                    QtGui.QMessageBox.warning(self, "MetExtract",
                                              "Selecting different result types is not supported. Please select only mz, mzbins or chromatographic peaks at a time",
                                              QtGui.QMessageBox.Ok)
                    continue
                t = item

                x = []
                y = []
                for j in range(t.childCount()):
                    child = t.child(j)
                    assert child.myType == "mz"
                    x.append(child.myData[3] / 60.)
                    y.append(child.myData[0])

                self.ui.pl2.x_vals.append(x)
                self.ui.pl2.y_vals.append(y)

        if self.ui.pl2.type == "feature":
            for i in range(len(self.ui.pl2.peaks)):
                maxInt = max(maxInt, max(self.ui.pl2.xics[i]))
                maxindex, maxvalue = max(enumerate(
                    xic[int(self.ui.pl2.peaks[i][0].NPeakCenter - 3):int(self.ui.pl2.peaks[i][0].NPeakCenter + 3)],
                    start=int(self.ui.pl2.peaks[i][0].NPeakCenter - 3)), key=itemgetter(1))
                if i == (len(self.ui.pl2.peaks) - 1):
                    self.drawPlot(self.ui.pl2, plotIndex=0, x=self.ui.pl2.times[i], y=self.ui.pl2.xics[i], fill=[
                        int(self.ui.pl2.peaks[i][0].NPeakCenter - self.ui.pl2.peaks[i][0].NBorderLeft),
                        int(self.ui.pl2.peaks[i][0].NPeakCenter + self.ui.pl2.peaks[i][0].NBorderRight)],
                                  useCol=i * 2)
                    if self.ui.plotAddLabels.checkState() == QtCore.Qt.Checked:
                        self.addAnnotation(self.ui.pl2, "%s\n%.5f\n@ %.2f" % (
                            self.ui.pl2.peaks[i][0].assignedName, self.ui.pl2.peaks[i][0].mz,
                            self.ui.pl2.peaks[i][0].NPeakCenterMin / 60.),
                                           (self.ui.pl2.times[i][maxindex], self.ui.pl2.xics[i][maxindex]),
                                           (self.ui.pl2.times[i][maxindex], self.ui.pl2.xics[i][maxindex]), 0)
            self.drawCanvas(self.ui.pl2, ylim=[0, maxInt * 1.2])
        if self.ui.pl2.type == "MZs":

            #matplotlib code taken from http://stackoverflow.com/questions/6508769/matplotlib-scatter-hist-with-stepfilled-histtype-in-histogram
            colour_LUT = ['#0000FF', '#FF0000', '#00FF00', '#FF00FF', '#00FFFF', '#FFFF00', '#FFFFFF']

            # the scatter plot:
            maxIntX = max([max(x) for x in self.ui.pl2.x_vals])
            maxIntY = max([max(y) for y in self.ui.pl2.y_vals])
            colors = []
            axScatter = self.ui.pl2.axes

            for i in range(len(self.ui.pl2.x_vals)):
                colour = colour_LUT[i % len(colour_LUT)]
                axScatter.scatter(self.ui.pl2.x_vals[i], self.ui.pl2.y_vals[i], color=colour)
                colors.append(colour)
            self.drawCanvas(self.ui.pl2, ylim=[0, maxIntY], xlim=[0, maxIntX])

    def res_doubleClick(self, item):
        self.multipleFilesPlotAdd()

    def resChildNode(self, node):
        node.setHidden(False)
        for c in range(node.childCount()):
            self.resChildNode(node.child(c))

    def resetFilter(self):
        for i in range(self.ui.res_ExtractedData.topLevelItemCount()):
            d = self.ui.res_ExtractedData.topLevelItem(i)
            self.resChildNode(d)

    def filterHelp1_matchMZ(self, ref, akt, maxDiff, type="ppm"):
        ref = str(ref)
        if " " in ref:
            ref = ref[0:ref.find(" ")]
        if "(" in ref:
            ref = ref[0:ref.find("(")]
        ref = float(str(ref))
        akt = float(str(akt))
        maxDiff = float(str(maxDiff))
        if type == "ppm":
            if (abs(akt - ref) * 1000000. / akt) <= maxDiff:
                return True
        if type == "amu":
            if abs(akt - ref) <= maxDiff:
                return True
        return False

    def filterEdited(self, text):
        if self.currentOpenResultsFile is None:
            return

        text = str(text)
        self.resetFilter()

        if len(text) > 0:
            text = text.split(" ")
            textl = []
            for w in range(len(text)):
                textl.append(lambda x: x.contains(text[w]))

            if "+-" in text[0]:
                h = text[0][0:text[0].rfind("+")]
                textl[0] = lambda x: self.filterHelp1_matchMZ(x, h, 5.)
                if "ppm" in text[0]:
                    ppm = text[0][(text[0].rfind("-") + 1):text[0].find("p")]
                    textl[0] = lambda x: self.filterHelp1_matchMZ(x, h, ppm, "ppm")
                elif "a" in text[0]:
                    amu = text[0][(text[0].rfind("-") + 1):text[0].find("a")]
                    textl[0] = lambda x: self.filterHelp1_matchMZ(x, h, amu, "amu")
                else:
                    amu = text[0][(text[0].rfind("-") + 1):]
                    if len(amu) > 0:
                        textl[0] = lambda x: self.filterHelp1_matchMZ(x, h, amu, "amu")
                    else:
                        textl[0] = lambda x: self.filterHelp1_matchMZ(x, h, 1, "amu")


            #process first 3 result categories (mzs, bins, features)
            for i in range(3):
                d = self.ui.res_ExtractedData.topLevelItem(i)
                for c in range(d.childCount()):
                    dd = d.child(c)
                    show = True
                    for w in range(len(text)):
                        if len(text[w]) > 0 and not (textl[w](dd.text(w))):  #not(dd.text(w).contains(text[w])):
                            show = False
                    dd.setHidden(not show)

            #process feature groups (perform search on features within groups)   
            d = self.ui.res_ExtractedData.topLevelItem(3)
            for v in range(d.childCount()):
                allShow = False
                fg = d.child(v)
                for c in range(fg.childCount()):
                    dd = fg.child(c)
                    show = True
                    for w in range(len(text)):
                        if len(text[w]) > 0 and not (textl[w](dd.text(w))):  #not(dd.text(w).contains(text[w])):
                            show = False
                    dd.setHidden(not show)
                    allShow = allShow or show
                fg.setHidden(not allShow)

    def setChromPeakName(self):
        cpName = str(self.ui.chromPeakName.text())

        selectedItems = self.ui.res_ExtractedData.selectedItems()

        for item in selectedItems:
            if item.myType == "Features" or item.myType == "feature":
                self.currentOpenResultsFile.curs.execute(
                    "UPDATE chrompeaks SET assignedName='%s' WHERE id=%d" % (cpName, item.myID))
                self.currentOpenResultsFile.conn.commit()
                item.setText(0, cpName)
            if item.myType == "FeatureGroup" or item.myType == "featureGroup":
                self.currentOpenResultsFile.curs.execute(
                    "UPDATE featureGroups SET featureName='%s' WHERE id=%d" % (cpName, item.myID))
                self.currentOpenResultsFile.conn.commit()
                item.myData = cpName
                item.setText(0, cpName)

    def openRawFile(self):
        if self.currentOpenResultsFile.file is not None:
            os.startfile(self.currentOpenResultsFile.file)

    def removeChromPeaksFrom(self, item, idToRem):
        rem = 0
        todel = set()

        for j in range(item.childCount()):
            child = item.child(j)
            if child.myType == "feature" and child.myID == idToRem:
                todel.add(j)
        todel = list(todel)
        for j in sorted(todel, reverse=True):
            item.removeChild(item.child(j))
            rem += 1

        return rem

    def sortTreeChildren(self, node, sortColumn):
        children = []
        for ci in range(node.childCount() - 1, -1, -1):
            child = node.child(ci)
            children.append(child)
            node.removeChild(child)
        children = natSort(children, key=lambda x: str(x.text(sortColumn)))
        for child in children:
            node.addChild(child)

    def delNodeFromResults(self, item):
        chromPeaksRem = 0
        featureGroupsRem = 0
        if item.myType == "feature":
            cp = item.myData
            myType = item.myType
            myID = item.myID
            self.currentOpenResultsFile.curs.execute("delete from ChromPeaks where id=%d" % cp.id)
            self.currentOpenResultsFile.curs.execute("delete from featureGroupFeatures where fID=%d" % cp.id)
            self.currentOpenResultsFile.curs.execute(
                "delete from featureGroups where id not in (select distinct fGroupID from featureGroupFeatures)")
            self.ui.res_ExtractedData.setItemSelected(item, False)

            chromPeaksRem += self.removeChromPeaksFrom(self.ui.res_ExtractedData.topLevelItem(2), myID)

            todel = set()
            for i in range(self.ui.res_ExtractedData.topLevelItem(3).childCount()):
                featureGroup = self.ui.res_ExtractedData.topLevelItem(3).child(i)
                assert featureGroup.myType == "featureGroup"
                self.removeChromPeaksFrom(featureGroup, myID)
                featureGroup.setText(1, "%d" % featureGroup.childCount())
                if featureGroup.childCount() == 0:
                    todel.add(i)
            todel = list(todel)
            for i in sorted(todel, reverse=True):
                self.ui.res_ExtractedData.topLevelItem(3).removeChild(
                    self.ui.res_ExtractedData.topLevelItem(3).child(i))
                featureGroupsRem += 1
        elif item.myType == "featureGroup":
            myType = item.myType
            myID = item.myID
            partChromPeaks = []
            for j in range(item.childCount()):
                partChromPeaks.append(item.child(j))
            for child in partChromPeaks:
                z, u = self.delNodeFromResults(child)
                chromPeaksRem += z
                featureGroupsRem = featureGroupsRem + u

        self.ui.res_ExtractedData.topLevelItem(2).setText(1, "%d" % self.ui.res_ExtractedData.topLevelItem(2).childCount())
        self.ui.res_ExtractedData.topLevelItem(3).setText(1, "%d" % self.ui.res_ExtractedData.topLevelItem(3).childCount())

        return chromPeaksRem, featureGroupsRem

    def createGroupFrom(self, items, tracerID):
        fIDs = []
        for item in items:
            fIDs.append(item.myID)

        i = 0
        while True:
            i += 1
            found = 0
            for row in self.currentOpenResultsFile.curs.execute(
                            "select * from featureGroups where featureName='fg_%d'" % i):
                found += 1
            if found == 0:
                break

        SQLInsert(self.currentOpenResultsFile.curs, "featureGroups", featureName="fgU_%d"%i, tracer=tracerID)
        found = 0
        newID = -1
        newName = ""
        tracerName = ""
        for row in self.currentOpenResultsFile.curs.execute(
                        "select featureGroups.id, featureName, name from featureGroups inner join tracerConfiguration on featureGroups.tracer=tracerConfiguration.id where featureName='fgU_%d' " % i):
            found += 1
            newID = int(row[0])
            newName = str(row[1])
            tracerName = str(row[2])

        d = QtGui.QTreeWidgetItem([str(newName), str(len(items)), "", str(tracerName)])
        d.myType = "featureGroup"
        d.myID = newID
        d.myData = Bunch(fgID=newID, featureName=newName, tracerName=str(tracerName) )
        self.ui.res_ExtractedData.topLevelItem(3).addChild(d)

        self.currentOpenResultsFile.curs.execute("update featureGroupFeatures set fGroupID=%d where fID in (%s)" % (
            newID, ", ".join([str(fId) for fId in fIDs])))
        sumRt = 0
        cpCount = 0
        for item in items:
            parent = item.parent()
            #if len(d.myData) == 3:
            #    d.myData.append(parent.myData[3])
            parent.removeChild(item)
            parent.setText(1, "%d" % parent.childCount())
            d.addChild(item)
            xp = item.myData
            sumRt = sumRt + xp.NPeakCenterMin
            cpCount += 1

        d.setText(1, "%d" % d.childCount())
        d.setText(2, "%.2f" % (sumRt / cpCount / 60.))
        self.ui.res_ExtractedData.topLevelItem(3).setText(1, "%d" % self.ui.res_ExtractedData.topLevelItem(3).childCount())

        self.currentOpenResultsFile.curs.execute(
            "delete from featureGroups where id not in (select distinct fGroupID from featureGroupFeatures)")

        self.currentOpenResultsFile.conn.commit()

        for ci in range(self.ui.res_ExtractedData.topLevelItem(3).childCount() - 1, -1, -1):
            child = self.ui.res_ExtractedData.topLevelItem(3).child(ci)
            if child.childCount() == 0:
                self.ui.res_ExtractedData.topLevelItem(3).removeChild(child)

        self.sortTreeChildren(self.ui.res_ExtractedData.topLevelItem(3), 2)

        return newName

    def showPopup(self, position):
        selectedItems = self.ui.res_ExtractedData.selectedItems()

        types = set()
        parentTypes = set()
        parentItems = []

        for item in selectedItems:
            if hasattr(item, "myType"):
                types.add(item.myType)
                if item.parent() is not None:
                    parentItems.append(item.parent())
                    if hasattr(item, "myType"):
                        parentTypes.add(item.parent().myType)
            else:
                types.add("generic item")
        types = list(types)

        menu = QtGui.QMenu()
        actionAvailable = False
        tracerActions = []

        newGroupAction = -1
        if all([si == "feature" for si in types]) and all([pt == "featureGroup" for pt in parentTypes]):
            newGroupAction = menu.addMenu("Extract to new Group")
            actionAvailable = True
            for tracer in self.currentOpenResultsFile.curs.execute("select id, name from tracerConfiguration"):
                tracerAction = newGroupAction.addAction(str(tracer[1]))
                tracerAction.tracerID = int(tracer[0])
                tracerActions.append(tracerAction)

        deleteAction = -1
        menu.addSeparator()

        if len(types) == 1 and (types[0] == "feature" or types[0] == "featureGroup"):
            deleteAction = menu.addAction("Delete")
            actionAvailable = True

        clipboardAction = -1
        menu.addSeparator()
        if len(types) == 1 and (types[0] == "feature" or types[0] == "featureGroup" or types[0] == "Features" or types[
            0] == "Feature Groups"):
            clipboardAction = menu.addAction("Copy")
            actionAvailable = True

        if actionAvailable:
            action = menu.exec_(self.ui.res_ExtractedData.mapToGlobal(position))

            if action == clipboardAction:
                clipboard = ["MZ\tRT\tXn\tZ\tIonMode\tFG\tTracer\tAdduct\tHeteroAtoms"]

                if len(selectedItems) == 1 and types[0] == "Features":
                    t = []
                    for j in range(selectedItems[0].childCount()):
                        child = item.child(j)
                        t.append(child)
                    selectedItems = t
                    types = ["feature"]

                if len(selectedItems) == 1 and types[0] == "Feature Groups":
                    t = []
                    for j in range(selectedItems[0].childCount()):
                        child = item.child(j)
                        t.append(child)
                    selectedItems = t
                    types = ["featureGroup"]

                if types[0] == "feature":
                    for item in selectedItems:
                        clipboard.append("%f\t%.2f\t%d\t%d\t%s\t%s\t%s\t%s\t%s" % (
                            item.myData.mz, item.myData.NPeakCenterMin / 60., item.myData.xCount, item.myData.loading,
                            item.myData.ionMode, "-", item.myData.tracer, item.myData.adducts, item.myData.heteroAtoms))
                if types[0] == "featureGroup":
                    for item in selectedItems:
                        for j in range(item.childCount()):
                            child = item.child(j)
                            clipboard.append("%f\t%.2f\t%d\t%d\t%s\t%s\t%s\t%s\t%s" % (
                                child.myData.mz, child.myData.NPeakCenterMin / 60., child.myData.xCount, child.myData.loading,
                                child.myData.ionMode, item.myData.fgID, child.myData.tracer, child.myData.adducts, child.myData.heteroAtoms))

                pyperclip.copy("\n".join(clipboard))

            elif action == deleteAction:
                if QtGui.QMessageBox.question(self, "MetExtract",
                                              "Are you sure you want to delete the selected result(s)?\nThis action cannot be undone",
                                              QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.Yes:

                    pw = ProgressWrapper(1)
                    pw.getCallingFunction()("max")(0)
                    pw.getCallingFunction()("header")("working..")
                    pw.show()

                    featureGroupsRem = 0
                    chromPeaksRem = 0
                    for item in selectedItems:
                        z, u = self.delNodeFromResults(item)
                        chromPeaksRem += z
                        featureGroupsRem += u

                    self.currentOpenResultsFile.conn.commit()

                    pw.hide()

                    QtGui.QMessageBox.information(self, "MetExtract",
                                                  "%d features and %d feature groups have been deleted" % (
                                                      chromPeaksRem, featureGroupsRem), QtGui.QMessageBox.Ok)

            elif action in tracerActions:
                if QtGui.QMessageBox.question(self, "MetExtract",
                                              "Are you sure you want to extract the selected features in a new group?\nThis action cannot be undone",
                                              QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.Yes:

                    pw = ProgressWrapper(1)
                    pw.getCallingFunction()("max")(0)
                    pw.getCallingFunction()("header")("working..")
                    pw.show()

                    for item in selectedItems:
                        assert item.myType == "feature"
                    groupName = self.createGroupFrom(selectedItems, action.tracerID)

                    pw.hide()

                    QtGui.QMessageBox.information(self, "MetExtract",
                                                  "Created group %s with %d features" % (groupName, len(selectedItems)),
                                                  QtGui.QMessageBox.Ok)

    def heteroAtomsConfiguration(self):
        t = heteroAtomEdit(heteroAtoms=deepcopy(self.heteroElements))
        if t.executeDialog() == QtGui.QDialog.Accepted:
            self.heteroElements = t.getHeteroAtoms()

            logging.info("Configured heteroatoms:")
            for ha in self.heteroElements:
                logging.info("  * " + ha.name)

    def relationShipConfiguration(self):
        t = adductsEdit(adds=deepcopy(self.adducts), nls=deepcopy(self.elementsForNL))
        if t.executeDialog() == QtGui.QDialog.Accepted:
            self.adducts = t.getAdducts()

            logging.info("Configured adducts")
            for add in self.adducts:
                logging.info("  *", add.name)

            self.elementsForNL = t.getNeutralLosses()

            logging.info("Configured neutral loss elements")
            for elem in self.elementsForNL:
                logging.info("  *", elem.name)

    def updateCores(self):
        curMax = self.ui.cpuCores.maximum()
        curVal = self.ui.cpuCores.value()

        cpus = cpu_count()
        if self.ui.workingCore.isChecked():
            cpus -= 1
        self.ui.cpuCores.setMaximum(cpus)
        if curMax == curVal:
            self.ui.cpuCores.setValue(cpus)

    #<editor-fold desc="### Tracer Dialog">
    def setModuleUI(self, val):
        if val==TRACER:
            self.ui.label_50.setText("TracExtract experiment")

            self.ui.groupBox_ISOA.setVisible(False)
            self.ui.groupBox_ISOB.setVisible(False)
            self.ui.useCValidation.setVisible(False)

            self.ui.useRatio.setVisible(False)
            self.ui.ratioGroupBox.setVisible(False)

            self.ui.setupTracers.setVisible(True)
            self.ui.tracerExperimentLabel.setVisible(True)
        elif val==METABOLOME:
            self.ui.label_50.setText("AllExtract experiment")

            self.ui.setupTracers.setVisible(False)
            self.ui.tracerExperimentLabel.setVisible(False)

            self.ui.groupBox_ISOA.setVisible(True)
            self.ui.groupBox_ISOB.setVisible(True)
            self.ui.useCValidation.setVisible(True)

            self.ui.useRatio.setVisible(True)
            self.ui.ratioGroupBox.setVisible(True)
        else:
            raise Exception("undefined module provided")

    def showTracerEditor(self):
        tracerDialog = tracerEdit()
        tracerDialog.setTracers(deepcopy(self.configuredTracers))
        if tracerDialog.executeDialog() == QtGui.QDialog.Accepted:
            self.lastOpenDir = tracerDialog.getOpenDir()

            self.configuredTracers = tracerDialog.getTracers()

            logging.info("Configured tracers:")
            for tracer in self.configuredTracers:
                logging.info(" * %s (%s/%s)" % (tracer.name, tracer.isotopeA, tracer.isotopeB))
        self.updateTracerInfo()

    def updateTracerInfo(self):
        if len(self.configuredTracers) == 0:
            self.ui.tracerExperimentLabel.setText("No tracers configured")
        else:
            trcs = []
            for tracer in self.configuredTracers:
                trcs.append(
                    " * %s (%s/%s)" % (tracer.name, tracer.isotopeA, tracer.isotopeB))
            self.ui.tracerExperimentLabel.setText("\n\n".join(trcs))
    #</editor-fold>

    def isotopeATextChanged(self, text):
        text = str(text)
        self.isotopeAmass = getIsotopeMass(text)[0]
        if self.isotopeAmass == -1:
            self.ui.isotopeAMassLabel.setText("Unknown Isotope %s" % text)
        else:
            self.ui.isotopeAMassLabel.setText("%s mass is %.5f" % (text, self.isotopeAmass))

    def isotopeBTextChanged(self, text):
        text = str(text)
        self.isotopeBmass = getIsotopeMass(text)[0]
        if self.isotopeBmass == -1:
            self.ui.isotopeBMassLabel.setText("Unknown Isotope %s" % text)
        else:
            self.ui.isotopeBMassLabel.setText("%s mass is %.5f" % (text, self.isotopeBmass))

    def useRatioValueStateChanged(self):
        self.ui.ratioGroupBox.setEnabled(self.ui.useRatio.checkState()==QtCore.Qt.Checked)

    def alignToggled(self, val):
        self.ui.polynomLabel.setEnabled(val)
        self.ui.polynomValue.setEnabled(val)


    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        if event.mimeData().hasUrls:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()

            links = []
            for url in event.mimeData().urls():
                links.append(str(url.toLocalFile()).replace("\\", "/"))

            if len(links) == 1 and (links[0].lower().endswith(".grp") or links[0].lower().endswith(".ini")):
                link = links[0]
                if link.lower().endswith(".grp") and \
                   QtGui.QMessageBox.question(self, "MetExtract", "Do you want to load the groups defined in this file?",
                                      QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.Yes:
                    self.loadGroups(link, askLoadSettings=False)
                if (link.lower().endswith(".ini") or link.lower().endswith(".grp")) and \
                   QtGui.QMessageBox.question(self, "MetExtract", "Do you want to load the settings saved in this file?",
                                      QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.Yes:
                    self.loadSettingsFile(link)
            else:
                self.showAddGroupDialog(initWithFiles=links)
        else:
            event.ignore()

    def loadAvailableSettingsFile(self, file):
        try:
            self.loadSettingsFile(file)
        except Exception as ex:
            pass



    # initialise main interface, triggers and command line parameters
    def __init__(self, module="TracExtract", parent=None, silent=False):
        super(Ui_MainWindow, self).__init__()
        QtGui.QMainWindow.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle("MetExtract II - %s"%module)
        self.checkedLCMSFiles={}

        self.silent=silent


        logging.info("MetExtract II started (module %s)"%module)


        self.preConfigured_adducts = [
                ConfiguredAdduct(name="[M+H]+",       mzoffset=  1.007276, charge=1, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+NH4]+",     mzoffset= 18.033823, charge=1, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+Na]+",      mzoffset= 22.989218, charge=1, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+CH3OH+H]+", mzoffset= 33.033489, charge=1, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+K]+",       mzoffset= 38.963158, charge=1, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+ACN+H]+",   mzoffset= 42.033823, charge=1, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+2Na-H]+",   mzoffset= 44.971160, charge=1, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+2K-H]+",    mzoffset= 76.919040, charge=1, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+CH3FeN]+",  mzoffset= 84.96094 , charge=1, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+2H]++",     mzoffset=  1.007276, charge=2, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+H+NH4]++",  mzoffset=  9.520550, charge=2, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+H+Na]++",   mzoffset= 11.998247, charge=2, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+H+K]++",    mzoffset= 19.985217, charge=2, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+2Na]++",    mzoffset= 22.989218, charge=2, polarity='+', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[2M+H]+",      mzoffset=  1.007276, charge=1, polarity='+', mCount=2, entryType="user"),
                ConfiguredAdduct(name="[2M+NH4]+",    mzoffset= 18.033823, charge=1, polarity='+', mCount=2, entryType="user"),
                ConfiguredAdduct(name="[2M+Na]+",     mzoffset= 22.989218, charge=1, polarity='+', mCount=2, entryType="user"),
                ConfiguredAdduct(name="[2M+K]+",      mzoffset= 38.963158, charge=1, polarity='+', mCount=2, entryType="user"),

                ConfiguredAdduct(name="[M-H2O-H]-",   mzoffset=-19.01839 , charge=1, polarity='-', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M-H]-",       mzoffset=- 1.007276, charge=1, polarity='-', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+Na-2H]-",   mzoffset= 20.974666, charge=1, polarity='-', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+Cl]-",      mzoffset= 34.969402, charge=1, polarity='-', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+K-2H]-",    mzoffset= 36.948606, charge=1, polarity='-', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[M+Br]-",      mzoffset= 78.918885, charge=1, polarity='-', mCount=1, entryType="user"),
                ConfiguredAdduct(name="[2M-H]-",      mzoffset=- 1.007276, charge=1, polarity='-', mCount=2, entryType="user"),
                ConfiguredAdduct(name="[2M+Fa-H]-",   mzoffset= 44.998201, charge=1, polarity='-', mCount=2, entryType="user"),
                ConfiguredAdduct(name="[2M+Hac-H]-",  mzoffset= 59.013851, charge=1, polarity='-', mCount=2, entryType="user")
        ]
        self.adducts=self.preConfigured_adducts

        self.preConfigured_elementsForNL = [
                ConfiguredElement(name="C",  weight=12.,       numberValenzElectrons=4, minCount=0, maxCount=3),
                ConfiguredElement(name="H",  weight=1.00783,   numberValenzElectrons=1, minCount=0, maxCount=30),
                ConfiguredElement(name="O",  weight=15.99491,  numberValenzElectrons=6, minCount=0, maxCount=20),
                ConfiguredElement(name="N",  weight=14.00307,  numberValenzElectrons=5, minCount=0, maxCount=2),
                ConfiguredElement(name="P",  weight=30.97376,  numberValenzElectrons=5, minCount=0, maxCount=2),
                ConfiguredElement(name="S",  weight=31.97207,  numberValenzElectrons=6, minCount=0, maxCount=2),
                ConfiguredElement(name="Cl", weight=34.968852, numberValenzElectrons=7, minCount=0, maxCount=1)
        ]
        self.elementsForNL=self.preConfigured_elementsForNL

        self.preConfigured_heteroElements = [
            ConfiguredHeteroAtom(name='34S',  mzOffset=1.9957960000000021, relativeAbundance=0.044306461797516308,               minCount=1, maxCount=4),
            ConfiguredHeteroAtom(name='37Cl', mzOffset=1.9970499999999944, relativeAbundance=0.31978355549689846,                minCount=1, maxCount=2),
            ConfiguredHeteroAtom(name='41K',  mzOffset= 1.99812,           relativeAbundance=0.07526881720430107526881720430108, minCount=1, maxCount=1)
        ]
        self.heteroElements=self.preConfigured_heteroElements

        self.lastOpenDir = "."
        if os.environ.has_key('USERPROFILE'):
            self.lastOpenDir = os.getenv('USERPROFILE')
        elif os.environ.has_key('HOME'):
            self.lastOpenDir = os.getenv('HOME')

        self.grpFile = None
        self.currentOpenResultsFile=None




        if module=="TracExtract":
            self.labellingExperiment=TRACER
            logging.info("Starting module TracExtract\n")
        elif module=="AllExtract":
            self.labellingExperiment=METABOLOME
            logging.info("Starting module AllExtract\n")
        else:
            logging.error("Error: invalid module '%s' selected.\nPlease specify either 'TracExtract' or 'AllExtract'\n"%module)
            QtGui.QMessageBox.warning(self, "MetExtract", "Error: invalid module '%s' selected.\nPlease specify either 'TracExtract' or 'AllExtract'\n"%module,
                                      QtGui.QMessageBox.Ok)
            raise Exception("Error: invalid module '%s' selected.\nPlease specify either 'TracExtract' or 'AllExtract'")

        self.setModuleUI(self.labellingExperiment)


        # display used R-Version in main interface
        try:
            v = r("R.Version()$version.string")
            module="AllExtract" if self.labellingExperiment==METABOLOME else ("TracExtract" if self.labellingExperiment==TRACER else "")
            self.ui.version.setText("%s II %s [%s]" % (module, MetExtractVersion, str(v)[5:(len(str(v)) - 1)]))
        except:
            self.ui.version.setText("%s [Error: R not available]" % version)

        # limit CPU usage to #cores-1 per default
        cpus = cpu_count()
        if cpus > 1:
            if self.ui.workingCore.isChecked():
                cpus -= 1
            self.ui.cpuCores.setMaximum(cpus)
            self.ui.cpuCores.setValue(cpus)
        else:
            self.ui.workingCore.setVisible(False)
            self.ui.cpuCores.setVisible(False)
            self.ui.label_43.setText("Calculations will be performed on one core only. System might be unresponsive")
            self.ui.cpuCores.setMaximum(1)
            self.ui.cpuCores.setValue(1)

        self.ui.workingCore.toggled.connect(self.updateCores)

        self.ui.tabWidget.setCurrentIndex(0)
        self.ui.tabWidget_2.setCurrentIndex(0)

        self.ui.isotopeAText.textChanged.connect(self.isotopeATextChanged)
        self.ui.isotopeBText.textChanged.connect(self.isotopeBTextChanged)

        self.ui.useRatio.stateChanged.connect(self.useRatioValueStateChanged)
        self.ui.useRatio.setCheckState(QtCore.Qt.Unchecked)
        self.ui.ratioGroupBox.setEnabled(False)

        self.configuredTracers = []
        self.updateTracerInfo()
        self.ui.setupTracers.clicked.connect(self.showTracerEditor)

        self.ui.processIndividualFiles.toggled.connect(self.procIndFilesChanges)
        self.ui.groupResults.toggled.connect(self.groupFilesChanges)
        self.ui.processMultipleFiles.toggled.connect(self.processMultipleFilesChanged)
        self.ui.saveMZXML.toggled.connect(self.saveMZXMLChanged)
        self.updateIndividualFileProcessing = True

        self.ui.alignChromatograms.toggled.connect(self.alignToggled)

        self.ui.groupsSave.setText("./results.tsv")

        self.ui.addGroup.clicked.connect(self.showAddGroupDialogClicked)
        self.ui.saveGroups.clicked.connect(self.saveGroupsClicked)
        self.ui.loadGroups.clicked.connect(self.loadGroupsClicked)
        self.ui.removeGroup.clicked.connect(self.remGrp)

        self.ui.forceFilterLineUpdate.clicked.connect(self.forceUpdateFL)

        self.ui.groupsList.doubleClicked.connect(self.editGroup)

        self.ui.startIdentification.clicked.connect(self.runProcess)

        self.ui.ppmRangeIdentification.valueChanged.connect(self.updateGroupPPM)

        self.ui.isotopeAText.setText("12C")
        self.ui.isotopeBText.setText("13C")

        self.ui.groupsSelectFile.clicked.connect(self.selectGroupsFile)

        self.ui.aboutMenue.triggered.connect(self.aboutMe)
        self.ui.helpMenue.triggered.connect(self.helpMe)
        self.ui.actionLoad_Settings.triggered.connect(self.loadSettings)
        self.ui.actionSave_Settings.triggered.connect(self.saveSettings)
        self.ui.exitMenue.triggered.connect(self.exitMe)

        self.ui.processedFilesComboBox.currentIndexChanged.connect(self.selectedResultChanged)

        self.ui.res_ExtractedData.itemSelectionChanged.connect(self.selectedResChanged)

        self.ui.resultsExperiment_TreeWidget.itemSelectionChanged.connect(self.resultsExperimentChanged)


        self.ui.eicSmoothingWindow.currentIndexChanged.connect(self.smoothingWindowChanged)
        self.smoothingWindowChanged()

        self.ui.isoAbundance.stateChanged.connect(self.isoSearchChanged)
        self.isoSearchChanged()

        # setup result plots
        #Setup first plot
        #http://eli.thegreenplace.net/2009/01/20/matplotlib-with-pyqt-guis/
        self.ui.pl_tic = QtCore.QObject()
        self.ui.pl_tic.dpi = 50
        self.ui.pl_tic.fig = Figure((5.0, 4.0), dpi=self.ui.pl_tic.dpi, facecolor='white')
        self.ui.pl_tic.fig.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.95)
        self.ui.pl_tic.canvas = FigureCanvas(self.ui.pl_tic.fig)
        self.ui.pl_tic.canvas.setParent(self.ui.visualizationWidget)
        self.ui.pl_tic.axes = self.ui.pl_tic.fig.add_subplot(111)
        simpleaxis(self.ui.pl_tic.axes)
        self.ui.pl_tic.twinxs = [self.ui.pl_tic.axes]
        self.ui.pl_tic.mpl_toolbar = NavigationToolbar(self.ui.pl_tic.canvas, self.ui.ticVisualisationWidget)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.ui.pl_tic.canvas)
        vbox.addWidget(self.ui.pl_tic.mpl_toolbar)
        self.ui.ticVisualisationWidget.setLayout(vbox)

        #Setup first plot
        #http://eli.thegreenplace.net/2009/01/20/matplotlib-with-pyqt-guis/
        self.ui.pl1 = QtCore.QObject()
        self.ui.pl1.dpi = 50
        self.ui.pl1.fig = Figure((5.0, 4.0), dpi=self.ui.pl1.dpi, facecolor='white')
        self.ui.pl1.fig.subplots_adjust(left=0.05, bottom=0.1, right=0.99, top=0.95)
        self.ui.pl1.canvas = FigureCanvas(self.ui.pl1.fig)
        self.ui.pl1.canvas.setParent(self.ui.visualizationWidget)
        self.ui.pl1.axes = self.ui.pl1.fig.add_subplot(111)
        simpleaxis(self.ui.pl1.axes)
        self.ui.pl1.twinxs = [self.ui.pl1.axes]
        self.ui.pl1.mpl_toolbar = NavigationToolbar(self.ui.pl1.canvas, self.ui.visualizationWidget)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.ui.pl1.canvas)
        vbox.addWidget(self.ui.pl1.mpl_toolbar)
        self.ui.visualizationWidget.setLayout(vbox)

        #setup second plot
        self.ui.pl2 = QtCore.QObject()
        self.ui.pl2.type = None
        self.ui.pl2.dpi = 50
        self.ui.pl2.fig = Figure((5.0, 4.0), dpi=self.ui.pl2.dpi, facecolor='white')
        self.ui.pl2.fig.subplots_adjust(left=0.15, bottom=0.05, right=0.99, top=0.85)
        self.ui.pl2.canvas = FigureCanvas(self.ui.pl2.fig)
        self.ui.pl2.canvas.setParent(self.ui.visualizationWidget2)
        self.ui.pl2.axes = self.ui.pl2.fig.add_subplot(111)
        #noaxis(self.ui.pl2.axes)
        self.ui.pl2.twinxs = [self.ui.pl2.axes]
        self.ui.pl2.mpl_toolbar = NavigationToolbar(self.ui.pl2.canvas, self.ui.visualizationWidget2)
        self.ui.pl2.pictureShown = False

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.ui.pl2.canvas)
        vbox.addWidget(self.ui.pl2.mpl_toolbar)
        self.ui.visualizationWidget2.setLayout(vbox)

        #setup third plot
        self.ui.pl3 = QtCore.QObject()
        self.ui.pl3.type = None
        self.ui.pl3.dpi = 50
        self.ui.pl3.fig = Figure((5.0, 4.0), dpi=self.ui.pl3.dpi, facecolor='white')
        self.ui.pl3.fig.subplots_adjust(left=0.05, bottom=0.1, right=0.99, top=0.95)
        self.ui.pl3.canvas = FigureCanvas(self.ui.pl3.fig)
        self.ui.pl3.canvas.setParent(self.ui.visualizationWidget3)
        self.ui.pl3.axes = self.ui.pl3.fig.add_subplot(111)
        simpleaxis(self.ui.pl3.axes)
        self.ui.pl3.twinxs = [self.ui.pl3.axes]
        self.ui.pl3.mpl_toolbar = NavigationToolbar(self.ui.pl3.canvas, self.ui.visualizationWidget3)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.ui.pl3.canvas)
        vbox.addWidget(self.ui.pl3.mpl_toolbar)
        self.ui.visualizationWidget3.setLayout(vbox)

        #Setup experiment plot - overlaid EICs plot
        #http://eli.thegreenplace.net/2009/01/20/matplotlib-with-pyqt-guis/
        self.ui.resultsExperiment_plot = QtCore.QObject()
        self.ui.resultsExperiment_plot.dpi = 50
        self.ui.resultsExperiment_plot.fig = Figure((5.0, 4.0), dpi=self.ui.resultsExperiment_plot.dpi, facecolor='white')
        self.ui.resultsExperiment_plot.fig.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.95)
        self.ui.resultsExperiment_plot.canvas = FigureCanvas(self.ui.resultsExperiment_plot.fig)
        self.ui.resultsExperiment_plot.canvas.setParent(self.ui.visualizationWidget)
        self.ui.resultsExperiment_plot.axes = self.ui.resultsExperiment_plot.fig.add_subplot(111)
        simpleaxis(self.ui.resultsExperiment_plot.axes)
        self.ui.resultsExperiment_plot.twinxs = [self.ui.resultsExperiment_plot.axes]
        self.ui.resultsExperiment_plot.mpl_toolbar = NavigationToolbar(self.ui.resultsExperiment_plot.canvas, self.ui.visualizationWidget)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.ui.resultsExperiment_plot.canvas)
        vbox.addWidget(self.ui.resultsExperiment_plot.mpl_toolbar)
        self.ui.resultsExperiment_widget.setLayout(vbox)

        #Setup experiment plot - separate chrom. peaks plot
        #http://eli.thegreenplace.net/2009/01/20/matplotlib-with-pyqt-guis/
        self.ui.resultsExperimentSeparatedPeaks_plot = QtCore.QObject()
        self.ui.resultsExperimentSeparatedPeaks_plot.dpi = 50
        self.ui.resultsExperimentSeparatedPeaks_plot.fig = Figure((5.0, 4.0), dpi=self.ui.resultsExperimentSeparatedPeaks_plot.dpi, facecolor='white')
        self.ui.resultsExperimentSeparatedPeaks_plot.fig.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.95)
        self.ui.resultsExperimentSeparatedPeaks_plot.canvas = FigureCanvas(self.ui.resultsExperimentSeparatedPeaks_plot.fig)
        self.ui.resultsExperimentSeparatedPeaks_plot.canvas.setParent(self.ui.visualizationWidget)
        self.ui.resultsExperimentSeparatedPeaks_plot.axes = self.ui.resultsExperimentSeparatedPeaks_plot.fig.add_subplot(111)
        simpleaxis(self.ui.resultsExperimentSeparatedPeaks_plot.axes)
        self.ui.resultsExperimentSeparatedPeaks_plot.twinxs = [self.ui.resultsExperimentSeparatedPeaks_plot.axes]
        self.ui.resultsExperimentSeparatedPeaks_plot.mpl_toolbar = NavigationToolbar(self.ui.resultsExperimentSeparatedPeaks_plot.canvas, self.ui.visualizationWidget)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.ui.resultsExperimentSeparatedPeaks_plot.canvas)
        vbox.addWidget(self.ui.resultsExperimentSeparatedPeaks_plot.mpl_toolbar)
        self.ui.resultsExperimentSeparatedPeaks_widget.setLayout(vbox)

        #Setup experiment plot - separate chrom. peaks plot
        #http://eli.thegreenplace.net/2009/01/20/matplotlib-with-pyqt-guis/
        self.ui.resultsExperimentMSScanPeaks_plot = QtCore.QObject()
        self.ui.resultsExperimentMSScanPeaks_plot.dpi = 50
        self.ui.resultsExperimentMSScanPeaks_plot.fig = Figure((5.0, 4.0), dpi=self.ui.resultsExperimentMSScanPeaks_plot.dpi, facecolor='white')
        self.ui.resultsExperimentMSScanPeaks_plot.fig.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.95)
        self.ui.resultsExperimentMSScanPeaks_plot.canvas = FigureCanvas(self.ui.resultsExperimentMSScanPeaks_plot.fig)
        self.ui.resultsExperimentMSScanPeaks_plot.canvas.setParent(self.ui.visualizationWidget)
        self.ui.resultsExperimentMSScanPeaks_plot.axes = self.ui.resultsExperimentMSScanPeaks_plot.fig.add_subplot(111)
        simpleaxis(self.ui.resultsExperimentMSScanPeaks_plot.axes)
        self.ui.resultsExperimentMSScanPeaks_plot.twinxs = [self.ui.resultsExperimentMSScanPeaks_plot.axes]
        self.ui.resultsExperimentMSScanPeaks_plot.mpl_toolbar = NavigationToolbar(self.ui.resultsExperimentMSScanPeaks_plot.canvas, self.ui.visualizationWidget)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.ui.resultsExperimentMSScanPeaks_plot.canvas)
        vbox.addWidget(self.ui.resultsExperimentMSScanPeaks_plot.mpl_toolbar)
        self.ui.resultsExperimentMSScan_widget.setLayout(vbox)

        self.ui.pl2.xics = []
        self.ui.pl2.times = []
        self.ui.pl2.peaks = []
        self.ui.pl2.x_vals = []
        self.ui.pl2.y_vals = []

        font = {'size': 18}
        matplotlib.rc('font', **font)


        self.ui.dataFilter.textChanged.connect(self.filterEdited)
        self.ui.res_ExtractedData.itemDoubleClicked.connect(self.res_doubleClick)


        self.ui.openRAW.clicked.connect(self.openRawFile)
        self.ui.setChromPeakName.clicked.connect(self.setChromPeakName)

        self.ui.res_ExtractedData.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.ui.res_ExtractedData.customContextMenuRequested.connect(self.showPopup)

        self.ui.relationshipConfig.clicked.connect(self.relationShipConfiguration)

        self.ui.defineHeteroAtoms.clicked.connect(self.heteroAtomsConfiguration)

        self.ui.visualConfig.setVisible(False)
        self.ui.scaleFeatures.setVisible(True)
        self.ui.correctcCount.setVisible(False)
        self.ui.snrThLabel.setVisible(False)
        self.ui.wavelet_SNRThreshold.setVisible(False)
        self.ui.scaleError.setVisible(False)
        self.ui.peak_scaleError.setVisible(False)

        #todo: implement results view with all brackated feature pairs
        self.ui.tabWidget.setTabText(3, "EXPERIMENTAL: "+self.ui.tabWidget.tabText(3))
        self.ui.tabWidget.setTabEnabled(3, True)

        # fetch provided settings and show them as a menu
        try:
            try:
                self.loadSettingsFile(get_main_dir()+"/Settings/defaultSettings.ini", checkExperimentType=False)
            except:
                logging.warning("Warning: Default settings are invalid. Skipping..")
            if os.path.exists(get_main_dir()+"/Settings"):
                menus = {}
                menus[get_main_dir()]="Predefined"
                for dirname, dirnames, filenames in os.walk(get_main_dir()+'/Settings/'):
                    for filename in filenames:

                        toMenu = self.ui.menuAvailable_Settings
                        cDir = get_main_dir()

                        for dir in filter(None, dirname[len(get_main_dir()):].replace("\\", "/").split("/")):

                            cDir = cDir + "/" + dir
                            if cDir not in menus:
                                menus[cDir] = QtGui.QMenu(toMenu)
                                menus[cDir].setTitle(dir)
                                toMenu.addMenu(menus[cDir])

                            toMenu = menus[cDir]

                        ac = QtGui.QAction(self.ui.menuAvailable_Settings)
                        ac.setText(filename.replace(".ini",""))

                        ac.triggered.connect(functools.partial(self.loadAvailableSettingsFile,os.path.join(dirname, filename).replace("\\", "/")))

                        toMenu.addAction(ac)
        except Exception as ex:
            import traceback
            traceback.print_exc()
            logging.error(str(traceback))

            logging.error("Error in %s: %s" % (self.file, str(ex)))
            QtGui.QMessageBox.information(self, "MetExtract", "Cannot load settings", QtGui.QMessageBox.Ok)





if __name__ == '__main__':
    # add freeze support (required for multiprocessing)
    freeze_support()

    # parse supplied options
    parser = OptionParser()
    parser.add_option("-m", "--module", dest="module", default="TracExtract", metavar="MODULE",
                      help="Select 'AllExtract' or 'TracExtract' module (default TracExtract)")
    parser.add_option("-g", "--groupFile", dest="groupdest", default=None, metavar="GROUP-FILE",
                      help="Load group compilation (use -l group to load associated settings automatically)")
    parser.add_option("-l", "--settingsFile", dest="settings", default=None, metavar="SETTINGS-FILE",
                      help="Load settings file. Specify 'group' to load settings from group file or provide path to a settings file ")

    parser.add_option("-u", "--processIndFiles", dest="processIndFiles", default=None,
                      help="Set if individual files should be processed")
    parser.add_option("-i", "--processMultipleFiles", dest="processMultFiles", default=None,
                      help="Set if multiple files should be processed")
    parser.add_option("-o", "--groupResults", dest="groupResults", default=None,
                      help="Set if individual results shall be grouped")
    parser.add_option("-p", "--reintegrateResults", dest="reintegrateResults", default=None,
                      help="Set if individual results shall be reintegrated")

    parser.add_option("-c", "--cores", dest="cores", default=1000, metavar="CORES", type="int",
                      help="Set maximum number of cores to use")
    parser.add_option("-s", "--startAutomatically", dest="start", default=False, action="store_true",
                      help="Start calculations for loaded LC-HRMS compilation automatically")
    parser.add_option("-e", "--exit", dest="exit", default=False, action="store_true",
                      help="Exit after calculations have been performed (only allowed in combination with -s)")
    parser.add_option("-x", "--silentStart", dest="silentStart", default=False, action="store_true",
                      help="Start MetExtract without any meta-output")

    (opts, args) = parser.parse_args()

    # check platform
    import platform
    import ctypes

    # inform user, what may not work on the current platform
    if not opts.silentStart:
        logging.info("Platform:")
        if platform.system() == "Darwin":  #MAC
            logging.info("  MAC OS detected")
        if platform.system() == "Windows":  #Windows
            logging.info("  Windows OS detected: Large LC-HRMS files are currently not supported (eg. Q-ToF in 4Ghz mode)")
        if platform.system() == "Linux":  #Linux
            logging.info("  Linux OS detected")

        if ctypes.sizeof(ctypes.c_voidp) == 4 and platform.system() != "Windows":  #32 bit system
            logging.info("  32-bit environment detected: Large files may not work correctly")

        logging.info("")

    #start PyQT GUI application
    app = QtGui.QApplication(sys.argv)

    #search for R-packages and install them if necessary
    if not opts.silentStart:
        logging.info("R-Configuration: ")
    rAvailable = False
    rPackagesAvailable = False
    try:
        try:
            import os
            logging.info("  R_HOME: %s"%(os.environ["R_HOME"].strip()))
            logging.info("  R_HOME_FROM: %s"%(os.environ["R_HOME_FROM"].strip()))
        except:
            pass
        import rpy2.robjects as ro

        r = ro.r

        v = r("R.Version()$version.string")
        if not opts.silentStart:
            logging.info("  %s" % (str(v)[5:(len(str(v)) - 1)]))

        rAvailable = True

        if checkRDependencies(r):  # missing dependencies
            QtGui.QMessageBox.warning(None, "MetExtract", "Error: Not all R-dependencies are available or can be installed", QtGui.QMessageBox.Ok)
            logging.warning("  Required R-packages are not available or installation failed")
        else:
            if not opts.silentStart:
                logging.info("  All necessary R-dependencies are installed/available")
            rPackagesAvailable = True

    except:
        logging.info("  Error: R could not be loaded..")
    logging.info("")

    if rAvailable and rPackagesAvailable:

        # show main window
        try:
            mainWin = mainWindow(module=opts.module, silent=opts.silentStart)

            if USEGRADIENTDESCENDPEAKPICKING:
                p = mainWin.palette()
                p.setColor(mainWin.backgroundRole(), QtCore.Qt.red)
                mainWin.setPalette(p)
                if not opts.start:
                    QtGui.QMessageBox.warning(None, "MetExtract",
                                  "WARNING: Gradient descend algorithm for chromatographic peak picking is selected",
                                  QtGui.QMessageBox.Ok)
        except:
            mainWin=None

        if mainWin is not None:
            mainWin.show()


            if opts.groupdest is not None:
                groupFile = opts.groupdest.replace("\\", "/")
                logging.info("Loading LC-HRMS data compilation from '%s'" % groupFile)
                mainWin.loadGroups(groupFile, forceLoadSettings=False, askLoadSettings=False)
                mainWin.lastOpenDir = groupFile[0:(groupFile.rfind("/") + 1)]

            if opts.settings is not None:
                settFile = opts.settings
                if opts.settings == 'group':
                    settFile = opts.groupdest
                settFile = settFile.replace("\\", "/")
                logging.info("Loading settings from '%s'" % settFile)
                mainWin.loadSettingsFile(settFile, silent=True)

            mainWin.ui.cpuCores.setValue(opts.cores)


            if opts.processIndFiles is not None:
                mainWin.ui.processIndividualFiles.setChecked(opts.processIndFiles.lower() in ["true", "t", "yes", "y", "1"])

            if opts.processMultFiles is not None:
                mainWin.ui.processMultipleFiles.setChecked(opts.processMultFiles.lower() in ["true", "t", "yes", "y", "1"])

            if opts.groupResults is not None:
                mainWin.ui.groupResults.setChecked(opts.groupResults.lower() in ["true", "t", "yes", "y", "1"])

            if opts.reintegrateResults is not None:
                mainWin.ui.integratedMissedPeaks.setChecked(opts.reintegrateResults.lower() in ["true", "t", "yes", "y", "1"])


            # print development messages
            if False:
                import textwrap

                messages={"Severe":[], "Error":[], "Warning":[], "Info":[]}
                messages["Severe"].append("Metabolite cluster calculation is in beta testing. Please review the process and report issues!")
                messages["Info"].append("Beta version")
                messages["Info"].append("Remove those messages (and the GUI-component)")
                messages["Info"].append("Add gradient descend peak picking")

                allMessages=[]

                logging.info("   -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")
                logging.info("   -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")

                for k in sorted(messages.keys()):
                    vs=messages[k]
                    if len(vs)>0:
                        allMessages.append("<span style=\"color: yellow\">%s</span>:"%k)
                    for v in sorted(vs):
                        logging.info("   -=-=-  "+"\033[91m"+"%-8s"%k+"\033[0m",)
                        k=""
                        logging.info(textwrap.fill(v).replace("\r", ""))
                        allMessages.append(v)

                logging.info("   -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")
                logging.info("   -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=")

                mainWin.ui.INFOLabel.setText("<div style=\"background-color: red\">%s</div>"%(", ".join(allMessages)))
            else:
                mainWin.ui.INFOLabel.setVisible(False)

            if opts.start:
                mainWin.runProcess(dontSave=True, askStarting=False)





            if opts.exit:
                QtGui.QApplication.exit()
            else:
                sys.exit(app.exec_())










