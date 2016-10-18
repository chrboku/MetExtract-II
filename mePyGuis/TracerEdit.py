import sys
import os
import pickle
import base64
from copy import copy, deepcopy

from PyQt4 import QtGui, QtCore

from mePyGuis.TracerEditor import Ui_Dialog
from utils import getRatio, getXCombinations
from formulaTools import formulaTools, getIsotopeMass



class ConfiguredTracer():
    def __init__(self, name="", elementCount=15, isotopeA="12C", isotopeB="13C", enrichmentA=.9893, enrichmentB=.995,
                 amountA=.9, amountB=1., monoisotopicRatio=1, maxRelNegBias=70, maxRelPosBias=130, tracerType="user",
                 id=-1, mzDelta=None):
        self.id = id

        self.name = name  # 0
        self.elementCount = elementCount  # 1
        self.isotopeA = isotopeA  # 2
        self.isotopeB = isotopeB  # 3
        self.enrichmentA = enrichmentA  # 4
        self.enrichmentB = enrichmentB  # 5
        self.amountA = amountA  # 6
        self.amountB = amountB  # 7
        self.monoisotopicRatio = monoisotopicRatio  # 8
        self.maxRelNegBias = maxRelNegBias  # 9
        self.maxRelPosBias = maxRelPosBias  # 10
        self.tracerType = tracerType  # 11

        if mzDelta is None:
            self.mzDelta = getIsotopeMass(self.isotopeB)[0] - getIsotopeMass(self.isotopeA)[0]
        else:
            self.mzDelta = mzDelta

    def __str__(self):
        return "ConfiguredTracer: %s %s %s %d  enrichment: %.3f %.3f amount %.1f %.1f monoisotopicRatio %.3f bias: %.1f %.1f tracerType %s" % (
            self.name, self.isotopeA, self.isotopeB, self.elementCount, self.enrichmentA, self.enrichmentB,
            self.amountA, self.amountB, self.monoisotopicRatio, self.maxRelNegBias, self.maxRelPosBias, self.tracerType)



def getShortName(element):
    fT = formulaTools()
    for x in fT.elemDetails:
        d = fT.elemDetails[x]
        if d[0] == element:
            return x
    return ""




class TracerTableModel(QtCore.QAbstractTableModel):
    def __init__(self, datain, headerdata, parent=None, formTool=None, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)
        self.arraydata = datain
        self.headerdata = headerdata

    def rowCount(self, parent):
        return len(self.arraydata)

    def columnCount(self, parent):
        return 12

    def getRawData(self, row):
        return self.arraydata[row]

    def data(self, index, role):
        if not index.isValid():
            return QtCore.QVariant()
        elif role != QtCore.Qt.DisplayRole and role != QtCore.Qt.EditRole:
            return QtCore.QVariant()

        if self.arraydata[index.row()].tracerType == "user":

            #print self.arraydata[index.row()]

            if index.row() == (len(self.arraydata) - 1):
                return QtCore.QString("")
            if index.column() == 0:
                return QtCore.QString(self.arraydata[index.row()].name)
            if index.column() == 2:
                return QtCore.QString(self.arraydata[index.row()].isotopeA)
            if index.column() == 3:
                return QtCore.QString(self.arraydata[index.row()].isotopeB)
            if index.column() == 11:
                return QtCore.QString(self.arraydata[index.row()].tracerType)

            if index.column() == 1:
                if self.arraydata[index.row()].elementCount == 0:
                    return QtCore.QString("")
                return QtCore.QVariant(self.arraydata[index.row()].elementCount)
            if index.column() == 6:
                if self.arraydata[index.row()].amountA == 0:
                    return QtCore.QString("")
                return QtCore.QVariant(self.arraydata[index.row()].amountA)
            if index.column() == 7:
                if self.arraydata[index.row()].amountB == 0:
                    return QtCore.QString("")
                return QtCore.QVariant(self.arraydata[index.row()].amountB)

            if index.column() in [4, 5, 8]:
                if index.column() == 8:
                    if len(self.arraydata[index.row()].isotopeA) > 0 and len(
                            self.arraydata[index.row()].isotopeB) > 0 and self.arraydata[
                        index.row()].elementCount > 0 and self.arraydata[index.row()].enrichmentA > 0 and \
                                    self.arraydata[index.row()].enrichmentB > 0 and self.arraydata[
                        index.row()].amountA > 0 and self.arraydata[index.row()].amountB > 0:
                        a = getRatio(self.arraydata[index.row()].enrichmentA, self.arraydata[index.row()].elementCount,
                                     0) * self.arraydata[index.row()].amountA
                        b = getRatio(self.arraydata[index.row()].enrichmentB, self.arraydata[index.row()].elementCount,
                                     0) * self.arraydata[index.row()].amountB
                        ratio = a / b
                        self.arraydata[index.row()].monoisotopicRatio = ratio

                        return QtCore.QString("%.3f" % self.arraydata[index.row()].monoisotopicRatio)

                if self.arraydata[index.row()].name == 0:
                    return QtCore.QString("")

                if role == QtCore.Qt.DisplayRole and index.column() == 4:
                    return QtCore.QString("%.2f %%" % (self.arraydata[index.row()].enrichmentA * 100.))
                if role == QtCore.Qt.DisplayRole and index.column() == 5:
                    return QtCore.QString("%.2f %%" % (self.arraydata[index.row()].enrichmentB * 100.))


                elif role == QtCore.Qt.DisplayRole and index.column() == 8:
                    return QtCore.QString("%.2f %%" % self.arraydata[index.row()].monoisotopicRatio)
                elif role == QtCore.Qt.EditRole and index.column() == 4:
                    return QtCore.QVariant(self.arraydata[index.row()].enrichmentA * 100.)
                elif role == QtCore.Qt.EditRole and index.column() == 5:
                    return QtCore.QVariant(self.arraydata[index.row()].enrichmentB * 100.)
                elif role == QtCore.Qt.EditRole and index.column() == 8:
                    return QtCore.QVariant()

            if index.column() == 4:
                if self.arraydata[index.row()].enrichmentA == 0:
                    return QtCore.QString("")
                if role == QtCore.Qt.DisplayRole:
                    return QtCore.QString("%d" % (self.arraydata[index.row()].enrichmentA * 100.))
                elif role == QtCore.Qt.EditRole:
                    return QtCore.QVariant(self.arraydata[index.row()].enrichmentA * 100.)

            if index.column() == 5:
                if self.arraydata[index.row()].enrichmentB == 0:
                    return QtCore.QString("")
                if role == QtCore.Qt.DisplayRole:
                    return QtCore.QString("%d" % (self.arraydata[index.row()].enrichmentB * 100.))
                elif role == QtCore.Qt.EditRole:
                    return QtCore.QVariant(self.arraydata[index.row()].enrichmentB * 100.)

            if index.column() == 9:
                if self.arraydata[index.row()].maxRelNegBias == 0:
                    return QtCore.QString("0 %")
                if role == QtCore.Qt.DisplayRole:
                    return QtCore.QString("%.3f %%" % (self.arraydata[index.row()].maxRelNegBias * 100.))
                elif role == QtCore.Qt.EditRole:
                    return QtCore.QVariant(self.arraydata[index.row()].maxRelNegBias * 100., )
            if index.column() == 10:
                if self.arraydata[index.row()].maxRelPosBias == 0:
                    return QtCore.QString("")
                if role == QtCore.Qt.DisplayRole:
                    return QtCore.QString("%.3f %%" % (self.arraydata[index.row()].maxRelPosBias * 100.))
                elif role == QtCore.Qt.EditRole:
                    return QtCore.QVariant(self.arraydata[index.row()].maxRelPosBias * 100.)

            elif self.arraydata[index.row()].tracerType == "calc":
                if index.column() == 0:
                    return QtCore.QString(self.arraydata[index.row()].name)
                if index.column() == 11:
                    return QtCore.QString(self.arraydata[index.row()].tracerType)
                if index.column() == 8:
                    return QtCore.QString("%.3f" % self.arraydata[index.row()].monoisotopicRatio)
                if index.column() == 9:
                    return QtCore.QString("%.1f %%" % (self.arraydata[index.row()].maxRelNegBias * 100.))
                if index.column() == 10:
                    return QtCore.QString("%.1f %%" % (self.arraydata[index.row()].maxRelPosBias * 100.))

        return QtCore.QString("")

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return QtCore.QVariant(self.headerdata[col])
        return QtCore.QVariant()

    def insertRows(self, row, count, parent=QtCore.QModelIndex(), type="empty", tracerToInsert=None):
        self.beginInsertRows(parent, row, row + count - 1)
        if tracerToInsert is None:
            tracerToInsert = ConfiguredTracer(name="", tracerType=type)
        self.arraydata.append(tracerToInsert)
        self.endInsertRows()
        return True

    def removeRows(self, row, count, parent=QtCore.QModelIndex()):
        if len(self.arraydata) > row and len(self.arraydata) > (row + count):
            self.beginRemoveRows(parent, row, row + count - 1)
            for i in range(count):
                self.arraydata.pop(row)
            self.endRemoveRows()
            return True
        else:
            return False

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        f = 0
        if index.row() == (len(self.arraydata) - 1) and not (value == ""):
            self.insertRows(len(self.arraydata) + 1, 1, type="empty")
            self.arraydata[len(self.arraydata) - 2].tracerType = "user"

        if index.column() == 0:
            self.arraydata[index.row()].name = str(value.toString())
            return True
        if index.column() == 2:
            self.arraydata[index.row()].isotopeA = str(value.toString())
            return True
        if index.column() == 3:
            self.arraydata[index.row()].isotopeB = str(value.toString())
            return True
        elif index.column() in [4, 5, 6, 7, 9, 10]:
            try:
                f, ok = value.toDouble()
                if not ok:
                    return False
                if index.column() == 4:
                    self.arraydata[index.row()].enrichmentA = f / 100.
                if index.column() == 5:
                    self.arraydata[index.row()].enrichmentB = f / 100.
                if index.column() == 6:
                    self.arraydata[index.row()].amountA = f
                if index.column() == 7:
                    self.arraydata[index.row()].amountB = f
                if index.column() == 9:
                    self.arraydata[index.row()].maxRelNegBias = f / 100.
                if index.column() == 10:
                    self.arraydata[index.row()].maxRelPosBias = f / 100.
                return True
            except:
                return False
        elif index.column() == 1:
            try:
                f, ok = value.toInt()
                if not ok:
                    return False
                self.arraydata[index.row()].elementCount = f
                return True
            except:
                return False

        return False

    def flags(self, index):

        if index.column() in [8, 11]:
            return QtCore.Qt.ItemIsSelectable

        if self.arraydata[index.row()].tracerType == "calc":
            if index.column() in [0, 9, 10]:
                return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable
            else:
                return QtCore.Qt.ItemIsSelectable
        if index.column() not in [8, 11]:
            return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable
        return QtCore.Qt.ItemIsEnabled


def calculateMultiTracerConjugateRatio(trcs):
    l = 1.
    r = 1.
    for i in range(len(trcs)):
        tr = trcs[i]

        l = l * getRatio(tr.enrichmentA / 100, tr.elementCount, 0) * tr.amountA / (
            getRatio(tr.enrichmentA / 100, tr.elementCount, 0) * tr.amountA + getRatio(tr.enrichmentB / 100,
                                                                                       tr.elementCount, 0) * tr.amountB)
        r = r * getRatio(tr.enrichmentB / 100, tr.elementCount, 0) * tr.amountB / (
            getRatio(tr.enrichmentA / 100, tr.elementCount, 0) * tr.amountA + getRatio(tr.enrichmentB / 100,
                                                                                       tr.elementCount, 0) * tr.amountB)

    return l / r


class tracerEdit(QtGui.QDialog, Ui_Dialog):
    def __init__(self, parent=None, initDir=None):
        QtGui.QDialog.__init__(self, parent)
        self.setupUi(self)
        self.setWindowTitle("Tracer editor")

        #self.data = [ConfiguredTracer(name="DON", elementCount=15, isotopeA="12C", isotopeB="13C", enrichmentA=.9893,
        #                              enrichmentB=.995, amountA=50, amountB=50, maxRelNegBias=.4, maxRelPosBias=.4)]
        #self.data.append(
        #    ConfiguredTracer(name="T2/HT2", elementCount=24, isotopeA="12C", isotopeB="13C", enrichmentA=.9893,
        #                     enrichmentB=.995, amountA=1.3, amountB=1., monoisotopicRatio=0, maxRelNegBias=.3,
        #                     maxRelPosBias=.3, tracerType="user"))
        self.data=[ConfiguredTracer(name="", tracerType="empty")]
        self.headers = ["Tracer", "Element count", "Isotope N", "Isotope L", "Natural abundance isotope N (%)",
                        "Enriched abundance isotope L (%)", "Amount N", "Amount L",
                        "Ratio (N:L)", "Min. ratio (%)",
                        "Max. ratio (%)", "Type"]

        self.tracerModel = TracerTableModel(self.data, self.headers)
        self.tracers.setModel(self.tracerModel)

        self.tracers.resizeColumnsToContents()

        self.tracers.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tracers.customContextMenuRequested.connect(self.showPopupTracer)

        self.discardButton.clicked.connect(self.dialogCan)
        self.acceptButton.clicked.connect(self.dialogFin)

        self.saveConfiguration.clicked.connect(self.save)
        self.loadConfiguration.clicked.connect(self.load)

        self.calculateMultiTracerConjugates.clicked.connect(self.calcMultiTracerConjugates)
        self.calculateMultiTracerConjugates.setVisible(False)

        self.checkRatios.setVisible(False)

        self.showMSScan.clicked.connect(self.plotMSScan)

        showSimulator = False
        if not showSimulator:
            self.line_5.setVisible(False)
            self.showMSScan.setVisible(False)
            self.normaliseRatio.setVisible(False)

        if initDir is None:
            self.initDir = "."
            if os.environ.has_key('USERPROFILE'):
                self.initDir = os.getenv('USERPROFILE')
            elif os.environ.has_key('HOME'):
                self.initDir = os.getenv('HOME')

        else:
            self.initDir = initDir

        self.acceptButton.setFocus(True)

    def calcMultiTracerConjugates(self):
        trs = []

        rm = set()
        for i in range(len(self.data)):
            d = self.data[i]
            if d.tracerType == "calc":
                rm.add(i)
        for i in sorted([f for f in rm], reverse=True):
            self.tracerModel.removeRows(i, 1)

        for d in self.data:
            if d.tracerType == "user":
                trs.append(d)

        if len(trs) == 1:
            QtGui.QMessageBox.information(self, "MetExtract",
                                          "Only one tracer has been specified. No multi-tracer conjugates are expected",
                                          QtGui.QMessageBox.Ok)
        else:
            x = getXCombinations(trs, 2)

            for i in x:
                print " & ".join([j.name for j in i])
                calcTracer = ConfiguredTracer(name=" & ".join([j.name for j in i]), tracerType="user")
                self.tracerModel.insertRows(row=len(self.data) + 1, count=1, tracerToInsert=calcTracer)

                calcTracer.monoisotopicRatio = calculateMultiTracerConjugateRatio(i)
                calcTracer.maxRelNegBias = 30
                calcTracer.maxRelPosBias = 30
                calcTracer.tracerType = "calc"


    def plotMSScan(self):
        import matplotlib
        import matplotlib.pyplot as plt

        def simpleaxis(ax):
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()

        fig = plt.figure(facecolor="white")
        ax = fig.add_subplot(111)
        simpleaxis(ax)

        m = 1
        for x in self.tracers.selectedIndexes():
            dat = self.tracerModel.getRawData(x.row())
            from utils import getNormRatio

            for i in range(dat.elementCount):
                x = i
                if self.normaliseRatio.isChecked():
                    y = getNormRatio(dat.enrichmentA / 100., dat.elementCount, i) * dat.amountA
                else:
                    y = getRatio(dat.enrichmentA / 100., dat.elementCount, i) * dat.amountA

                ax.plot([x, x], [0, y], color="red")
                x = i
                if self.normaliseRatio.isChecked():
                    y = getNormRatio(dat.enrichmentB / 100., dat.elementCount, i) * dat.amountB
                else:
                    y = getRatio(dat.enrichmentB / 100., dat.elementCount, i) * dat.amountB
                ax.plot([dat.elementCount - i, dat.elementCount - i], [0, y], color="blue")

            m = max(m, dat.elementCount)

        ax.set_xlim([-1, m + 1])

        fig.show()

    def showPopupTracer(self, position):
        menu = QtGui.QMenu()
        quitAction = menu.addAction("Delete")
        action = menu.exec_(self.tracers.mapToGlobal(position))
        if action == quitAction:
            for x in self.tracers.selectedIndexes():
                self.tracerModel.removeRows(x.row(), 1)

    def load(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, "Select tracer file", self.initDir,
                                                     "tracer file (*.trc);;All files (*.*)")
        if len(filename) > 0:
            self.initDir = str(filename).replace("\\", "/")
            self.initDir = self.initDir[:self.initDir.rfind("/")]
        if len(filename) > 0:
            self.loadTracerFile(filename)

    def save(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, "Select tracer file", self.initDir,
                                                     filter="tracer file (*.trc);;All files(*.*)")
        if len(filename) > 0:
            self.initDir = str(filename).replace("\\", "/")
            self.initDir = self.initDir[:self.initDir.rfind("/")]
        if len(filename) > 0:
            self.saveTracerFile(filename)

    def saveTracerFile(self, filename):
        conf = QtCore.QSettings(filename, QtCore.QSettings.IniFormat)
        conf.setValue("Tracers", base64.b64encode(pickle.dumps(self.getTracers())))

    def loadTracerFile(self, filename):
        conf = QtCore.QSettings(filename, QtCore.QSettings.IniFormat)
        tracers = pickle.loads(base64.b64decode(str(conf.value("Tracers").toString())))

        self.data = tracers
        self.data.append(ConfiguredTracer(name="", tracerType="empty"))
        self.tracerModel = TracerTableModel(self.data, self.headers)
        self.tracers.setModel(self.tracerModel)
        return tracers

    def setTracers(self, data):
        self.data = []
        for d in data:
            self.data.append(d)

        self.data = data
        self.data.append(ConfiguredTracer(name="", tracerType="empty"))
        self.tracerModel = TracerTableModel(self.data, self.headers)
        self.tracers.setModel(self.tracerModel)

    def getTracers(self):
        x = deepcopy([d for d in self.data if d.tracerType in ["user", "calc"]])

        return x

    def checkConfiguration(self):
        ok = True
        severe = False
        tracers = self.getTracers()
        last = len(tracers) - 1
        if self.checkRatios.checkState() == QtCore.Qt.Checked:
            for i, tracer in enumerate(tracers):
                if i <= last:
                    ma, ea = getIsotopeMass(tracer.isotopeA)
                    mb, eb = getIsotopeMass(tracer.isotopeB)
                    if tracer.name == "":
                        severe = True
                        QtGui.QMessageBox.question(self, "MetExtract",
                                                   "Error in tracer number %d \nYou did not specify a name. \nPlease provide a name or delete the row" % (
                                                       i + 1), QtGui.QMessageBox.Ok)
                    elif len(ea) == 0 or ea != eb:
                        ok = False
                        QtGui.QMessageBox.question(self, "MetExtract",
                                                   "Error in tracer %s \nYou cannot use two different elements for the labelling process. \nPlease enter isotopes of the same element" % tracer.name,
                                                   QtGui.QMessageBox.Ok)
                    for j, tr in enumerate(tracers):
                        if last >= j > i:
                            if (tracer.monoisotopicRatio >= tr.monoisotopicRatio and (
                                        tracer.monoisotopicRatio - tracer.maxRelNegBias) <= (
                                        tr.monoisotopicRatio + tr.maxRelPosBias)) or (
                                            tracer.monoisotopicRatio <= tr.monoisotopicRatio and (
                                                tracer.monoisotopicRatio + tracer.maxRelNegBias) >= (
                                                tr.monoisotopicRatio - tr.maxRelPosBias)):
                                ok = False
                                QtGui.QMessageBox.question(self, "MetExtract",
                                                           "Error in tracers %s and %s \nAllowed intensities overlap. \nPlease use non overlapping intensities" % (
                                                               tracer.name, tr.name), QtGui.QMessageBox.Ok)
                            if i < j and tracer.name == tr.name:
                                ok = False
                                QtGui.QMessageBox.question(self, "MetExtract",
                                                           "Error in tracers %s and %s (Rows: %d and %d)\nBoth have the same name" % (
                                                               tracer.name, tr.name, i, j), QtGui.QMessageBox.Ok)
            if severe:
                QtGui.QMessageBox.information(self, "MetExtract",
                                              "Settings are invalid. Please correct them or discard the changes",
                                              QtGui.QMessageBox.Ok)
                return False

            if not ok:
                if QtGui.QMessageBox.question(self, "MetExtract",
                                              "Problems were encountered. Continue with non optimal settings?",
                                              QtGui.QMessageBox.Yes | QtGui.QMessageBox.No) == QtGui.QMessageBox.Yes:
                    ok = True

        return ok and not severe

    def getOpenDir(self):
        return self.initDir

    def dialogCan(self):
        self.reject()

    def dialogFin(self):
        if self.checkConfiguration():
            self.accept()

    def executeDialog(self):
        x = self.exec_()

        d = []
        for tr in self.data:
            if tr.tracerType in ["user", "calc"]:
                d.append(tr)

        self.data = d

        return x


if __name__ == "__main__":
    import sys

    app = QtGui.QApplication(sys.argv)
    Dialog = tracerEdit()

    Dialog.show()
    x = app.exec_()

    print "final tracers:", len(Dialog.getTracers())
    for tracer in Dialog.getTracers():
        print tracer

    sys.exit(x)
    