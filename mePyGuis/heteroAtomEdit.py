import sys
import os
import pickle
import base64
from copy import deepcopy

from PyQt4 import QtGui, QtCore

from mePyGuis.heteroAtomEditor import Ui_Dialog

from formulaTools import formulaTools


class ConfiguredHeteroAtom():
    def __init__(self, name="34S", mzOffset=1.995796, relativeAbundance=.0443, minCount=0, maxCount=4,
                 entryType="user"):
        self.name = name
        self.mzOffset = mzOffset
        self.relativeAbundance = relativeAbundance
        self.minCount = minCount
        self.maxCount = maxCount
        self.entryType = entryType

    def __str__(self):
        return "ConfiguredHeteroAtom: (%s %.5f)" %( self.name, self.mzOffset)


defaultHeteroAtoms = [
    ConfiguredHeteroAtom(name='34S',  mzOffset=1.9957960000000021, relativeAbundance=0.044306461797516308, minCount=0, maxCount=4),
    ConfiguredHeteroAtom(name='41K',  mzOffset=1.9981170000000006, relativeAbundance=0.07221030042918454,  minCount=0, maxCount=1),
    ConfiguredHeteroAtom(name='37Cl', mzOffset=1.9970499999999944, relativeAbundance=0.31978355549689846,  minCount=0, maxCount=2)
]


class HeteroAtomsTableModel(QtCore.QAbstractTableModel):
    def __init__(self, datain, headerdata, parent=None, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)
        self.arraydata = datain
        self.headerdata = headerdata

        self.fT = formulaTools();

    def rowCount(self, parent):
        return len(self.arraydata)

    def columnCount(self, parent):
        return 5

    def data(self, index, role):
        if not index.isValid():
            return QtCore.QVariant()
        elif role != QtCore.Qt.DisplayRole and role != QtCore.Qt.EditRole:
            return QtCore.QVariant()
        if index.row() == (len(self.arraydata) - 1):
            return QtCore.QString("")
        if index.column() == 0:
            return QtCore.QString(self.arraydata[index.row()].name)
        if index.column() == 1:
            if isinstance(self.arraydata[index.row()].mzOffset, float):
                return QtCore.QString("%.5f" % self.arraydata[index.row()].mzOffset)
            return QtCore.QString(self.arraydata[index.row()].mzOffset)
        if index.column() == 2:
            a = self.arraydata[index.row()].relativeAbundance
            a = a * 100
            return QtCore.QVariant(a)
        if index.column() == 3:
            return QtCore.QVariant(self.arraydata[index.row()].minCount)
        if index.column() == 4:
            return QtCore.QVariant(self.arraydata[index.row()].maxCount)
        return QtCore.QVariant()

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return QtCore.QVariant(self.headerdata[col])
        return QtCore.QVariant()

    def insertRows(self, row, count, parent=QtCore.QModelIndex()):
        self.beginInsertRows(parent, row, row + count - 1)
        self.arraydata[-1].entryType = "user"
        x = ConfiguredHeteroAtom(entryType="empty")
        self.arraydata.append(x)
        self.endInsertRows()
        return True

    def removeRows(self, row, count, parent=QtCore.QModelIndex()):
        assert len(self.arraydata) > row and len(self.arraydata) > (row + count)
        self.beginRemoveRows(parent, row, row + count - 1)
        for i in range(count):
            self.arraydata.pop(row)
        self.endRemoveRows()
        return True

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        f = 0
        if index.row() == (len(self.arraydata) - 1) and not (value == ""):
            self.insertRows(len(self.arraydata) + 1, 1)

        if index.column() in [1, 2]:
            try:
                f, ok = value.toDouble()
                if index.column() == 2:
                    f = f / 100.
                if not ok:
                    return False
                if index.column() == 1:
                    self.arraydata[index.row()].mzOffset = f
                if index.column() == 2:
                    self.arraydata[index.row()].relativeAbundance = f
                return True
            except:
                return False
        elif index.column() in [3, 4]:
            try:
                f, ok = value.toInt()
                if not ok:
                    return False
                if index.column() == 3:
                    self.arraydata[index.row()].minCount = f
                if index.column() == 4:
                    self.arraydata[index.row()].maxCount = f
                return True
            except:
                return False
        elif index.column() == 0:
            val = str(value.toString())
            if val in self.fT.getIsotopes(minInt=0.0).keys():
                self.arraydata[index.row()].mzOffset = self.fT.getIsotopes(minInt=0.0)[val][0]
                self.arraydata[index.row()].relativeAbundance = self.fT.getIsotopes(minInt=0.0)[val][1]
                self.arraydata[index.row()].minCount = 0
                self.arraydata[index.row()].maxCount = 4

            self.arraydata[index.row()].name = val
            return True
        return False

    def flags(self, index):
        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable


class heteroAtomEdit(QtGui.QDialog, Ui_Dialog):
    def __init__(self, parent=None, initDir=None, heteroAtoms=None):
        if heteroAtoms is None:
            heteroAtoms = defaultHeteroAtoms
        QtGui.QDialog.__init__(self, parent)
        self.setWindowTitle("Hetero atom editor")
        self.setupUi(self)

        self.hAtoms = deepcopy(heteroAtoms)
        self.hAtoms.append(ConfiguredHeteroAtom(entryType="empty"))
        self.haModel = HeteroAtomsTableModel(self.hAtoms,
                                             ["Isotope", u"\u0394 m/z", "Rel. ab. [%]", "Min. atoms", "Max. atoms"])

        self.heteroAtoms.setModel(self.haModel)

        self.heteroAtoms.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.heteroAtoms.customContextMenuRequested.connect(self.showPopupHA)

        self.discardButton.clicked.connect(self.dialogCan)
        self.acceptButton.clicked.connect(self.dialogFin)

        self.loadDefaults.clicked.connect(self.loadDefaultHeteroAtoms)

        self.acceptButton.setFocus(True)

    def showPopupHA(self, position):
        menu = QtGui.QMenu()
        quitAction = menu.addAction("Delete")
        action = menu.exec_(self.heteroAtoms.mapToGlobal(position))
        if action == quitAction:
            for x in self.heteroAtoms.selectedIndexes():
                self.haModel.removeRows(x.row(), 1)

    def loadDefaultHeteroAtoms(self):
        while self.haModel.rowCount(self)>1:
            self.haModel.removeRows(0, 1)
        for i in range(len(defaultHeteroAtoms)):
            ha=deepcopy(defaultHeteroAtoms[i])
            self.haModel.insertRows(i, 1)
            self.hAtoms[i]=ha


    def dialogCan(self):
        self.reject()

    def dialogFin(self):
        self.accept()

    def getHeteroAtoms(self):
        return [ha for ha in self.hAtoms if ha.entryType != "empty"]

    def executeDialog(self):
        x = self.exec_()
        self.hAtoms.pop(len(self.hAtoms) - 1)
        return x


if __name__ == "__main__":
    import sys

    app = QtGui.QApplication(sys.argv)
    Dialog = heteroAtomEdit()

    Dialog.show()
    x = app.exec_()

    has = Dialog.getHeteroAtoms()

    import pickle
    import base64

    print base64.b64encode(pickle.dumps(has)), has
    for ha in has:
        print ha
    sys.exit(app.exec_())











