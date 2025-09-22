from __future__ import print_function, division, absolute_import
import sys
import os
import pickle
import base64
from copy import deepcopy

from PySide6 import QtCore, QtGui, QtWidgets

from ..formulaTools import formulaTools

from .adductsEditor import Ui_Dialog


class ConfiguredAdduct:
    def __init__(
        self,
        name="",
        mzoffset=1.99705,
        charge=1,
        polarity="+",
        mCount=1,
        entryType="user",
    ):
        self.name = name
        self.mzoffset = mzoffset
        self.charge = charge
        self.polarity = polarity
        self.mCount = mCount
        self.entryType = entryType

    def __str__(self):
        return "ConfiguredAdduct (%s %.5f %s %d)" % (
            self.name,
            self.mzoffset,
            self.polarity,
            self.charge,
        )


class ConfiguredElement:
    def __init__(
        self,
        name="",
        weight=0.0,
        numberValenzElectrons=1,
        minCount=1,
        maxCount=1,
        entryType="user",
    ):
        self.name = name
        self.weight = weight
        self.numberValenzElectrons = numberValenzElectrons
        self.minCount = minCount
        self.maxCount = maxCount
        self.entryType = entryType

    def __str__(self):
        return "ConfiguredElement (%s %.5f %d %d %d)" % (
            self.name,
            self.weight,
            self.numberValenzElectrons,
            self.minCount,
            self.maxCount,
        )


defaultAdducts = [
    ConfiguredAdduct(
        name="[M+H]+",
        mzoffset=1.007276,
        charge=1,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+NH4]+",
        mzoffset=18.033823,
        charge=1,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+Na]+",
        mzoffset=22.989218,
        charge=1,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+CH3OH+H]+",
        mzoffset=33.033489,
        charge=1,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+K]+",
        mzoffset=38.963158,
        charge=1,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+ACN+H]+",
        mzoffset=42.033823,
        charge=1,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+2Na-H]+",
        mzoffset=44.971160,
        charge=1,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+2K-H]+",
        mzoffset=76.919040,
        charge=1,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+CH3FeN]+",
        mzoffset=84.96094,
        charge=1,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+2H]++",
        mzoffset=1.007276,
        charge=2,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+H+NH4]++",
        mzoffset=9.520550,
        charge=2,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+H+Na]++",
        mzoffset=11.998247,
        charge=2,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+H+K]++",
        mzoffset=19.985217,
        charge=2,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+2Na]++",
        mzoffset=22.989218,
        charge=2,
        polarity="+",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[2M+H]+",
        mzoffset=1.007276,
        charge=1,
        polarity="+",
        mCount=2,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[2M+NH4]+",
        mzoffset=18.033823,
        charge=1,
        polarity="+",
        mCount=2,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[2M+Na]+",
        mzoffset=22.989218,
        charge=1,
        polarity="+",
        mCount=2,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[2M+K]+",
        mzoffset=38.963158,
        charge=1,
        polarity="+",
        mCount=2,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M-H2O-H]-",
        mzoffset=-19.01839,
        charge=1,
        polarity="-",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M-H]-",
        mzoffset=-1.007276,
        charge=1,
        polarity="-",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+Na-2H]-",
        mzoffset=20.974666,
        charge=1,
        polarity="-",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+Cl]-",
        mzoffset=34.969402,
        charge=1,
        polarity="-",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+K-2H]-",
        mzoffset=36.948606,
        charge=1,
        polarity="-",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[M+Br]-",
        mzoffset=78.918885,
        charge=1,
        polarity="-",
        mCount=1,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[2M-H]-",
        mzoffset=-1.007276,
        charge=1,
        polarity="-",
        mCount=2,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[2M+Fa-H]-",
        mzoffset=44.998201,
        charge=1,
        polarity="-",
        mCount=2,
        entryType="user",
    ),
    ConfiguredAdduct(
        name="[2M+Hac-H]-",
        mzoffset=59.013851,
        charge=1,
        polarity="-",
        mCount=2,
        entryType="user",
    ),
]

defaultNeutralLosses = [
    ConfiguredElement(name="C", weight=12.0, numberValenzElectrons=4, minCount=0, maxCount=3),
    ConfiguredElement(name="H", weight=1.00783, numberValenzElectrons=1, minCount=0, maxCount=30),
    ConfiguredElement(name="O", weight=15.99491, numberValenzElectrons=6, minCount=0, maxCount=20),
    ConfiguredElement(name="N", weight=14.00307, numberValenzElectrons=5, minCount=0, maxCount=2),
    ConfiguredElement(name="P", weight=30.97376, numberValenzElectrons=5, minCount=0, maxCount=2),
    ConfiguredElement(name="S", weight=31.97207, numberValenzElectrons=6, minCount=0, maxCount=2),
    ConfiguredElement(name="Cl", weight=34.968852, numberValenzElectrons=7, minCount=0, maxCount=1),
]


class AdductsTableModel(QtCore.QAbstractTableModel):
    def __init__(self, datain, headerdata, parent=None, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)
        self.arraydata = datain
        self.headerdata = headerdata

    def rowCount(self, parent):
        return len(self.arraydata)

    def columnCount(self, parent):
        return 5

    def data(self, index, role):
        if not index.isValid() or index.row() == (len(self.arraydata) - 1):
            return None
        elif role != QtCore.Qt.DisplayRole and role != QtCore.Qt.EditRole:
            return None
        if index.column() == 0:
            return self.arraydata[index.row()].name
        if index.column() == 1:
            return "%.5f" % self.arraydata[index.row()].mzoffset
        if index.column() == 2:
            return self.arraydata[index.row()].polarity
        if index.column() == 3:
            return "%d" % self.arraydata[index.row()].charge
        if index.column() == 4:
            return "%d" % self.arraydata[index.row()].mCount

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self.headerdata[col]
        return None

    def insertRows(self, row, count, parent=QtCore.QModelIndex()):
        self.beginInsertRows(parent, row, row + count - 1)
        x = ConfiguredAdduct(entryType="empty")
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
            self.arraydata[-1].entryType = "user"
            self.insertRows(len(self.arraydata) + 1, 1)

        if index.column() == 1:
            try:
                f, ok = value.toDouble()
                if not ok:
                    return False
                self.arraydata[index.row()].mzoffset = f
                return True
            except:
                return False
        elif index.column() == 0:
            self.arraydata[index.row()].name = str(value.toString())
            return True
        elif index.column() == 2:
            val = str(value.toString())
            if val == "+" or val == "-":
                self.arraydata[index.row()].polarity = str(value.toString())
                return True
            else:
                return False
        elif index.column() == 3:
            try:
                f, ok = value.toInt()
                if not ok:
                    return False
                self.arraydata[index.row()].charge = f
                return True
            except:
                return False
        elif index.column() == 4:
            try:
                f, ok = value.toInt()
                if not ok:
                    return False
                self.arraydata[index.row()].mCount = f
                return True
            except:
                return False
        return False

    def flags(self, index):
        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable


class ElementsTableModel(QtCore.QAbstractTableModel):
    def __init__(self, datain, headerdata, parent=None, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)
        self.arraydata = datain
        self.headerdata = headerdata

        self.ft = formulaTools()

    def rowCount(self, parent):
        return len(self.arraydata)

    def columnCount(self, parent):
        return 5

    def data(self, index, role):
        if not index.isValid() or index.row() == (len(self.arraydata) - 1):
            return None
        elif role != QtCore.Qt.DisplayRole and role != QtCore.Qt.EditRole:
            return None
        if index.column() == 0:
            return self.arraydata[index.row()].name
        if index.column() == 1:
            return "%.5f" % self.arraydata[index.row()].weight
        if index.column() == 2:
            return "%d" % self.arraydata[index.row()].numberValenzElectrons
        if index.column() == 3:
            return "%d" % self.arraydata[index.row()].minCount
        if index.column() == 4:
            return "%d" % self.arraydata[index.row()].maxCount
        return None

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self.headerdata[col]
        return None

    def insertRows(self, row, count, parent=QtCore.QModelIndex()):
        self.beginInsertRows(parent, row, row + count - 1)
        x = ConfiguredElement(entryType="empty")
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
            self.arraydata[-1].entryType = "user"
            self.insertRows(len(self.arraydata) + 1, 1)

        if index.column() == 1:
            try:
                f, ok = value.toDouble()
                if not ok:
                    return False
                self.arraydata[index.row()].weight = f
                return True
            except:
                return False
        elif index.column() in [2, 3, 4]:
            try:
                f, ok = value.toInt()
                if not ok:
                    return False
                if index.column() == 2:
                    self.arraydata[index.row()].numberValenzElectrons = f
                if index.column() == 3:
                    self.arraydata[index.row()].minCount = f
                if index.column() == 4:
                    self.arraydata[index.row()].maxCount = f
                return True
            except:
                return False
        elif index.column() == 0:
            elem = str(value.toString())
            if elem in self.ft.elemDetails.keys():
                self.arraydata[index.row()].weight = self.ft.elemDetails[elem][3]
            self.arraydata[index.row()].name = elem

            return True
        return False

    def flags(self, index):
        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable


def saveANFile(filename, a, n):
    conf = QtCore.QSettings(filename, QtCore.QSettings.IniFormat)
    conf.setValue("adducts", base64.b64encode(pickle.dumps(a)))
    conf.setValue("neutralloss", base64.b64encode(pickle.dumps(n)))


def loadANFile(filename):
    conf = QtCore.QSettings(filename, QtCore.QSettings.IniFormat)
    a = pickle.loads(base64.b64decode(str(conf.value("adducts").toString())))
    n = pickle.loads(base64.b64decode(str(conf.value("neutralloss").toString())))
    return a, n


class adductsEdit(QtWidgets.QDialog, Ui_Dialog):
    def __init__(
        self,
        parent=None,
        initDir=None,
        adds=None,
        nls=None,
        showAdductsConfiguration=True,
        showRelationshipConfiguratio=True,
    ):
        if adds is None:
            adds = defaultAdducts
        if nls is None:
            nls = defaultNeutralLosses
        QtWidgets.QDialog.__init__(self, parent)
        self.setWindowTitle("Adducts editor")
        self.setupUi(self)

        if not (showAdductsConfiguration):
            # remove adducts controls from dialog
            self.adductsFrame.setVisible(False)

        if not (showRelationshipConfiguratio):
            # remove relationship controls from dialog
            self.neutralLossFrame.setVisible(False)

        self.adds = deepcopy(adds)
        self.adds.append(ConfiguredAdduct(entryType="empty"))
        self.addsModel = AdductsTableModel(self.adds, ["Adduct", "MZ offset", "Polarity", "Charge", "M count"])
        self.adducts.setModel(self.addsModel)

        self.nls = deepcopy(nls)
        self.nls.append(ConfiguredElement(entryType="empty"))
        self.nlModel = ElementsTableModel(
            self.nls,
            ["Element", "Weight", "Valenz electrons", "Min count", "Max count"],
        )
        self.neutralLoss.setModel(self.nlModel)

        self.initDir = "."
        if "USERPROFILE" in os.environ:
            self.initDir = os.getenv("USERPROFILE")
        elif "HOME" in os.environ:
            self.initDir = os.getenv("HOME")

        self.saveConfiguration.clicked.connect(self.save)
        self.loadConfiguration.clicked.connect(self.load)
        # self.saveConfiguration.setVisible(False)
        # self.loadConfiguration.setVisible(False)

        self.adducts.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.adducts.customContextMenuRequested.connect(self.showPopupA)
        self.neutralLoss.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.neutralLoss.customContextMenuRequested.connect(self.showPopupN)

        self.discardButton.clicked.connect(self.dialogCan)
        self.acceptButton.clicked.connect(self.dialogFin)

        self.loadDefaults.clicked.connect(self.loadDefaultAdductsAndNeutralLosses)

        self.acceptButton.setFocus()

    def showPopupA(self, position):
        menu = QtWidgets.QMenu()
        quitAction = menu.addAction("Delete")
        action = menu.exec_(self.adducts.mapToGlobal(position))
        if action == quitAction:
            for x in self.adducts.selectedIndexes():
                self.addsModel.removeRows(x.row(), 1)

    def showPopupN(self, position):
        menu = QtWidgets.QMenu()
        quitAction = menu.addAction("Delete")
        action = menu.exec_(self.neutralLoss.mapToGlobal(position))
        if action == quitAction:
            for x in self.neutralLoss.selectedIndexes():
                self.nlModel.removeRows(x.row(), 1)

    def load(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Select configuration file",
            self.initDir,
            "confif file (*.conf);;All files (*.*)",
        )
        if len(filename) > 0:
            a, n = loadANFile(filename)
            self.adds = a
            self.adds.append(ConfiguredAdduct(entryType="empty"))
            self.addsModel = AdductsTableModel(self.adds, ["Adduct", "MZ offset", "Polarity", "Charge"])
            self.adducts.setModel(self.addsModel)
            self.nls = n
            self.nls.append(ConfiguredElement(entryType="empty"))
            self.nlModel = ElementsTableModel(
                self.nls,
                ["Element", "Weight", "Valenz electrons", "Min count", "Max count"],
            )
            self.neutralLoss.setModel(self.nlModel)

    def save(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(
            caption="Select configuration file",
            filter="config file (*.conf);;All files(*.*)",
        )
        if len(filename) > 0:
            saveANFile(filename, self.getAdducts(), self.getNeutralLosses())

    def loadDefaultAdductsAndNeutralLosses(self):
        while self.addsModel.rowCount(self) > 1:
            self.addsModel.removeRows(0, 1)
        for i in range(len(defaultAdducts)):
            add = deepcopy(defaultAdducts[i])
            self.addsModel.insertRows(i, 1)
            self.adds[i] = add

        while self.nlModel.rowCount(self) > 1:
            self.nlModel.removeRows(0, 1)
        for i in range(len(defaultNeutralLosses)):
            nl = deepcopy(defaultNeutralLosses[i])
            self.nlModel.insertRows(i, 1)
            self.nls[i] = nl

    def getAdducts(self):
        return [add for add in self.adds if add.entryType != "empty"]

    def getNeutralLosses(self):
        return [nl for nl in self.nls if nl.entryType != "empty"]

    def dialogCan(self):
        self.reject()

    def dialogFin(self):
        self.accept()

    def executeDialog(self):
        x = self.exec()
        self.adds.pop(len(self.adds) - 1)
        self.nls.pop(len(self.nls) - 1)
        return x


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    Dialog = adductsEdit()

    Dialog.show()
    x = app.exec()

    nls = Dialog.getNeutralLosses()
    for nl in nls:
        print(nl)

    adds = Dialog.getAdducts()
    for add in adds:
        print(add)

    sys.exit(app.exec())
