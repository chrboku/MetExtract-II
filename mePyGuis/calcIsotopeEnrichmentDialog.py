# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\calcIsotopeEnrichmentDialog.ui'
#
# Created: Tue Jan 03 12:31:01 2017
#      by: PyQt4 UI code generator 4.10
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(1195, 1109)
        Dialog.setStyleSheet(_fromUtf8("background-color: white;"))
        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setGeometry(QtCore.QRect(950, 1060, 221, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.label = QtGui.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(19, 25, 1161, 501))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setStyleSheet(_fromUtf8("background-image: url(:/EnrichmentDialog/resources/EnrichmentDialog.png)"))
        self.label.setText(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setGeometry(QtCore.QRect(589, 85, 41, 21))
        self.label_3.setStyleSheet(_fromUtf8("font-size: 16pt;"))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.label_4 = QtGui.QLabel(Dialog)
        self.label_4.setGeometry(QtCore.QRect(90, 10, 21, 21))
        self.label_4.setStyleSheet(_fromUtf8("font-size: 16pt;"))
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.label_5 = QtGui.QLabel(Dialog)
        self.label_5.setGeometry(QtCore.QRect(130, 245, 61, 21))
        self.label_5.setStyleSheet(_fromUtf8("font-size: 16pt;"))
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.label_6 = QtGui.QLabel(Dialog)
        self.label_6.setGeometry(QtCore.QRect(192, 405, 61, 21))
        self.label_6.setStyleSheet(_fromUtf8("font-size: 16pt;"))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.label_7 = QtGui.QLabel(Dialog)
        self.label_7.setGeometry(QtCore.QRect(1088, 12, 31, 21))
        self.label_7.setStyleSheet(_fromUtf8("font-size: 16pt;"))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.label_8 = QtGui.QLabel(Dialog)
        self.label_8.setGeometry(QtCore.QRect(1010, 320, 61, 21))
        self.label_8.setStyleSheet(_fromUtf8("font-size: 16pt;"))
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.label_9 = QtGui.QLabel(Dialog)
        self.label_9.setGeometry(QtCore.QRect(951, 464, 61, 21))
        self.label_9.setStyleSheet(_fromUtf8("font-size: 16pt;"))
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.tableWidget = QtGui.QTableWidget(Dialog)
        self.tableWidget.setGeometry(QtCore.QRect(20, 530, 1151, 521))
        self.tableWidget.setObjectName(_fromUtf8("tableWidget"))
        self.tableWidget.setColumnCount(7)
        self.tableWidget.setRowCount(50)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(1, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(2, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(3, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(4, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(5, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(6, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(7, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(8, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(9, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(10, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(11, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(12, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(13, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(14, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(15, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(16, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(17, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(18, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(19, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(20, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(21, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(22, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(23, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(24, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(25, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(26, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(27, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(28, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(29, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(30, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(31, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(32, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(33, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(34, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(35, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(36, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(37, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(38, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(39, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(40, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(41, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(42, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(43, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(44, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(45, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(46, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(47, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(48, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(49, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(3, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(4, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(5, item)
        item = QtGui.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(6, item)

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), Dialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label_3.setText(_translate("Dialog", "Cn", None))
        self.label_4.setText(_translate("Dialog", "M", None))
        self.label_5.setText(_translate("Dialog", "M+1", None))
        self.label_6.setText(_translate("Dialog", "M+2", None))
        self.label_7.setText(_translate("Dialog", "M\'", None))
        self.label_8.setText(_translate("Dialog", "M\'-1", None))
        self.label_9.setText(_translate("Dialog", "M\'-2", None))
        item = self.tableWidget.verticalHeaderItem(0)
        item.setText(_translate("Dialog", "#1", None))
        item = self.tableWidget.verticalHeaderItem(1)
        item.setText(_translate("Dialog", "#2", None))
        item = self.tableWidget.verticalHeaderItem(2)
        item.setText(_translate("Dialog", "#3", None))
        item = self.tableWidget.verticalHeaderItem(3)
        item.setText(_translate("Dialog", "#4", None))
        item = self.tableWidget.verticalHeaderItem(4)
        item.setText(_translate("Dialog", "#5", None))
        item = self.tableWidget.verticalHeaderItem(5)
        item.setText(_translate("Dialog", "#6", None))
        item = self.tableWidget.verticalHeaderItem(6)
        item.setText(_translate("Dialog", "#7", None))
        item = self.tableWidget.verticalHeaderItem(7)
        item.setText(_translate("Dialog", "#8", None))
        item = self.tableWidget.verticalHeaderItem(8)
        item.setText(_translate("Dialog", "#9", None))
        item = self.tableWidget.verticalHeaderItem(9)
        item.setText(_translate("Dialog", "#10", None))
        item = self.tableWidget.verticalHeaderItem(10)
        item.setText(_translate("Dialog", "#11", None))
        item = self.tableWidget.verticalHeaderItem(11)
        item.setText(_translate("Dialog", "#12", None))
        item = self.tableWidget.verticalHeaderItem(12)
        item.setText(_translate("Dialog", "#13", None))
        item = self.tableWidget.verticalHeaderItem(13)
        item.setText(_translate("Dialog", "#14", None))
        item = self.tableWidget.verticalHeaderItem(14)
        item.setText(_translate("Dialog", "#15", None))
        item = self.tableWidget.verticalHeaderItem(15)
        item.setText(_translate("Dialog", "#16", None))
        item = self.tableWidget.verticalHeaderItem(16)
        item.setText(_translate("Dialog", "#17", None))
        item = self.tableWidget.verticalHeaderItem(17)
        item.setText(_translate("Dialog", "#18", None))
        item = self.tableWidget.verticalHeaderItem(18)
        item.setText(_translate("Dialog", "#19", None))
        item = self.tableWidget.verticalHeaderItem(19)
        item.setText(_translate("Dialog", "#20", None))
        item = self.tableWidget.verticalHeaderItem(20)
        item.setText(_translate("Dialog", "#21", None))
        item = self.tableWidget.verticalHeaderItem(21)
        item.setText(_translate("Dialog", "#22", None))
        item = self.tableWidget.verticalHeaderItem(22)
        item.setText(_translate("Dialog", "#23", None))
        item = self.tableWidget.verticalHeaderItem(23)
        item.setText(_translate("Dialog", "#24", None))
        item = self.tableWidget.verticalHeaderItem(24)
        item.setText(_translate("Dialog", "#25", None))
        item = self.tableWidget.verticalHeaderItem(25)
        item.setText(_translate("Dialog", "#26", None))
        item = self.tableWidget.verticalHeaderItem(26)
        item.setText(_translate("Dialog", "#27", None))
        item = self.tableWidget.verticalHeaderItem(27)
        item.setText(_translate("Dialog", "#28", None))
        item = self.tableWidget.verticalHeaderItem(28)
        item.setText(_translate("Dialog", "#29", None))
        item = self.tableWidget.verticalHeaderItem(29)
        item.setText(_translate("Dialog", "#30", None))
        item = self.tableWidget.verticalHeaderItem(30)
        item.setText(_translate("Dialog", "#31", None))
        item = self.tableWidget.verticalHeaderItem(31)
        item.setText(_translate("Dialog", "#32", None))
        item = self.tableWidget.verticalHeaderItem(32)
        item.setText(_translate("Dialog", "#33", None))
        item = self.tableWidget.verticalHeaderItem(33)
        item.setText(_translate("Dialog", "#34", None))
        item = self.tableWidget.verticalHeaderItem(34)
        item.setText(_translate("Dialog", "#35", None))
        item = self.tableWidget.verticalHeaderItem(35)
        item.setText(_translate("Dialog", "#36", None))
        item = self.tableWidget.verticalHeaderItem(36)
        item.setText(_translate("Dialog", "#37", None))
        item = self.tableWidget.verticalHeaderItem(37)
        item.setText(_translate("Dialog", "#38", None))
        item = self.tableWidget.verticalHeaderItem(38)
        item.setText(_translate("Dialog", "#39", None))
        item = self.tableWidget.verticalHeaderItem(39)
        item.setText(_translate("Dialog", "#40", None))
        item = self.tableWidget.verticalHeaderItem(40)
        item.setText(_translate("Dialog", "#41", None))
        item = self.tableWidget.verticalHeaderItem(41)
        item.setText(_translate("Dialog", "#42", None))
        item = self.tableWidget.verticalHeaderItem(42)
        item.setText(_translate("Dialog", "#43", None))
        item = self.tableWidget.verticalHeaderItem(43)
        item.setText(_translate("Dialog", "#44", None))
        item = self.tableWidget.verticalHeaderItem(44)
        item.setText(_translate("Dialog", "#45", None))
        item = self.tableWidget.verticalHeaderItem(45)
        item.setText(_translate("Dialog", "#46", None))
        item = self.tableWidget.verticalHeaderItem(46)
        item.setText(_translate("Dialog", "#47", None))
        item = self.tableWidget.verticalHeaderItem(47)
        item.setText(_translate("Dialog", "#48", None))
        item = self.tableWidget.verticalHeaderItem(48)
        item.setText(_translate("Dialog", "#49", None))
        item = self.tableWidget.verticalHeaderItem(49)
        item.setText(_translate("Dialog", "#50", None))
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("Dialog", "Abundance M", None))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("Dialog", "Abundance M+1", None))
        item = self.tableWidget.horizontalHeaderItem(2)
        item.setText(_translate("Dialog", "Abundance M\'-1", None))
        item = self.tableWidget.horizontalHeaderItem(3)
        item.setText(_translate("Dialog", "Abundance M\'", None))
        item = self.tableWidget.horizontalHeaderItem(4)
        item.setText(_translate("Dialog", "Cn", None))
        item = self.tableWidget.horizontalHeaderItem(5)
        item.setText(_translate("Dialog", "12C-enrichment", None))
        item = self.tableWidget.horizontalHeaderItem(6)
        item.setText(_translate("Dialog", "13C-enrichment", None))

import resources_rc

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

