# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\TracerEditor.ui'
#
# Created: Wed Dec 21 15:58:27 2016
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
        Dialog.resize(1236, 552)
        Dialog.setStyleSheet(_fromUtf8(""))
        self.gridLayout = QtGui.QGridLayout(Dialog)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout_2 = QtGui.QGridLayout()
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.line_4 = QtGui.QFrame(Dialog)
        self.line_4.setFrameShape(QtGui.QFrame.VLine)
        self.line_4.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout_2.addWidget(self.line_4, 3, 2, 1, 1)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.label_2 = QtGui.QLabel(Dialog)
        self.label_2.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout_2.addWidget(self.label_2)
        spacerItem = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        self.verticalLayout_2.addItem(spacerItem)
        self.gridLayout_2.addLayout(self.verticalLayout_2, 3, 0, 1, 1)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.title_tracer = QtGui.QLabel(Dialog)
        self.title_tracer.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.title_tracer.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.title_tracer.setObjectName(_fromUtf8("title_tracer"))
        self.verticalLayout.addWidget(self.title_tracer)
        spacerItem1 = QtGui.QSpacerItem(0, 15, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout.addItem(spacerItem1)
        self.label_3 = QtGui.QLabel(Dialog)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setMaximumSize(QtCore.QSize(200, 16777215))
        self.label_3.setStyleSheet(_fromUtf8("color: rgb(90, 90, 90);\n"
"font: 7pt;"))
        self.label_3.setAlignment(QtCore.Qt.AlignJustify|QtCore.Qt.AlignVCenter)
        self.label_3.setWordWrap(True)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout.addWidget(self.label_3)
        spacerItem2 = QtGui.QSpacerItem(200, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem2)
        self.gridLayout_2.addLayout(self.verticalLayout, 1, 0, 1, 1)
        self.tracers = QtGui.QTableView(Dialog)
        self.tracers.setObjectName(_fromUtf8("tracers"))
        self.gridLayout_2.addWidget(self.tracers, 1, 4, 1, 1)
        spacerItem3 = QtGui.QSpacerItem(3, 0, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem3, 1, 1, 1, 1)
        spacerItem4 = QtGui.QSpacerItem(3, 0, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem4, 1, 3, 1, 1)
        self.line_3 = QtGui.QFrame(Dialog)
        self.line_3.setFrameShape(QtGui.QFrame.HLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.gridLayout_2.addWidget(self.line_3, 2, 0, 1, 5)
        self.line = QtGui.QFrame(Dialog)
        self.line.setFrameShape(QtGui.QFrame.VLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout_2.addWidget(self.line, 1, 2, 1, 1)
        self.line_2 = QtGui.QFrame(Dialog)
        self.line_2.setFrameShape(QtGui.QFrame.HLine)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.gridLayout_2.addWidget(self.line_2, 0, 0, 1, 5)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.loadConfiguration = QtGui.QPushButton(Dialog)
        self.loadConfiguration.setObjectName(_fromUtf8("loadConfiguration"))
        self.horizontalLayout.addWidget(self.loadConfiguration)
        self.saveConfiguration = QtGui.QPushButton(Dialog)
        self.saveConfiguration.setObjectName(_fromUtf8("saveConfiguration"))
        self.horizontalLayout.addWidget(self.saveConfiguration)
        self.line_6 = QtGui.QFrame(Dialog)
        self.line_6.setFrameShape(QtGui.QFrame.VLine)
        self.line_6.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_6.setObjectName(_fromUtf8("line_6"))
        self.horizontalLayout.addWidget(self.line_6)
        self.checkRatios = QtGui.QCheckBox(Dialog)
        self.checkRatios.setChecked(True)
        self.checkRatios.setObjectName(_fromUtf8("checkRatios"))
        self.horizontalLayout.addWidget(self.checkRatios)
        self.calculateMultiTracerConjugates = QtGui.QPushButton(Dialog)
        self.calculateMultiTracerConjugates.setObjectName(_fromUtf8("calculateMultiTracerConjugates"))
        self.horizontalLayout.addWidget(self.calculateMultiTracerConjugates)
        self.line_5 = QtGui.QFrame(Dialog)
        self.line_5.setFrameShape(QtGui.QFrame.VLine)
        self.line_5.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_5.setObjectName(_fromUtf8("line_5"))
        self.horizontalLayout.addWidget(self.line_5)
        self.showMSScan = QtGui.QPushButton(Dialog)
        self.showMSScan.setObjectName(_fromUtf8("showMSScan"))
        self.horizontalLayout.addWidget(self.showMSScan)
        self.normaliseRatio = QtGui.QCheckBox(Dialog)
        self.normaliseRatio.setObjectName(_fromUtf8("normaliseRatio"))
        self.horizontalLayout.addWidget(self.normaliseRatio)
        spacerItem5 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem5)
        self.discardButton = QtGui.QPushButton(Dialog)
        self.discardButton.setObjectName(_fromUtf8("discardButton"))
        self.horizontalLayout.addWidget(self.discardButton)
        self.acceptButton = QtGui.QPushButton(Dialog)
        self.acceptButton.setObjectName(_fromUtf8("acceptButton"))
        self.horizontalLayout.addWidget(self.acceptButton)
        self.gridLayout_2.addLayout(self.horizontalLayout, 3, 4, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_2, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label_2.setText(_translate("Dialog", "Actions", None))
        self.title_tracer.setText(_translate("Dialog", "Tracers", None))
        self.label_3.setText(_translate("Dialog", "<html><head/><body><p>Please specify the tracers added to the experiment. Each row specifies a tracer</p><p><br/></p></body></html>", None))
        self.loadConfiguration.setText(_translate("Dialog", "Load", None))
        self.saveConfiguration.setText(_translate("Dialog", "Save", None))
        self.checkRatios.setText(_translate("Dialog", "Check ratios", None))
        self.calculateMultiTracerConjugates.setText(_translate("Dialog", "calculate multi tracer conjugates", None))
        self.showMSScan.setText(_translate("Dialog", "Show MS scan", None))
        self.normaliseRatio.setText(_translate("Dialog", "Normalise ratios", None))
        self.discardButton.setText(_translate("Dialog", "Discard", None))
        self.acceptButton.setText(_translate("Dialog", "Accept", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

