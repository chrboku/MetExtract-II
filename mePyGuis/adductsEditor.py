# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\adductsEditor.ui'
#
# Created: Tue Jul 19 15:19:54 2016
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
        Dialog.resize(819, 610)
        Dialog.setStyleSheet(_fromUtf8(""))
        self.gridLayout = QtGui.QGridLayout(Dialog)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout_4 = QtGui.QGridLayout()
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.adductsFrame = QtGui.QFrame(Dialog)
        self.adductsFrame.setFrameShape(QtGui.QFrame.StyledPanel)
        self.adductsFrame.setFrameShadow(QtGui.QFrame.Raised)
        self.adductsFrame.setObjectName(_fromUtf8("adductsFrame"))
        self.gridLayout_2 = QtGui.QGridLayout(self.adductsFrame)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.adducts = QtGui.QTableView(self.adductsFrame)
        self.adducts.setObjectName(_fromUtf8("adducts"))
        self.gridLayout_2.addWidget(self.adducts, 0, 4, 1, 1)
        spacerItem = QtGui.QSpacerItem(3, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem, 0, 3, 1, 1)
        self.adductsLine = QtGui.QFrame(self.adductsFrame)
        self.adductsLine.setFrameShape(QtGui.QFrame.VLine)
        self.adductsLine.setFrameShadow(QtGui.QFrame.Sunken)
        self.adductsLine.setObjectName(_fromUtf8("adductsLine"))
        self.gridLayout_2.addWidget(self.adductsLine, 0, 2, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(3, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem1, 0, 1, 1, 1)
        self.adductsHelp = QtGui.QVBoxLayout()
        self.adductsHelp.setObjectName(_fromUtf8("adductsHelp"))
        self.adductsHelpTitle = QtGui.QLabel(self.adductsFrame)
        self.adductsHelpTitle.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.adductsHelpTitle.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.adductsHelpTitle.setObjectName(_fromUtf8("adductsHelpTitle"))
        self.adductsHelp.addWidget(self.adductsHelpTitle)
        spacerItem2 = QtGui.QSpacerItem(0, 5, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.adductsHelp.addItem(spacerItem2)
        self.adductsHelpText = QtGui.QLabel(self.adductsFrame)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.adductsHelpText.sizePolicy().hasHeightForWidth())
        self.adductsHelpText.setSizePolicy(sizePolicy)
        self.adductsHelpText.setStyleSheet(_fromUtf8("color: rgb(90, 90, 90);\n"
"font: 7pt;"))
        self.adductsHelpText.setAlignment(QtCore.Qt.AlignJustify|QtCore.Qt.AlignVCenter)
        self.adductsHelpText.setWordWrap(True)
        self.adductsHelpText.setObjectName(_fromUtf8("adductsHelpText"))
        self.adductsHelp.addWidget(self.adductsHelpText)
        spacerItem3 = QtGui.QSpacerItem(200, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.adductsHelp.addItem(spacerItem3)
        self.gridLayout_2.addLayout(self.adductsHelp, 0, 0, 1, 1)
        self.gridLayout_4.addWidget(self.adductsFrame, 3, 0, 1, 5)
        self.line = QtGui.QFrame(Dialog)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout_4.addWidget(self.line, 0, 0, 1, 5)
        self.line_5 = QtGui.QFrame(Dialog)
        self.line_5.setFrameShape(QtGui.QFrame.VLine)
        self.line_5.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_5.setObjectName(_fromUtf8("line_5"))
        self.gridLayout_4.addWidget(self.line_5, 9, 2, 1, 1)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_3.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout_3.addWidget(self.label_3)
        spacerItem4 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.verticalLayout_3.addItem(spacerItem4)
        self.gridLayout_4.addLayout(self.verticalLayout_3, 9, 0, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.loadDefaults = QtGui.QPushButton(Dialog)
        self.loadDefaults.setObjectName(_fromUtf8("loadDefaults"))
        self.horizontalLayout.addWidget(self.loadDefaults)
        self.loadConfiguration = QtGui.QPushButton(Dialog)
        self.loadConfiguration.setObjectName(_fromUtf8("loadConfiguration"))
        self.horizontalLayout.addWidget(self.loadConfiguration)
        self.saveConfiguration = QtGui.QPushButton(Dialog)
        self.saveConfiguration.setObjectName(_fromUtf8("saveConfiguration"))
        self.horizontalLayout.addWidget(self.saveConfiguration)
        spacerItem5 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem5)
        self.discardButton = QtGui.QPushButton(Dialog)
        self.discardButton.setObjectName(_fromUtf8("discardButton"))
        self.horizontalLayout.addWidget(self.discardButton)
        self.acceptButton = QtGui.QPushButton(Dialog)
        self.acceptButton.setObjectName(_fromUtf8("acceptButton"))
        self.horizontalLayout.addWidget(self.acceptButton)
        self.gridLayout_4.addLayout(self.horizontalLayout, 9, 4, 1, 1)
        self.line_4 = QtGui.QFrame(Dialog)
        self.line_4.setFrameShape(QtGui.QFrame.HLine)
        self.line_4.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout_4.addWidget(self.line_4, 8, 0, 1, 5)
        self.line_3 = QtGui.QFrame(Dialog)
        self.line_3.setFrameShape(QtGui.QFrame.HLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.gridLayout_4.addWidget(self.line_3, 5, 0, 1, 5)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.gridLayout_4.addLayout(self.horizontalLayout_2, 4, 0, 1, 5)
        self.neutralLossFrame = QtGui.QFrame(Dialog)
        self.neutralLossFrame.setFrameShape(QtGui.QFrame.StyledPanel)
        self.neutralLossFrame.setFrameShadow(QtGui.QFrame.Raised)
        self.neutralLossFrame.setObjectName(_fromUtf8("neutralLossFrame"))
        self.gridLayout_3 = QtGui.QGridLayout(self.neutralLossFrame)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        spacerItem6 = QtGui.QSpacerItem(3, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem6, 0, 3, 1, 1)
        self.neutralLossLine = QtGui.QFrame(self.neutralLossFrame)
        self.neutralLossLine.setFrameShape(QtGui.QFrame.VLine)
        self.neutralLossLine.setFrameShadow(QtGui.QFrame.Sunken)
        self.neutralLossLine.setObjectName(_fromUtf8("neutralLossLine"))
        self.gridLayout_3.addWidget(self.neutralLossLine, 0, 2, 1, 1)
        self.neutralLoss = QtGui.QTableView(self.neutralLossFrame)
        self.neutralLoss.setObjectName(_fromUtf8("neutralLoss"))
        self.gridLayout_3.addWidget(self.neutralLoss, 0, 4, 1, 1)
        spacerItem7 = QtGui.QSpacerItem(3, 20, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem7, 0, 1, 1, 1)
        self.neutralLossHelp = QtGui.QVBoxLayout()
        self.neutralLossHelp.setObjectName(_fromUtf8("neutralLossHelp"))
        self.label_2 = QtGui.QLabel(self.neutralLossFrame)
        self.label_2.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.neutralLossHelp.addWidget(self.label_2)
        spacerItem8 = QtGui.QSpacerItem(0, 5, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.neutralLossHelp.addItem(spacerItem8)
        self.help_2 = QtGui.QLabel(self.neutralLossFrame)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.help_2.sizePolicy().hasHeightForWidth())
        self.help_2.setSizePolicy(sizePolicy)
        self.help_2.setMaximumSize(QtCore.QSize(200, 16777215))
        self.help_2.setStyleSheet(_fromUtf8("color: rgb(90, 90, 90);\n"
"font: 7pt;"))
        self.help_2.setAlignment(QtCore.Qt.AlignJustify|QtCore.Qt.AlignVCenter)
        self.help_2.setWordWrap(True)
        self.help_2.setObjectName(_fromUtf8("help_2"))
        self.neutralLossHelp.addWidget(self.help_2)
        spacerItem9 = QtGui.QSpacerItem(200, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.neutralLossHelp.addItem(spacerItem9)
        self.gridLayout_3.addLayout(self.neutralLossHelp, 0, 0, 1, 1)
        self.gridLayout_4.addWidget(self.neutralLossFrame, 6, 0, 1, 5)
        self.gridLayout.addLayout(self.gridLayout_4, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Relationship configuration", None))
        self.adductsHelpTitle.setText(_translate("Dialog", "Adducts", None))
        self.adductsHelpText.setText(_translate("Dialog", "<html><head/><body><p>Please specify the adducts you want to use during the non-targeted feature grouping. Each row in the table represents one adduct</p><p>The first column is the name of the adduct<br/>The second column specifies the adducts m/z difference to a neutral molecule<br/>The third column specifies in which ionisation mode the adduct may occour. Valid values are + and -</p></body></html>", None))
        self.label_3.setText(_translate("Dialog", "Actions", None))
        self.loadDefaults.setText(_translate("Dialog", "Load defaults", None))
        self.loadConfiguration.setText(_translate("Dialog", "Load", None))
        self.saveConfiguration.setText(_translate("Dialog", "Save", None))
        self.discardButton.setText(_translate("Dialog", "Discard", None))
        self.acceptButton.setText(_translate("Dialog", "Accept", None))
        self.label_2.setText(_translate("Dialog", "Neutral loss (elements)", None))
        self.help_2.setText(_translate("Dialog", "<html><head/><body><p>Please specify the elements you want to use for the non-targeted feature group annotation. Each row represents one element</p><p>The first column specifies the chemical symbol of the element<br/>The second column specifies the neutral weight of the most abundant isotope of the elment<br/>The third column specifies the number of valenz electron this element has<br/>The fourth and fifth column specify the minimal and maximal number this element can occour in non-targetd group annotation</p></body></html>", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

