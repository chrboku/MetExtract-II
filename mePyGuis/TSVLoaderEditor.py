# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\TSVLoaderEditor.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
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
        Dialog.resize(1036, 266)
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
        self.titleText = QtGui.QLabel(Dialog)
        self.titleText.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.titleText.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.titleText.setObjectName(_fromUtf8("titleText"))
        self.verticalLayout.addWidget(self.titleText)
        spacerItem1 = QtGui.QSpacerItem(0, 15, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout.addItem(spacerItem1)
        self.descriptionText = QtGui.QLabel(Dialog)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.descriptionText.sizePolicy().hasHeightForWidth())
        self.descriptionText.setSizePolicy(sizePolicy)
        self.descriptionText.setMaximumSize(QtCore.QSize(200, 16777215))
        self.descriptionText.setStyleSheet(_fromUtf8("color: rgb(90, 90, 90);\n"
"font: 7pt;"))
        self.descriptionText.setAlignment(QtCore.Qt.AlignJustify|QtCore.Qt.AlignVCenter)
        self.descriptionText.setWordWrap(True)
        self.descriptionText.setObjectName(_fromUtf8("descriptionText"))
        self.verticalLayout.addWidget(self.descriptionText)
        spacerItem2 = QtGui.QSpacerItem(200, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem2)
        self.gridLayout_2.addLayout(self.verticalLayout, 1, 0, 1, 1)
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
        self.loadTable = QtGui.QPushButton(Dialog)
        self.loadTable.setObjectName(_fromUtf8("loadTable"))
        self.horizontalLayout.addWidget(self.loadTable)
        spacerItem5 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem5)
        self.discardButton = QtGui.QPushButton(Dialog)
        self.discardButton.setObjectName(_fromUtf8("discardButton"))
        self.horizontalLayout.addWidget(self.discardButton)
        self.acceptButton = QtGui.QPushButton(Dialog)
        self.acceptButton.setObjectName(_fromUtf8("acceptButton"))
        self.horizontalLayout.addWidget(self.acceptButton)
        self.gridLayout_2.addLayout(self.horizontalLayout, 3, 4, 1, 1)
        self.scrollArea = QtGui.QScrollArea(Dialog)
        self.scrollArea.setStyleSheet(_fromUtf8(""))
        self.scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName(_fromUtf8("scrollArea"))
        self.scrollAreaWidgetContents = QtGui.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 762, 195))
        self.scrollAreaWidgetContents.setObjectName(_fromUtf8("scrollAreaWidgetContents"))
        self.gridLayout_3 = QtGui.QGridLayout(self.scrollAreaWidgetContents)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        spacerItem6 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem6)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.l1 = QtGui.QLabel(self.scrollAreaWidgetContents)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.l1.sizePolicy().hasHeightForWidth())
        self.l1.setSizePolicy(sizePolicy)
        self.l1.setObjectName(_fromUtf8("l1"))
        self.horizontalLayout_5.addWidget(self.l1)
        self.c1 = QtGui.QComboBox(self.scrollAreaWidgetContents)
        self.c1.setObjectName(_fromUtf8("c1"))
        self.horizontalLayout_5.addWidget(self.c1)
        self.verticalLayout_3.addLayout(self.horizontalLayout_5)
        spacerItem7 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem7)
        self.gridLayout_3.addLayout(self.verticalLayout_3, 0, 0, 1, 1)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.gridLayout_2.addWidget(self.scrollArea, 1, 4, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_2, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label_2.setText(_translate("Dialog", "Actions", None))
        self.titleText.setText(_translate("Dialog", "some title", None))
        self.descriptionText.setText(_translate("Dialog", "<html><head/><body><p>some description</p></body></html>", None))
        self.loadTable.setText(_translate("Dialog", "Load table", None))
        self.discardButton.setText(_translate("Dialog", "Discard", None))
        self.acceptButton.setText(_translate("Dialog", "Accept", None))
        self.l1.setText(_translate("Dialog", "<html><head/><body><p><span style=\" text-decoration: underline;\">ID</span> map to</p></body></html>", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

