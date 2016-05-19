# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\heteroAtomEditor.ui'
#
# Created: Wed May 04 10:21:24 2016
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
        Dialog.resize(765, 474)
        Dialog.setStyleSheet(_fromUtf8(""))
        self.gridLayout = QtGui.QGridLayout(Dialog)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout_4 = QtGui.QGridLayout()
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.label_2 = QtGui.QLabel(Dialog)
        self.label_2.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout_2.addWidget(self.label_2)
        spacerItem = QtGui.QSpacerItem(0, 5, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout_2.addItem(spacerItem)
        self.help_2 = QtGui.QLabel(Dialog)
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
        self.verticalLayout_2.addWidget(self.help_2)
        spacerItem1 = QtGui.QSpacerItem(200, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem1)
        self.gridLayout_4.addLayout(self.verticalLayout_2, 2, 0, 1, 1)
        self.line = QtGui.QFrame(Dialog)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout_4.addWidget(self.line, 0, 0, 1, 3)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.loadDefaults = QtGui.QPushButton(Dialog)
        self.loadDefaults.setObjectName(_fromUtf8("loadDefaults"))
        self.horizontalLayout.addWidget(self.loadDefaults)
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem2)
        self.discardButton = QtGui.QPushButton(Dialog)
        self.discardButton.setObjectName(_fromUtf8("discardButton"))
        self.horizontalLayout.addWidget(self.discardButton)
        self.acceptButton = QtGui.QPushButton(Dialog)
        self.acceptButton.setObjectName(_fromUtf8("acceptButton"))
        self.horizontalLayout.addWidget(self.acceptButton)
        self.gridLayout_4.addLayout(self.horizontalLayout, 4, 2, 1, 1)
        self.heteroAtoms = QtGui.QTableView(Dialog)
        self.heteroAtoms.setSortingEnabled(True)
        self.heteroAtoms.setObjectName(_fromUtf8("heteroAtoms"))
        self.gridLayout_4.addWidget(self.heteroAtoms, 2, 2, 1, 1)
        self.line_4 = QtGui.QFrame(Dialog)
        self.line_4.setFrameShape(QtGui.QFrame.HLine)
        self.line_4.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout_4.addWidget(self.line_4, 3, 0, 1, 3)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.label_3 = QtGui.QLabel(Dialog)
        self.label_3.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_3.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout_3.addWidget(self.label_3)
        spacerItem3 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.verticalLayout_3.addItem(spacerItem3)
        self.gridLayout_4.addLayout(self.verticalLayout_3, 4, 0, 1, 1)
        self.line_5 = QtGui.QFrame(Dialog)
        self.line_5.setFrameShape(QtGui.QFrame.VLine)
        self.line_5.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_5.setObjectName(_fromUtf8("line_5"))
        self.gridLayout_4.addWidget(self.line_5, 4, 1, 1, 1)
        self.line_6 = QtGui.QFrame(Dialog)
        self.line_6.setFrameShape(QtGui.QFrame.VLine)
        self.line_6.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_6.setObjectName(_fromUtf8("line_6"))
        self.gridLayout_4.addWidget(self.line_6, 2, 1, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_4, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Relationship configuration", None))
        self.label_2.setText(_translate("Dialog", "Heteroatoms (elements)", None))
        self.help_2.setText(_translate("Dialog", "<html><head/><body><p>Please specify the elements you want to use for the targeted search for elements other then the labeling element. Each row represents one element</p><p>The first column specifies the chemical symbol of the respective isotope. Use &lt;protonNumber&gt;&lt;ElementSymbol&gt; (e.g. <sup>34</sup>S)<br/>The second column specifies the mass offset to the most abundant isotope of the elment (e.g.  <sup>34</sup>S and  <sup>32</sup>S have a mass difference of 1.9958)<br/>The third column specifies the expected intensity of this isotope in respect to the most abundant isotope of the respective element (e.g.  <sup>34</sup>S occours in nature with a probability of 4.21% and  <sup>32</sup>S with a probability of 95.02%. Therefore the calculated ratio for a molecule with a singe sulphur atom is 4.43%)<br/>The fourth and fifth column specify the minimal and maximal number this element can occour in non-targetd group annotation. </p></body></html>", None))
        self.loadDefaults.setText(_translate("Dialog", "Load defaults", None))
        self.discardButton.setText(_translate("Dialog", "Discard", None))
        self.acceptButton.setText(_translate("Dialog", "Accept", None))
        self.label_3.setText(_translate("Dialog", "Actions", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())

