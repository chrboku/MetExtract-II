# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\groupEditor.ui'
#
# Created: Tue Aug 16 17:25:57 2016
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

class Ui_GroupEditor(object):
    def setupUi(self, GroupEditor):
        GroupEditor.setObjectName(_fromUtf8("GroupEditor"))
        GroupEditor.resize(598, 552)
        GroupEditor.setStyleSheet(_fromUtf8(""))
        self.gridLayout = QtGui.QGridLayout(GroupEditor)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout_2 = QtGui.QGridLayout()
        self.gridLayout_2.setMargin(9)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.label_3 = QtGui.QLabel(GroupEditor)
        self.label_3.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_3.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout_2.addWidget(self.label_3)
        spacerItem = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.verticalLayout_2.addItem(spacerItem)
        self.gridLayout_2.addLayout(self.verticalLayout_2, 7, 0, 1, 1)
        self.line_2 = QtGui.QFrame(GroupEditor)
        self.line_2.setFrameShape(QtGui.QFrame.HLine)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.gridLayout_2.addWidget(self.line_2, 2, 0, 1, 5)
        spacerItem1 = QtGui.QSpacerItem(3, 0, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem1, 1, 3, 1, 1)
        self.line_3 = QtGui.QFrame(GroupEditor)
        self.line_3.setFrameShape(QtGui.QFrame.VLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.gridLayout_2.addWidget(self.line_3, 1, 2, 1, 1)
        spacerItem2 = QtGui.QSpacerItem(3, 0, QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem2, 1, 1, 1, 1)
        self.groupName = QtGui.QLineEdit(GroupEditor)
        self.groupName.setAcceptDrops(False)
        self.groupName.setPlaceholderText(_fromUtf8(""))
        self.groupName.setObjectName(_fromUtf8("groupName"))
        self.gridLayout_2.addWidget(self.groupName, 1, 4, 1, 1)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        spacerItem3 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem3)
        self.dialogCanceled = QtGui.QPushButton(GroupEditor)
        self.dialogCanceled.setObjectName(_fromUtf8("dialogCanceled"))
        self.horizontalLayout_2.addWidget(self.dialogCanceled)
        self.dialogFinished = QtGui.QPushButton(GroupEditor)
        self.dialogFinished.setObjectName(_fromUtf8("dialogFinished"))
        self.horizontalLayout_2.addWidget(self.dialogFinished)
        self.gridLayout_2.addLayout(self.horizontalLayout_2, 9, 4, 1, 1)
        self.gridLayout_3 = QtGui.QGridLayout()
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        spacerItem4 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_3.addItem(spacerItem4, 0, 0, 1, 1)
        self.groupMinimumFound = QtGui.QSpinBox(GroupEditor)
        self.groupMinimumFound.setMinimum(1)
        self.groupMinimumFound.setMaximum(999)
        self.groupMinimumFound.setProperty("value", 1)
        self.groupMinimumFound.setObjectName(_fromUtf8("groupMinimumFound"))
        self.gridLayout_3.addWidget(self.groupMinimumFound, 0, 3, 1, 1)
        self.colorsComboBox = QtGui.QComboBox(GroupEditor)
        self.colorsComboBox.setObjectName(_fromUtf8("colorsComboBox"))
        self.gridLayout_3.addWidget(self.colorsComboBox, 3, 2, 1, 2)
        self.label_7 = QtGui.QLabel(GroupEditor)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout_3.addWidget(self.label_7, 3, 1, 1, 1)
        self.label_4 = QtGui.QLabel(GroupEditor)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout_3.addWidget(self.label_4, 0, 1, 1, 2)
        self.omitFeatures = QtGui.QCheckBox(GroupEditor)
        self.omitFeatures.setChecked(True)
        self.omitFeatures.setObjectName(_fromUtf8("omitFeatures"))
        self.gridLayout_3.addWidget(self.omitFeatures, 1, 1, 1, 3)
        self.useForMetaboliteGrouping = QtGui.QCheckBox(GroupEditor)
        self.useForMetaboliteGrouping.setChecked(True)
        self.useForMetaboliteGrouping.setObjectName(_fromUtf8("useForMetaboliteGrouping"))
        self.gridLayout_3.addWidget(self.useForMetaboliteGrouping, 2, 1, 1, 3)
        self.gridLayout_2.addLayout(self.gridLayout_3, 7, 4, 1, 1)
        self.line = QtGui.QFrame(GroupEditor)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout_2.addWidget(self.line, 0, 0, 1, 5)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_2 = QtGui.QLabel(GroupEditor)
        self.label_2.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout.addWidget(self.label_2)
        spacerItem5 = QtGui.QSpacerItem(0, 5, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.verticalLayout.addItem(spacerItem5)
        self.label_6 = QtGui.QLabel(GroupEditor)
        self.label_6.setStyleSheet(_fromUtf8("color: rgb(90, 90, 90);\n"
"font: 7pt;"))
        self.label_6.setAlignment(QtCore.Qt.AlignJustify|QtCore.Qt.AlignVCenter)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.verticalLayout.addWidget(self.label_6)
        spacerItem6 = QtGui.QSpacerItem(200, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem6)
        self.gridLayout_2.addLayout(self.verticalLayout, 4, 0, 2, 1)
        self.groupFiles = QtGui.QListWidget(GroupEditor)
        self.groupFiles.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.groupFiles.setObjectName(_fromUtf8("groupFiles"))
        self.gridLayout_2.addWidget(self.groupFiles, 4, 4, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem7 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem7)
        self.addFolder = QtGui.QPushButton(GroupEditor)
        self.addFolder.setObjectName(_fromUtf8("addFolder"))
        self.horizontalLayout.addWidget(self.addFolder)
        self.addFiles = QtGui.QPushButton(GroupEditor)
        self.addFiles.setObjectName(_fromUtf8("addFiles"))
        self.horizontalLayout.addWidget(self.addFiles)
        self.removeSelected = QtGui.QPushButton(GroupEditor)
        self.removeSelected.setObjectName(_fromUtf8("removeSelected"))
        self.horizontalLayout.addWidget(self.removeSelected)
        self.gridLayout_2.addLayout(self.horizontalLayout, 5, 4, 1, 1)
        self.line_5 = QtGui.QFrame(GroupEditor)
        self.line_5.setFrameShape(QtGui.QFrame.HLine)
        self.line_5.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_5.setObjectName(_fromUtf8("line_5"))
        self.gridLayout_2.addWidget(self.line_5, 6, 0, 1, 5)
        self.line_6 = QtGui.QFrame(GroupEditor)
        self.line_6.setFrameShape(QtGui.QFrame.HLine)
        self.line_6.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_6.setObjectName(_fromUtf8("line_6"))
        self.gridLayout_2.addWidget(self.line_6, 8, 0, 1, 5)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.label_5 = QtGui.QLabel(GroupEditor)
        self.label_5.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_5.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.verticalLayout_3.addWidget(self.label_5)
        spacerItem8 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.verticalLayout_3.addItem(spacerItem8)
        self.gridLayout_2.addLayout(self.verticalLayout_3, 9, 0, 1, 1)
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.label = QtGui.QLabel(GroupEditor)
        self.label.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout_4.addWidget(self.label)
        spacerItem9 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        self.verticalLayout_4.addItem(spacerItem9)
        self.gridLayout_2.addLayout(self.verticalLayout_4, 1, 0, 1, 1)
        self.line_4 = QtGui.QFrame(GroupEditor)
        self.line_4.setFrameShape(QtGui.QFrame.VLine)
        self.line_4.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout_2.addWidget(self.line_4, 4, 2, 2, 1)
        self.line_7 = QtGui.QFrame(GroupEditor)
        self.line_7.setFrameShape(QtGui.QFrame.VLine)
        self.line_7.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_7.setObjectName(_fromUtf8("line_7"))
        self.gridLayout_2.addWidget(self.line_7, 7, 2, 1, 1)
        self.line_8 = QtGui.QFrame(GroupEditor)
        self.line_8.setFrameShape(QtGui.QFrame.VLine)
        self.line_8.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_8.setObjectName(_fromUtf8("line_8"))
        self.gridLayout_2.addWidget(self.line_8, 9, 2, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_2, 0, 0, 1, 3)

        self.retranslateUi(GroupEditor)
        QtCore.QMetaObject.connectSlotsByName(GroupEditor)

    def retranslateUi(self, GroupEditor):
        GroupEditor.setWindowTitle(_translate("GroupEditor", "Dialog", None))
        self.label_3.setText(_translate("GroupEditor", "Options", None))
        self.dialogCanceled.setText(_translate("GroupEditor", "Cancel", None))
        self.dialogFinished.setText(_translate("GroupEditor", "Ok", None))
        self.label_7.setText(_translate("GroupEditor", "Color", None))
        self.label_4.setText(_translate("GroupEditor", "Minimum found ", None))
        self.omitFeatures.setText(_translate("GroupEditor", "Omit features", None))
        self.useForMetaboliteGrouping.setText(_translate("GroupEditor", "Use samples for metabolite grouping", None))
        self.label_2.setText(_translate("GroupEditor", "Files", None))
        self.label_6.setText(_translate("GroupEditor", "<html><head/><body><p>Specify the measurement files for this group</p></body></html>", None))
        self.addFolder.setText(_translate("GroupEditor", "Add folder", None))
        self.addFiles.setText(_translate("GroupEditor", "Add file(s)", None))
        self.removeSelected.setText(_translate("GroupEditor", "Remove selected", None))
        self.label_5.setText(_translate("GroupEditor", "Actions", None))
        self.label.setText(_translate("GroupEditor", "Group name", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    GroupEditor = QtGui.QDialog()
    ui = Ui_GroupEditor()
    ui.setupUi(GroupEditor)
    GroupEditor.show()
    sys.exit(app.exec_())

