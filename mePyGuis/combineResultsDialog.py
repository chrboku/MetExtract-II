# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\combineResultsDialog.ui'
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
        Dialog.resize(602, 750)
        self.gridLayout = QtGui.QGridLayout(Dialog)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_16 = QtGui.QLabel(Dialog)
        self.label_16.setWordWrap(True)
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.gridLayout.addWidget(self.label_16, 0, 0, 1, 1)
        self.groupBox = QtGui.QGroupBox(Dialog)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout_2 = QtGui.QGridLayout(self.groupBox)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.line_4 = QtGui.QFrame(self.groupBox)
        self.line_4.setFrameShape(QtGui.QFrame.HLine)
        self.line_4.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout_2.addWidget(self.line_4, 14, 1, 1, 3)
        spacerItem = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem, 15, 1, 1, 1)
        self.expDLoad = QtGui.QPushButton(self.groupBox)
        self.expDLoad.setObjectName(_fromUtf8("expDLoad"))
        self.gridLayout_2.addWidget(self.expDLoad, 12, 3, 1, 1)
        self.expELoad = QtGui.QPushButton(self.groupBox)
        self.expELoad.setObjectName(_fromUtf8("expELoad"))
        self.gridLayout_2.addWidget(self.expELoad, 16, 3, 1, 1)
        self.expFLoad = QtGui.QPushButton(self.groupBox)
        self.expFLoad.setObjectName(_fromUtf8("expFLoad"))
        self.gridLayout_2.addWidget(self.expFLoad, 20, 3, 1, 1)
        self.label = QtGui.QLabel(self.groupBox)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout_2.addWidget(self.label, 0, 1, 1, 1)
        self.expCPrefix = QtGui.QLineEdit(self.groupBox)
        self.expCPrefix.setObjectName(_fromUtf8("expCPrefix"))
        self.gridLayout_2.addWidget(self.expCPrefix, 9, 2, 1, 2)
        self.expDPrefix = QtGui.QLineEdit(self.groupBox)
        self.expDPrefix.setObjectName(_fromUtf8("expDPrefix"))
        self.gridLayout_2.addWidget(self.expDPrefix, 13, 2, 1, 2)
        self.expESave = QtGui.QLineEdit(self.groupBox)
        self.expESave.setObjectName(_fromUtf8("expESave"))
        self.gridLayout_2.addWidget(self.expESave, 16, 2, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem1, 7, 1, 1, 1)
        self.expBPrefix = QtGui.QLineEdit(self.groupBox)
        self.expBPrefix.setObjectName(_fromUtf8("expBPrefix"))
        self.gridLayout_2.addWidget(self.expBPrefix, 5, 2, 1, 2)
        self.expFSave = QtGui.QLineEdit(self.groupBox)
        self.expFSave.setObjectName(_fromUtf8("expFSave"))
        self.gridLayout_2.addWidget(self.expFSave, 20, 2, 1, 1)
        self.expFPrefix = QtGui.QLineEdit(self.groupBox)
        self.expFPrefix.setObjectName(_fromUtf8("expFPrefix"))
        self.gridLayout_2.addWidget(self.expFPrefix, 21, 2, 1, 2)
        self.expEPrefix = QtGui.QLineEdit(self.groupBox)
        self.expEPrefix.setObjectName(_fromUtf8("expEPrefix"))
        self.gridLayout_2.addWidget(self.expEPrefix, 17, 2, 1, 2)
        spacerItem2 = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem2, 19, 1, 1, 1)
        self.expASave = QtGui.QLineEdit(self.groupBox)
        self.expASave.setObjectName(_fromUtf8("expASave"))
        self.gridLayout_2.addWidget(self.expASave, 0, 2, 1, 1)
        spacerItem3 = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem3, 11, 1, 1, 1)
        self.label_5 = QtGui.QLabel(self.groupBox)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_2.addWidget(self.label_5, 8, 1, 1, 1)
        self.line = QtGui.QFrame(self.groupBox)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout_2.addWidget(self.line, 2, 1, 1, 3)
        spacerItem4 = QtGui.QSpacerItem(20, 10, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_2.addItem(spacerItem4, 3, 1, 1, 1)
        self.label_3 = QtGui.QLabel(self.groupBox)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_2.addWidget(self.label_3, 4, 1, 1, 1)
        self.expAPrefix = QtGui.QLineEdit(self.groupBox)
        self.expAPrefix.setObjectName(_fromUtf8("expAPrefix"))
        self.gridLayout_2.addWidget(self.expAPrefix, 1, 2, 1, 2)
        self.label_4 = QtGui.QLabel(self.groupBox)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout_2.addWidget(self.label_4, 5, 1, 1, 1)
        self.label_2 = QtGui.QLabel(self.groupBox)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout_2.addWidget(self.label_2, 1, 1, 1, 1)
        self.expALoad = QtGui.QPushButton(self.groupBox)
        self.expALoad.setObjectName(_fromUtf8("expALoad"))
        self.gridLayout_2.addWidget(self.expALoad, 0, 3, 1, 1)
        self.label_6 = QtGui.QLabel(self.groupBox)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout_2.addWidget(self.label_6, 9, 1, 1, 1)
        self.label_7 = QtGui.QLabel(self.groupBox)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout_2.addWidget(self.label_7, 12, 1, 1, 1)
        self.label_8 = QtGui.QLabel(self.groupBox)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridLayout_2.addWidget(self.label_8, 13, 1, 1, 1)
        self.label_9 = QtGui.QLabel(self.groupBox)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.gridLayout_2.addWidget(self.label_9, 16, 1, 1, 1)
        self.label_10 = QtGui.QLabel(self.groupBox)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.gridLayout_2.addWidget(self.label_10, 17, 1, 1, 1)
        self.line_2 = QtGui.QFrame(self.groupBox)
        self.line_2.setFrameShape(QtGui.QFrame.HLine)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.gridLayout_2.addWidget(self.line_2, 6, 1, 1, 3)
        self.expCLoad = QtGui.QPushButton(self.groupBox)
        self.expCLoad.setObjectName(_fromUtf8("expCLoad"))
        self.gridLayout_2.addWidget(self.expCLoad, 8, 3, 1, 1)
        self.label_11 = QtGui.QLabel(self.groupBox)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.gridLayout_2.addWidget(self.label_11, 20, 1, 1, 1)
        self.expBSave = QtGui.QLineEdit(self.groupBox)
        self.expBSave.setObjectName(_fromUtf8("expBSave"))
        self.gridLayout_2.addWidget(self.expBSave, 4, 2, 1, 1)
        self.label_12 = QtGui.QLabel(self.groupBox)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.gridLayout_2.addWidget(self.label_12, 21, 1, 1, 1)
        self.expBLoad = QtGui.QPushButton(self.groupBox)
        self.expBLoad.setObjectName(_fromUtf8("expBLoad"))
        self.gridLayout_2.addWidget(self.expBLoad, 4, 3, 1, 1)
        self.line_3 = QtGui.QFrame(self.groupBox)
        self.line_3.setFrameShape(QtGui.QFrame.HLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.gridLayout_2.addWidget(self.line_3, 10, 1, 1, 3)
        self.line_5 = QtGui.QFrame(self.groupBox)
        self.line_5.setFrameShape(QtGui.QFrame.HLine)
        self.line_5.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_5.setObjectName(_fromUtf8("line_5"))
        self.gridLayout_2.addWidget(self.line_5, 18, 1, 1, 3)
        self.expCSave = QtGui.QLineEdit(self.groupBox)
        self.expCSave.setObjectName(_fromUtf8("expCSave"))
        self.gridLayout_2.addWidget(self.expCSave, 8, 2, 1, 1)
        self.expDSave = QtGui.QLineEdit(self.groupBox)
        self.expDSave.setObjectName(_fromUtf8("expDSave"))
        self.gridLayout_2.addWidget(self.expDSave, 12, 2, 1, 1)
        self.label_17 = QtGui.QLabel(self.groupBox)
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.gridLayout_2.addWidget(self.label_17, 0, 0, 2, 1)
        self.label_18 = QtGui.QLabel(self.groupBox)
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.gridLayout_2.addWidget(self.label_18, 4, 0, 2, 1)
        self.label_19 = QtGui.QLabel(self.groupBox)
        self.label_19.setObjectName(_fromUtf8("label_19"))
        self.gridLayout_2.addWidget(self.label_19, 8, 0, 2, 1)
        self.label_20 = QtGui.QLabel(self.groupBox)
        self.label_20.setObjectName(_fromUtf8("label_20"))
        self.gridLayout_2.addWidget(self.label_20, 12, 0, 2, 1)
        self.label_21 = QtGui.QLabel(self.groupBox)
        self.label_21.setObjectName(_fromUtf8("label_21"))
        self.gridLayout_2.addWidget(self.label_21, 16, 0, 2, 1)
        self.label_22 = QtGui.QLabel(self.groupBox)
        self.label_22.setObjectName(_fromUtf8("label_22"))
        self.gridLayout_2.addWidget(self.label_22, 20, 0, 2, 1)
        self.gridLayout.addWidget(self.groupBox, 2, 0, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setSizeConstraint(QtGui.QLayout.SetMinimumSize)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.label_15 = QtGui.QLabel(Dialog)
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.horizontalLayout.addWidget(self.label_15)
        self.saveResults = QtGui.QLineEdit(Dialog)
        self.saveResults.setObjectName(_fromUtf8("saveResults"))
        self.horizontalLayout.addWidget(self.saveResults)
        self.loadResults = QtGui.QPushButton(Dialog)
        self.loadResults.setObjectName(_fromUtf8("loadResults"))
        self.horizontalLayout.addWidget(self.loadResults)
        self.gridLayout.addLayout(self.horizontalLayout, 4, 0, 1, 1)
        self.groupBox_2 = QtGui.QGroupBox(Dialog)
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.gridLayout_3 = QtGui.QGridLayout(self.groupBox_2)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.label_13 = QtGui.QLabel(self.groupBox_2)
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.gridLayout_3.addWidget(self.label_13, 0, 0, 1, 1)
        self.label_14 = QtGui.QLabel(self.groupBox_2)
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.gridLayout_3.addWidget(self.label_14, 1, 0, 1, 1)
        self.maxPPMDev = QtGui.QDoubleSpinBox(self.groupBox_2)
        self.maxPPMDev.setMinimum(0.0)
        self.maxPPMDev.setProperty("value", 5.0)
        self.maxPPMDev.setObjectName(_fromUtf8("maxPPMDev"))
        self.gridLayout_3.addWidget(self.maxPPMDev, 0, 1, 1, 1)
        self.maxRTDev = QtGui.QDoubleSpinBox(self.groupBox_2)
        self.maxRTDev.setProperty("value", 0.15)
        self.maxRTDev.setObjectName(_fromUtf8("maxRTDev"))
        self.gridLayout_3.addWidget(self.maxRTDev, 1, 1, 1, 1)
        self.gridLayout.addWidget(self.groupBox_2, 3, 0, 1, 1)
        spacerItem5 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem5, 0, 1, 4, 1)
        self.line_6 = QtGui.QFrame(Dialog)
        self.line_6.setFrameShape(QtGui.QFrame.HLine)
        self.line_6.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_6.setObjectName(_fromUtf8("line_6"))
        self.gridLayout.addWidget(self.line_6, 5, 0, 1, 1)
        spacerItem6 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        self.gridLayout.addItem(spacerItem6, 1, 0, 1, 1)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        spacerItem7 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem7)
        self.run = QtGui.QPushButton(Dialog)
        self.run.setObjectName(_fromUtf8("run"))
        self.horizontalLayout_3.addWidget(self.run)
        self.gridLayout.addLayout(self.horizontalLayout_3, 6, 0, 1, 1)
        self.label_18.setBuddy(self.expBSave)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)
        Dialog.setTabOrder(self.expASave, self.expALoad)
        Dialog.setTabOrder(self.expALoad, self.expAPrefix)
        Dialog.setTabOrder(self.expAPrefix, self.expBSave)
        Dialog.setTabOrder(self.expBSave, self.expBLoad)
        Dialog.setTabOrder(self.expBLoad, self.expBPrefix)
        Dialog.setTabOrder(self.expBPrefix, self.expCSave)
        Dialog.setTabOrder(self.expCSave, self.expCLoad)
        Dialog.setTabOrder(self.expCLoad, self.expCPrefix)
        Dialog.setTabOrder(self.expCPrefix, self.expDSave)
        Dialog.setTabOrder(self.expDSave, self.expDLoad)
        Dialog.setTabOrder(self.expDLoad, self.expDPrefix)
        Dialog.setTabOrder(self.expDPrefix, self.expESave)
        Dialog.setTabOrder(self.expESave, self.expELoad)
        Dialog.setTabOrder(self.expELoad, self.expEPrefix)
        Dialog.setTabOrder(self.expEPrefix, self.expFSave)
        Dialog.setTabOrder(self.expFSave, self.expFLoad)
        Dialog.setTabOrder(self.expFLoad, self.expFPrefix)
        Dialog.setTabOrder(self.expFPrefix, self.maxPPMDev)
        Dialog.setTabOrder(self.maxPPMDev, self.maxRTDev)
        Dialog.setTabOrder(self.maxRTDev, self.saveResults)
        Dialog.setTabOrder(self.saveResults, self.loadResults)
        Dialog.setTabOrder(self.loadResults, self.run)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "MetExtract II", None))
        self.label_16.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:12pt; font-weight:600;\">Combine results</span></p><p>Combine the data processing results from several MetExtract II experiments. For example, if three experiments (U-13C-labeling, U-15N-labeling and a tracer-biotransformation experiment) have been carried out, the results can be combined into a singe data matrix using this tool. The tool is a very basic tool and does not perform any sort of chromatographic alignment or a correction for m/z shifts across the chromatograms. </p></body></html>", None))
        self.groupBox.setTitle(_translate("Dialog", "Input experiments", None))
        self.expDLoad.setText(_translate("Dialog", "Load", None))
        self.expELoad.setText(_translate("Dialog", "Load", None))
        self.expFLoad.setText(_translate("Dialog", "Load", None))
        self.label.setText(_translate("Dialog", "Results experiment A", None))
        self.label_5.setText(_translate("Dialog", "Results experiment C", None))
        self.label_3.setText(_translate("Dialog", "Results experiment B", None))
        self.label_4.setText(_translate("Dialog", "Prefix experiment B", None))
        self.label_2.setText(_translate("Dialog", "Prefix experiment A", None))
        self.expALoad.setText(_translate("Dialog", "Load", None))
        self.label_6.setText(_translate("Dialog", "Prefix experiment C", None))
        self.label_7.setText(_translate("Dialog", "Results experiment D", None))
        self.label_8.setText(_translate("Dialog", "Prefix experiment D", None))
        self.label_9.setText(_translate("Dialog", "Results experiment E", None))
        self.label_10.setText(_translate("Dialog", "Prefix experiment E", None))
        self.expCLoad.setText(_translate("Dialog", "Load", None))
        self.label_11.setText(_translate("Dialog", "Results experiment F", None))
        self.label_12.setText(_translate("Dialog", "Prefix experiment F", None))
        self.expBLoad.setText(_translate("Dialog", "Load", None))
        self.label_17.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:32pt;\">A</span></p></body></html>", None))
        self.label_18.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:32pt;\">B</span></p></body></html>", None))
        self.label_19.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:32pt; color:#8b8b8b;\">C</span></p></body></html>", None))
        self.label_20.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:32pt; color:#8b8b8b;\">D</span></p></body></html>", None))
        self.label_21.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:32pt; color:#8b8b8b;\">E</span></p></body></html>", None))
        self.label_22.setText(_translate("Dialog", "<html><head/><body><p><span style=\" font-size:32pt; color:#8b8b8b;\">F</span></p></body></html>", None))
        self.label_15.setText(_translate("Dialog", "Save combined results to", None))
        self.loadResults.setText(_translate("Dialog", "Select", None))
        self.groupBox_2.setTitle(_translate("Dialog", "Processing parameters", None))
        self.label_13.setText(_translate("Dialog", "Maximum m/z deviation (± ppm)", None))
        self.label_14.setText(_translate("Dialog", "Maximum retentiontime deviation (± minutes)", None))
        self.run.setText(_translate("Dialog", "Run", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
