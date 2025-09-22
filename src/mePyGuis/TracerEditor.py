# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\TracerEditor.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PySide6 import QtCore, QtGui, QtWidgets

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:

    def _fromUtf8(s):
        return s


try:
    _encoding = QtWidgets.QApplication.UnicodeUTF8

    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:

    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig)


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName(_fromUtf8("Dialog"))
        Dialog.resize(759, 506)
        Dialog.setStyleSheet(_fromUtf8(""))
        self.gridLayout = QtWidgets.QGridLayout(Dialog)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.gridLayout_2 = QtWidgets.QGridLayout()
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.line_4 = QtWidgets.QFrame(Dialog)
        self.line_4.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout_2.addWidget(self.line_4, 3, 2, 1, 1)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.label_2 = QtWidgets.QLabel(Dialog)
        self.label_2.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_2.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.verticalLayout_2.addWidget(self.label_2)
        spacerItem = QtWidgets.QSpacerItem(0, 0, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        self.verticalLayout_2.addItem(spacerItem)
        self.gridLayout_2.addLayout(self.verticalLayout_2, 3, 0, 1, 1)
        self.line_3 = QtWidgets.QFrame(Dialog)
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.gridLayout_2.addWidget(self.line_3, 2, 0, 1, 5)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.title_tracer = QtWidgets.QLabel(Dialog)
        self.title_tracer.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.title_tracer.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter)
        self.title_tracer.setObjectName(_fromUtf8("title_tracer"))
        self.verticalLayout.addWidget(self.title_tracer)
        spacerItem1 = QtWidgets.QSpacerItem(0, 15, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        self.verticalLayout.addItem(spacerItem1)
        self.label_3 = QtWidgets.QLabel(Dialog)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setMaximumSize(QtCore.QSize(200, 16777215))
        self.label_3.setStyleSheet(_fromUtf8("color: rgb(90, 90, 90);\nfont: 7pt;"))
        self.label_3.setAlignment(QtCore.Qt.AlignJustify | QtCore.Qt.AlignVCenter)
        self.label_3.setWordWrap(True)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.verticalLayout.addWidget(self.label_3)
        spacerItem2 = QtWidgets.QSpacerItem(200, 0, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem2)
        self.gridLayout_2.addLayout(self.verticalLayout, 1, 0, 1, 1)
        self.line = QtWidgets.QFrame(Dialog)
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout_2.addWidget(self.line, 1, 2, 1, 1)
        spacerItem3 = QtWidgets.QSpacerItem(3, 0, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem3, 1, 3, 1, 1)
        self.line_2 = QtWidgets.QFrame(Dialog)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.gridLayout_2.addWidget(self.line_2, 0, 0, 1, 5)
        spacerItem4 = QtWidgets.QSpacerItem(3, 0, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem4, 1, 1, 1, 1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem5)
        self.discardButton = QtWidgets.QPushButton(Dialog)
        self.discardButton.setObjectName(_fromUtf8("discardButton"))
        self.horizontalLayout.addWidget(self.discardButton)
        self.acceptButton = QtWidgets.QPushButton(Dialog)
        self.acceptButton.setObjectName(_fromUtf8("acceptButton"))
        self.horizontalLayout.addWidget(self.acceptButton)
        self.gridLayout_2.addLayout(self.horizontalLayout, 3, 4, 1, 1)
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label)
        self.tName = QtWidgets.QLineEdit(Dialog)
        self.tName.setObjectName(_fromUtf8("tName"))
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.tName)
        self.label_4 = QtWidgets.QLabel(Dialog)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_4)
        self.tAtomCount = QtWidgets.QSpinBox(Dialog)
        self.tAtomCount.setMinimum(1)
        self.tAtomCount.setProperty("value", 15)
        self.tAtomCount.setObjectName(_fromUtf8("tAtomCount"))
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.tAtomCount)
        self.label_5 = QtWidgets.QLabel(Dialog)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_5)
        self.label_6 = QtWidgets.QLabel(Dialog)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.label_6)
        self.label_7 = QtWidgets.QLabel(Dialog)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.label_7)
        self.label_8 = QtWidgets.QLabel(Dialog)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.label_8)
        self.tNativeEnrichment = QtWidgets.QDoubleSpinBox(Dialog)
        self.tNativeEnrichment.setProperty("value", 98.93)
        self.tNativeEnrichment.setObjectName(_fromUtf8("tNativeEnrichment"))
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.tNativeEnrichment)
        self.tLabeledEnrichment = QtWidgets.QDoubleSpinBox(Dialog)
        self.tLabeledEnrichment.setProperty("value", 99.5)
        self.tLabeledEnrichment.setObjectName(_fromUtf8("tLabeledEnrichment"))
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.tLabeledEnrichment)
        self.tNativeIsotope = QtWidgets.QLineEdit(Dialog)
        self.tNativeIsotope.setObjectName(_fromUtf8("tNativeIsotope"))
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.tNativeIsotope)
        self.tLabeledIsotope = QtWidgets.QLineEdit(Dialog)
        self.tLabeledIsotope.setObjectName(_fromUtf8("tLabeledIsotope"))
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.tLabeledIsotope)
        self.tRatioGroupBox = QtWidgets.QGroupBox(Dialog)
        self.tRatioGroupBox.setObjectName(_fromUtf8("tRatioGroupBox"))
        self.gridLayout_3 = QtWidgets.QGridLayout(self.tRatioGroupBox)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.label_9 = QtWidgets.QLabel(self.tRatioGroupBox)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.gridLayout_3.addWidget(self.label_9, 0, 0, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.tRatioGroupBox)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.gridLayout_3.addWidget(self.label_10, 1, 0, 1, 1)
        self.tAmountA = QtWidgets.QDoubleSpinBox(self.tRatioGroupBox)
        self.tAmountA.setMinimum(1.0)
        self.tAmountA.setProperty("value", 50.0)
        self.tAmountA.setObjectName(_fromUtf8("tAmountA"))
        self.gridLayout_3.addWidget(self.tAmountA, 0, 1, 1, 1)
        self.tAmountB = QtWidgets.QDoubleSpinBox(self.tRatioGroupBox)
        self.tAmountB.setProperty("value", 50.0)
        self.tAmountB.setObjectName(_fromUtf8("tAmountB"))
        self.gridLayout_3.addWidget(self.tAmountB, 1, 1, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.tRatioGroupBox)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.gridLayout_3.addWidget(self.label_11, 3, 0, 1, 1)
        self.label_12 = QtWidgets.QLabel(self.tRatioGroupBox)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.gridLayout_3.addWidget(self.label_12, 4, 0, 1, 1)
        self.tMaxRatio = QtWidgets.QDoubleSpinBox(self.tRatioGroupBox)
        self.tMaxRatio.setMaximum(100000000.0)
        self.tMaxRatio.setProperty("value", 160.0)
        self.tMaxRatio.setObjectName(_fromUtf8("tMaxRatio"))
        self.gridLayout_3.addWidget(self.tMaxRatio, 4, 1, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.tRatioGroupBox)
        self.label_14.setWordWrap(True)
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.gridLayout_3.addWidget(self.label_14, 5, 0, 1, 2)
        self.rRatioInfo = QtWidgets.QLabel(self.tRatioGroupBox)
        self.rRatioInfo.setObjectName(_fromUtf8("rRatioInfo"))
        self.gridLayout_3.addWidget(self.rRatioInfo, 6, 0, 1, 2)
        self.tMinRatio = QtWidgets.QDoubleSpinBox(self.tRatioGroupBox)
        self.tMinRatio.setProperty("value", 62.5)
        self.tMinRatio.setObjectName(_fromUtf8("tMinRatio"))
        self.gridLayout_3.addWidget(self.tMinRatio, 3, 1, 1, 1)
        self.label_13 = QtWidgets.QLabel(self.tRatioGroupBox)
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.gridLayout_3.addWidget(self.label_13, 2, 0, 1, 1)
        self.tTheoreticalSignalRatio = QtWidgets.QLabel(self.tRatioGroupBox)
        self.tTheoreticalSignalRatio.setObjectName(_fromUtf8("tTheoreticalSignalRatio"))
        self.gridLayout_3.addWidget(self.tTheoreticalSignalRatio, 2, 1, 1, 1)
        self.formLayout.setWidget(12, QtWidgets.QFormLayout.SpanningRole, self.tRatioGroupBox)
        self.tCheckRatio = QtWidgets.QCheckBox(Dialog)
        self.tCheckRatio.setObjectName(_fromUtf8("tCheckRatio"))
        self.formLayout.setWidget(10, QtWidgets.QFormLayout.FieldRole, self.tCheckRatio)
        self.label_15 = QtWidgets.QLabel(Dialog)
        self.label_15.setWordWrap(True)
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.formLayout.setWidget(11, QtWidgets.QFormLayout.SpanningRole, self.label_15)
        self.gridLayout_2.addLayout(self.formLayout, 1, 4, 1, 1)
        self.gridLayout.addLayout(self.gridLayout_2, 0, 0, 1, 1)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)
        Dialog.setTabOrder(self.tName, self.tAtomCount)
        Dialog.setTabOrder(self.tAtomCount, self.tNativeIsotope)
        Dialog.setTabOrder(self.tNativeIsotope, self.tNativeEnrichment)
        Dialog.setTabOrder(self.tNativeEnrichment, self.tLabeledIsotope)
        Dialog.setTabOrder(self.tLabeledIsotope, self.tLabeledEnrichment)
        Dialog.setTabOrder(self.tLabeledEnrichment, self.tMinRatio)
        Dialog.setTabOrder(self.tMinRatio, self.tMaxRatio)
        Dialog.setTabOrder(self.tMaxRatio, self.acceptButton)
        Dialog.setTabOrder(self.acceptButton, self.discardButton)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(_translate("Dialog", "Dialog", None))
        self.label_2.setText(_translate("Dialog", "Actions", None))
        self.title_tracer.setText(_translate("Dialog", "Tracer", None))
        self.label_3.setText(
            _translate(
                "Dialog",
                "<html><head/><body><p>Please specify the tracer added to the experiment</p><p><br/></p></body></html>",
                None,
            )
        )
        self.discardButton.setText(_translate("Dialog", "Discard", None))
        self.acceptButton.setText(_translate("Dialog", "Accept", None))
        self.label.setText(_translate("Dialog", "Name", None))
        self.tName.setText(_translate("Dialog", "DON", None))
        self.label_4.setText(_translate("Dialog", "Atom count", None))
        self.label_5.setText(_translate("Dialog", "Native isotopolog", None))
        self.label_6.setText(_translate("Dialog", "Native isotopolog enrichment", None))
        self.label_7.setText(_translate("Dialog", "Labeled isotopolog", None))
        self.label_8.setText(_translate("Dialog", "Labeled isotopolog enrichment", None))
        self.tNativeEnrichment.setSuffix(_translate("Dialog", " %", None))
        self.tLabeledEnrichment.setSuffix(_translate("Dialog", " %", None))
        self.tNativeIsotope.setText(_translate("Dialog", "12C", None))
        self.tLabeledIsotope.setText(_translate("Dialog", "13C", None))
        self.tRatioGroupBox.setTitle(_translate("Dialog", "Virtually allowed ratio for exogenous tracers only", None))
        self.label_9.setText(_translate("Dialog", "Amount native compound", None))
        self.label_10.setText(_translate("Dialog", "Amount labeled compound", None))
        self.tAmountA.setPrefix(_translate("Dialog", "= ", None))
        self.tAmountB.setPrefix(_translate("Dialog", "= ", None))
        self.label_11.setText(
            _translate(
                "Dialog",
                "Minimum ratio deviation relative to theoretical signal ratio",
                None,
            )
        )
        self.label_12.setText(
            _translate(
                "Dialog",
                "Maximum ratio deviation relative to theoretical signal ratio",
                None,
            )
        )
        self.tMaxRatio.setPrefix(_translate("Dialog", "≤ ", None))
        self.tMaxRatio.setSuffix(_translate("Dialog", " %", None))
        self.label_14.setText(
            _translate(
                "Dialog",
                "Allowed ratios in the LC-HRMS data calculated from the isotopic enrichments and the spiked amounts:",
                None,
            )
        )
        self.rRatioInfo.setText(_translate("Dialog", "TextLabel", None))
        self.tMinRatio.setPrefix(_translate("Dialog", "≥ ", None))
        self.tMinRatio.setSuffix(_translate("Dialog", " %", None))
        self.label_13.setText(_translate("Dialog", "Theoretical signal ratio", None))
        self.tTheoreticalSignalRatio.setText(_translate("Dialog", "-", None))
        self.tCheckRatio.setText(
            _translate(
                "Dialog",
                "Exogenous tracer (i.e. check ratio of native and labeled tracer form)",
                None,
            )
        )
        self.label_15.setText(
            _translate(
                "Dialog",
                "Note: Activate the ratio check for exogenous tracer compounds and deactivate it for endogenous tracer compounds",
                None,
            )
        )


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec())
