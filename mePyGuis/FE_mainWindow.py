# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\FE_mainWindow.ui'
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

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1391, 840)
        MainWindow.setAcceptDrops(True)
        MainWindow.setTabShape(QtGui.QTabWidget.Triangular)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.tabWidget = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget.setAcceptDrops(True)
        self.tabWidget.setTabPosition(QtGui.QTabWidget.West)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.inputTab = QtGui.QWidget()
        self.inputTab.setObjectName(_fromUtf8("inputTab"))
        self.gridLayout_2 = QtGui.QGridLayout(self.inputTab)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.gridLayout_6 = QtGui.QGridLayout()
        self.gridLayout_6.setObjectName(_fromUtf8("gridLayout_6"))
        self.gridLayout_7 = QtGui.QGridLayout()
        self.gridLayout_7.setObjectName(_fromUtf8("gridLayout_7"))
        self.label_4 = QtGui.QLabel(self.inputTab)
        self.label_4.setMinimumSize(QtCore.QSize(200, 0))
        self.label_4.setStyleSheet(_fromUtf8("color: rgb(90, 90, 90);\n"
"font: 7pt;"))
        self.label_4.setWordWrap(True)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout_7.addWidget(self.label_4, 2, 0, 1, 1)
        self.label_3 = QtGui.QLabel(self.inputTab)
        self.label_3.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_7.addWidget(self.label_3, 0, 0, 1, 1)
        spacerItem = QtGui.QSpacerItem(0, 5, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Maximum)
        self.gridLayout_7.addItem(spacerItem, 1, 0, 1, 1)
        spacerItem1 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_7.addItem(spacerItem1, 3, 0, 1, 1)
        self.gridLayout_6.addLayout(self.gridLayout_7, 0, 0, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.saveCompilation = QtGui.QPushButton(self.inputTab)
        self.saveCompilation.setObjectName(_fromUtf8("saveCompilation"))
        self.horizontalLayout.addWidget(self.saveCompilation)
        self.loadCompilation = QtGui.QPushButton(self.inputTab)
        self.loadCompilation.setObjectName(_fromUtf8("loadCompilation"))
        self.horizontalLayout.addWidget(self.loadCompilation)
        spacerItem2 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem2)
        self.addMSMSTarget = QtGui.QPushButton(self.inputTab)
        self.addMSMSTarget.setObjectName(_fromUtf8("addMSMSTarget"))
        self.horizontalLayout.addWidget(self.addMSMSTarget)
        self.deleteMSMSTarget = QtGui.QPushButton(self.inputTab)
        self.deleteMSMSTarget.setObjectName(_fromUtf8("deleteMSMSTarget"))
        self.horizontalLayout.addWidget(self.deleteMSMSTarget)
        self.gridLayout_6.addLayout(self.horizontalLayout, 3, 4, 1, 1)
        spacerItem3 = QtGui.QSpacerItem(3, 0, QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum)
        self.gridLayout_6.addItem(spacerItem3, 0, 3, 1, 1)
        spacerItem4 = QtGui.QSpacerItem(3, 0, QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Minimum)
        self.gridLayout_6.addItem(spacerItem4, 0, 1, 1, 1)
        self.label_5 = QtGui.QLabel(self.inputTab)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_6.addWidget(self.label_5, 3, 0, 1, 1)
        self.line = QtGui.QFrame(self.inputTab)
        self.line.setFrameShape(QtGui.QFrame.VLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout_6.addWidget(self.line, 0, 2, 1, 1)
        self.line_2 = QtGui.QFrame(self.inputTab)
        self.line_2.setFrameShape(QtGui.QFrame.HLine)
        self.line_2.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.gridLayout_6.addWidget(self.line_2, 1, 0, 1, 5)
        self.line_3 = QtGui.QFrame(self.inputTab)
        self.line_3.setFrameShape(QtGui.QFrame.VLine)
        self.line_3.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.gridLayout_6.addWidget(self.line_3, 3, 2, 1, 1)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.processFilesTable = QtGui.QTableView(self.inputTab)
        self.processFilesTable.setAlternatingRowColors(False)
        self.processFilesTable.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
        self.processFilesTable.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.processFilesTable.setObjectName(_fromUtf8("processFilesTable"))
        self.verticalLayout_2.addWidget(self.processFilesTable)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.configAdductsButton = QtGui.QPushButton(self.inputTab)
        self.configAdductsButton.setObjectName(_fromUtf8("configAdductsButton"))
        self.horizontalLayout_2.addWidget(self.configAdductsButton)
        spacerItem5 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem5)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.gridLayout_6.addLayout(self.verticalLayout_2, 0, 4, 1, 1)
        self.gridLayout_2.addLayout(self.gridLayout_6, 0, 0, 1, 1)
        self.tabWidget.addTab(self.inputTab, _fromUtf8(""))
        self.processTab = QtGui.QWidget()
        self.processTab.setObjectName(_fromUtf8("processTab"))
        self.gridLayout_8 = QtGui.QGridLayout(self.processTab)
        self.gridLayout_8.setObjectName(_fromUtf8("gridLayout_8"))
        self.gridLayout_9 = QtGui.QGridLayout()
        self.gridLayout_9.setObjectName(_fromUtf8("gridLayout_9"))
        self.line_5 = QtGui.QFrame(self.processTab)
        self.line_5.setFrameShape(QtGui.QFrame.HLine)
        self.line_5.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_5.setObjectName(_fromUtf8("line_5"))
        self.gridLayout_9.addWidget(self.line_5, 1, 0, 1, 6)
        spacerItem6 = QtGui.QSpacerItem(3, 0, QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum)
        self.gridLayout_9.addItem(spacerItem6, 0, 3, 1, 1)
        spacerItem7 = QtGui.QSpacerItem(3, 0, QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum)
        self.gridLayout_9.addItem(spacerItem7, 0, 1, 1, 1)
        self.gridLayout_10 = QtGui.QGridLayout()
        self.gridLayout_10.setObjectName(_fromUtf8("gridLayout_10"))
        spacerItem8 = QtGui.QSpacerItem(0, 5, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Maximum)
        self.gridLayout_10.addItem(spacerItem8, 1, 0, 1, 1)
        self.label_6 = QtGui.QLabel(self.processTab)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_6.sizePolicy().hasHeightForWidth())
        self.label_6.setSizePolicy(sizePolicy)
        self.label_6.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout_10.addWidget(self.label_6, 0, 0, 1, 1)
        self.label = QtGui.QLabel(self.processTab)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setMinimumSize(QtCore.QSize(200, 0))
        self.label.setStyleSheet(_fromUtf8("color: rgb(90, 90, 90);\n"
"font: 7pt;"))
        self.label.setWordWrap(True)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout_10.addWidget(self.label, 2, 0, 1, 1)
        spacerItem9 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        self.gridLayout_10.addItem(spacerItem9, 3, 0, 1, 1)
        self.gridLayout_9.addLayout(self.gridLayout_10, 0, 0, 1, 1)
        self.line_4 = QtGui.QFrame(self.processTab)
        self.line_4.setFrameShape(QtGui.QFrame.VLine)
        self.line_4.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout_9.addWidget(self.line_4, 0, 2, 1, 1)
        self.gridLayout_11 = QtGui.QGridLayout()
        self.gridLayout_11.setObjectName(_fromUtf8("gridLayout_11"))
        self.label_7 = QtGui.QLabel(self.processTab)
        self.label_7.setStyleSheet(_fromUtf8("font: 10pt;"))
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout_11.addWidget(self.label_7, 0, 0, 1, 1)
        self.gridLayout_9.addLayout(self.gridLayout_11, 2, 0, 1, 1)
        self.gridLayout_12 = QtGui.QGridLayout()
        self.gridLayout_12.setObjectName(_fromUtf8("gridLayout_12"))
        self.label_8 = QtGui.QLabel(self.processTab)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridLayout_12.addWidget(self.label_8, 0, 2, 1, 1)
        self.useCPUCoresSpinner = QtGui.QSpinBox(self.processTab)
        self.useCPUCoresSpinner.setMinimum(1)
        self.useCPUCoresSpinner.setObjectName(_fromUtf8("useCPUCoresSpinner"))
        self.gridLayout_12.addWidget(self.useCPUCoresSpinner, 0, 3, 1, 1)
        spacerItem10 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_12.addItem(spacerItem10, 0, 0, 1, 1)
        self.keepCPUCoreUnusedCheckbox = QtGui.QCheckBox(self.processTab)
        self.keepCPUCoreUnusedCheckbox.setObjectName(_fromUtf8("keepCPUCoreUnusedCheckbox"))
        self.gridLayout_12.addWidget(self.keepCPUCoreUnusedCheckbox, 0, 1, 1, 1)
        self.startProcessingButton = QtGui.QPushButton(self.processTab)
        self.startProcessingButton.setObjectName(_fromUtf8("startProcessingButton"))
        self.gridLayout_12.addWidget(self.startProcessingButton, 0, 4, 1, 1)
        self.gridLayout_9.addLayout(self.gridLayout_12, 2, 5, 1, 1)
        self.gridLayout_14 = QtGui.QGridLayout()
        self.gridLayout_14.setObjectName(_fromUtf8("gridLayout_14"))
        self.groupBox = QtGui.QGroupBox(self.processTab)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout_15 = QtGui.QGridLayout(self.groupBox)
        self.gridLayout_15.setObjectName(_fromUtf8("gridLayout_15"))
        self.label_9 = QtGui.QLabel(self.groupBox)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.gridLayout_15.addWidget(self.label_9, 0, 0, 1, 1)
        self.eicPPMSpinner = QtGui.QDoubleSpinBox(self.groupBox)
        self.eicPPMSpinner.setProperty("value", 5.0)
        self.eicPPMSpinner.setObjectName(_fromUtf8("eicPPMSpinner"))
        self.gridLayout_15.addWidget(self.eicPPMSpinner, 1, 0, 1, 1)
        self.label_10 = QtGui.QLabel(self.groupBox)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.gridLayout_15.addWidget(self.label_10, 2, 0, 1, 1)
        self.intensityThresoldSpinner = QtGui.QDoubleSpinBox(self.groupBox)
        self.intensityThresoldSpinner.setDecimals(0)
        self.intensityThresoldSpinner.setMaximum(999999999.0)
        self.intensityThresoldSpinner.setProperty("value", 10000.0)
        self.intensityThresoldSpinner.setObjectName(_fromUtf8("intensityThresoldSpinner"))
        self.gridLayout_15.addWidget(self.intensityThresoldSpinner, 3, 0, 1, 1)
        spacerItem11 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        self.gridLayout_15.addItem(spacerItem11, 4, 0, 1, 1)
        self.gridLayout_14.addWidget(self.groupBox, 0, 0, 1, 1)
        self.groupBox_3 = QtGui.QGroupBox(self.processTab)
        self.groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
        self.gridLayout_16 = QtGui.QGridLayout(self.groupBox_3)
        self.gridLayout_16.setObjectName(_fromUtf8("gridLayout_16"))
        self.label_12 = QtGui.QLabel(self.groupBox_3)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.gridLayout_16.addWidget(self.label_12, 0, 0, 1, 1)
        self.label_13 = QtGui.QLabel(self.groupBox_3)
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.gridLayout_16.addWidget(self.label_13, 2, 0, 1, 1)
        self.maxPPMErrorMatching_Spinner = QtGui.QDoubleSpinBox(self.groupBox_3)
        self.maxPPMErrorMatching_Spinner.setMaximum(9999.0)
        self.maxPPMErrorMatching_Spinner.setProperty("value", 30.0)
        self.maxPPMErrorMatching_Spinner.setObjectName(_fromUtf8("maxPPMErrorMatching_Spinner"))
        self.gridLayout_16.addWidget(self.maxPPMErrorMatching_Spinner, 1, 0, 1, 1)
        self.maxRelIntensityError_spinner = QtGui.QDoubleSpinBox(self.groupBox_3)
        self.maxRelIntensityError_spinner.setMaximum(100.0)
        self.maxRelIntensityError_spinner.setProperty("value", 20.0)
        self.maxRelIntensityError_spinner.setObjectName(_fromUtf8("maxRelIntensityError_spinner"))
        self.gridLayout_16.addWidget(self.maxRelIntensityError_spinner, 3, 0, 1, 1)
        spacerItem12 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        self.gridLayout_16.addItem(spacerItem12, 4, 0, 1, 1)
        self.gridLayout_14.addWidget(self.groupBox_3, 0, 3, 1, 1)
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.groupBox_2 = QtGui.QGroupBox(self.processTab)
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.gridLayout_13 = QtGui.QGridLayout(self.groupBox_2)
        self.gridLayout_13.setObjectName(_fromUtf8("gridLayout_13"))
        self.label_11 = QtGui.QLabel(self.groupBox_2)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.gridLayout_13.addWidget(self.label_11, 0, 0, 1, 1)
        self.minScaledPeakIntensity_Spinner = QtGui.QDoubleSpinBox(self.groupBox_2)
        self.minScaledPeakIntensity_Spinner.setDecimals(1)
        self.minScaledPeakIntensity_Spinner.setMinimum(0.0)
        self.minScaledPeakIntensity_Spinner.setProperty("value", 0.5)
        self.minScaledPeakIntensity_Spinner.setObjectName(_fromUtf8("minScaledPeakIntensity_Spinner"))
        self.gridLayout_13.addWidget(self.minScaledPeakIntensity_Spinner, 1, 0, 1, 1)
        self.scalePrecursorMZ = QtGui.QCheckBox(self.groupBox_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.scalePrecursorMZ.sizePolicy().hasHeightForWidth())
        self.scalePrecursorMZ.setSizePolicy(sizePolicy)
        self.scalePrecursorMZ.setObjectName(_fromUtf8("scalePrecursorMZ"))
        self.gridLayout_13.addWidget(self.scalePrecursorMZ, 2, 0, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_2)
        spacerItem13 = QtGui.QSpacerItem(0, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Ignored)
        self.verticalLayout.addItem(spacerItem13)
        self.gridLayout_14.addLayout(self.verticalLayout, 0, 1, 1, 1)
        self.groupBox_4 = QtGui.QGroupBox(self.processTab)
        self.groupBox_4.setObjectName(_fromUtf8("groupBox_4"))
        self.gridLayout_17 = QtGui.QGridLayout(self.groupBox_4)
        self.gridLayout_17.setObjectName(_fromUtf8("gridLayout_17"))
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        self.fragmentAnnotationButton = QtGui.QPushButton(self.groupBox_4)
        self.fragmentAnnotationButton.setObjectName(_fromUtf8("fragmentAnnotationButton"))
        self.horizontalLayout_3.addWidget(self.fragmentAnnotationButton)
        spacerItem14 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem14)
        self.gridLayout_17.addLayout(self.horizontalLayout_3, 2, 0, 1, 1)
        self.label_14 = QtGui.QLabel(self.groupBox_4)
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.gridLayout_17.addWidget(self.label_14, 3, 0, 1, 1)
        self.useTracExtractAnnotation = QtGui.QCheckBox(self.groupBox_4)
        self.useTracExtractAnnotation.setObjectName(_fromUtf8("useTracExtractAnnotation"))
        self.gridLayout_17.addWidget(self.useTracExtractAnnotation, 6, 0, 1, 1)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        self.annotationPPMErrorSpinner = QtGui.QDoubleSpinBox(self.groupBox_4)
        self.annotationPPMErrorSpinner.setDecimals(1)
        self.annotationPPMErrorSpinner.setMaximum(1000.0)
        self.annotationPPMErrorSpinner.setProperty("value", 5.0)
        self.annotationPPMErrorSpinner.setObjectName(_fromUtf8("annotationPPMErrorSpinner"))
        self.horizontalLayout_4.addWidget(self.annotationPPMErrorSpinner)
        spacerItem15 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem15)
        self.gridLayout_17.addLayout(self.horizontalLayout_4, 5, 0, 1, 1)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        self.label_15 = QtGui.QLabel(self.groupBox_4)
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.horizontalLayout_5.addWidget(self.label_15)
        self.minXn = QtGui.QSpinBox(self.groupBox_4)
        self.minXn.setObjectName(_fromUtf8("minXn"))
        self.horizontalLayout_5.addWidget(self.minXn)
        self.gridLayout_17.addLayout(self.horizontalLayout_5, 0, 0, 1, 1)
        self.useZeroLabelingAtoms = QtGui.QCheckBox(self.groupBox_4)
        self.useZeroLabelingAtoms.setObjectName(_fromUtf8("useZeroLabelingAtoms"))
        self.gridLayout_17.addWidget(self.useZeroLabelingAtoms, 1, 0, 1, 1)
        self.gridLayout_14.addWidget(self.groupBox_4, 1, 0, 1, 2)
        self.groupBox_5 = QtGui.QGroupBox(self.processTab)
        self.groupBox_5.setObjectName(_fromUtf8("groupBox_5"))
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.groupBox_5)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.applyParentFragmentConsistencyRuleCheckBox = QtGui.QCheckBox(self.groupBox_5)
        self.applyParentFragmentConsistencyRuleCheckBox.setObjectName(_fromUtf8("applyParentFragmentConsistencyRuleCheckBox"))
        self.verticalLayout_3.addWidget(self.applyParentFragmentConsistencyRuleCheckBox)
        spacerItem16 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Preferred)
        self.verticalLayout_3.addItem(spacerItem16)
        self.gridLayout_14.addWidget(self.groupBox_5, 1, 3, 1, 1)
        self.groupBox_6 = QtGui.QGroupBox(self.processTab)
        self.groupBox_6.setObjectName(_fromUtf8("groupBox_6"))
        self.gridLayout_18 = QtGui.QGridLayout(self.groupBox_6)
        self.gridLayout_18.setObjectName(_fromUtf8("gridLayout_18"))
        self.saveResultsAsTSVCheckBox = QtGui.QCheckBox(self.groupBox_6)
        self.saveResultsAsTSVCheckBox.setChecked(True)
        self.saveResultsAsTSVCheckBox.setObjectName(_fromUtf8("saveResultsAsTSVCheckBox"))
        self.gridLayout_18.addWidget(self.saveResultsAsTSVCheckBox, 0, 0, 1, 1)
        self.saveResultsAsPDFCheckBox = QtGui.QCheckBox(self.groupBox_6)
        self.saveResultsAsPDFCheckBox.setObjectName(_fromUtf8("saveResultsAsPDFCheckBox"))
        self.gridLayout_18.addWidget(self.saveResultsAsPDFCheckBox, 1, 0, 1, 1)
        self.gridLayout_14.addWidget(self.groupBox_6, 2, 0, 1, 2)
        self.gridLayout_9.addLayout(self.gridLayout_14, 0, 5, 1, 1)
        self.line_6 = QtGui.QFrame(self.processTab)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.line_6.sizePolicy().hasHeightForWidth())
        self.line_6.setSizePolicy(sizePolicy)
        self.line_6.setFrameShape(QtGui.QFrame.VLine)
        self.line_6.setFrameShadow(QtGui.QFrame.Sunken)
        self.line_6.setObjectName(_fromUtf8("line_6"))
        self.gridLayout_9.addWidget(self.line_6, 2, 2, 1, 1)
        self.gridLayout_8.addLayout(self.gridLayout_9, 0, 0, 1, 1)
        spacerItem17 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_8.addItem(spacerItem17, 1, 0, 1, 1)
        spacerItem18 = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_8.addItem(spacerItem18, 0, 1, 1, 1)
        self.gridLayout_8.setColumnStretch(1, 2)
        self.gridLayout_8.setRowStretch(1, 2)
        self.tabWidget.addTab(self.processTab, _fromUtf8(""))
        self.resultsTab = QtGui.QWidget()
        self.resultsTab.setObjectName(_fromUtf8("resultsTab"))
        self.gridLayout_3 = QtGui.QGridLayout(self.resultsTab)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.gridLayout_4 = QtGui.QGridLayout()
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.gridLayout_5 = QtGui.QGridLayout()
        self.gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
        self.label_2 = QtGui.QLabel(self.resultsTab)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout_5.addWidget(self.label_2, 0, 0, 1, 1)
        self.processedFilesComboBox = QtGui.QComboBox(self.resultsTab)
        self.processedFilesComboBox.setMinimumSize(QtCore.QSize(150, 0))
        self.processedFilesComboBox.setObjectName(_fromUtf8("processedFilesComboBox"))
        self.gridLayout_5.addWidget(self.processedFilesComboBox, 1, 0, 1, 1)
        self.openExternallyButton = QtGui.QPushButton(self.resultsTab)
        self.openExternallyButton.setObjectName(_fromUtf8("openExternallyButton"))
        self.gridLayout_5.addWidget(self.openExternallyButton, 2, 0, 1, 1)
        spacerItem19 = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.gridLayout_5.addItem(spacerItem19, 3, 0, 1, 1)
        self.gridLayout_4.addLayout(self.gridLayout_5, 0, 0, 1, 1)
        self.resultsTreeWidget = QtGui.QTreeWidget(self.resultsTab)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.resultsTreeWidget.sizePolicy().hasHeightForWidth())
        self.resultsTreeWidget.setSizePolicy(sizePolicy)
        self.resultsTreeWidget.setObjectName(_fromUtf8("resultsTreeWidget"))
        self.resultsTreeWidget.headerItem().setText(1, _fromUtf8("MZ"))
        self.gridLayout_4.addWidget(self.resultsTreeWidget, 0, 1, 1, 1)
        self.gridLayout_3.addLayout(self.gridLayout_4, 3, 0, 1, 2)
        self.gridLayout_20 = QtGui.QGridLayout()
        self.gridLayout_20.setObjectName(_fromUtf8("gridLayout_20"))
        self.plotLabelledScanCheckBox = QtGui.QCheckBox(self.resultsTab)
        self.plotLabelledScanCheckBox.setChecked(True)
        self.plotLabelledScanCheckBox.setObjectName(_fromUtf8("plotLabelledScanCheckBox"))
        self.gridLayout_20.addWidget(self.plotLabelledScanCheckBox, 0, 1, 1, 1)
        spacerItem20 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.gridLayout_20.addItem(spacerItem20, 0, 3, 1, 1)
        self.plotCleanedScanCheckBox = QtGui.QCheckBox(self.resultsTab)
        self.plotCleanedScanCheckBox.setChecked(True)
        self.plotCleanedScanCheckBox.setObjectName(_fromUtf8("plotCleanedScanCheckBox"))
        self.gridLayout_20.addWidget(self.plotCleanedScanCheckBox, 0, 2, 1, 1)
        self.plotNativeScanCheckbox = QtGui.QCheckBox(self.resultsTab)
        self.plotNativeScanCheckbox.setChecked(True)
        self.plotNativeScanCheckbox.setObjectName(_fromUtf8("plotNativeScanCheckbox"))
        self.gridLayout_20.addWidget(self.plotNativeScanCheckbox, 0, 0, 1, 1)
        self.gridLayout_3.addLayout(self.gridLayout_20, 0, 0, 1, 2)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        self.visualizationWidget = QtGui.QWidget(self.resultsTab)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.visualizationWidget.sizePolicy().hasHeightForWidth())
        self.visualizationWidget.setSizePolicy(sizePolicy)
        self.visualizationWidget.setSizeIncrement(QtCore.QSize(3, 0))
        self.visualizationWidget.setObjectName(_fromUtf8("visualizationWidget"))
        self.horizontalLayout_6.addWidget(self.visualizationWidget)
        self.eicWidget = QtGui.QWidget(self.resultsTab)
        self.eicWidget.setSizeIncrement(QtCore.QSize(1, 0))
        self.eicWidget.setObjectName(_fromUtf8("eicWidget"))
        self.horizontalLayout_6.addWidget(self.eicWidget)
        self.horizontalLayout_6.setStretch(0, 2)
        self.horizontalLayout_6.setStretch(1, 1)
        self.gridLayout_3.addLayout(self.horizontalLayout_6, 2, 0, 1, 2)
        self.tabWidget.addTab(self.resultsTab, _fromUtf8(""))
        self.gridLayout.addWidget(self.tabWidget, 0, 0, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menuBar = QtGui.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 1391, 26))
        self.menuBar.setObjectName(_fromUtf8("menuBar"))
        self.menuFile = QtGui.QMenu(self.menuBar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuHelp = QtGui.QMenu(self.menuBar)
        self.menuHelp.setObjectName(_fromUtf8("menuHelp"))
        MainWindow.setMenuBar(self.menuBar)
        self.actionLoad_Settings = QtGui.QAction(MainWindow)
        self.actionLoad_Settings.setObjectName(_fromUtf8("actionLoad_Settings"))
        self.actionSave_Settings = QtGui.QAction(MainWindow)
        self.actionSave_Settings.setObjectName(_fromUtf8("actionSave_Settings"))
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionHelp = QtGui.QAction(MainWindow)
        self.actionHelp.setObjectName(_fromUtf8("actionHelp"))
        self.actionAbout = QtGui.QAction(MainWindow)
        self.actionAbout.setObjectName(_fromUtf8("actionAbout"))
        self.menuFile.addAction(self.actionLoad_Settings)
        self.menuFile.addAction(self.actionSave_Settings)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menuHelp.addAction(self.actionHelp)
        self.menuHelp.addSeparator()
        self.menuHelp.addAction(self.actionAbout)
        self.menuBar.addAction(self.menuFile.menuAction())
        self.menuBar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.label_4.setText(_translate("MainWindow", "<html><head/><body><p align=\"justify\">Please specify the performed LC-HRMSMS mesurements of M and M\'. Specify each target in single row using the button \'Add MS/MS target(s)\'. Load an LC-HRMSMS measurement multiple times if it contains more than one target. Select the appropriate scan events for the native and the labeled metabolite ions.</p></body></html>", None))
        self.label_3.setText(_translate("MainWindow", "<html><head/><body><p align=\"right\">Define MS/MS targets</p></body></html>", None))
        self.saveCompilation.setText(_translate("MainWindow", "Save compilation", None))
        self.loadCompilation.setText(_translate("MainWindow", "Load compilation", None))
        self.addMSMSTarget.setText(_translate("MainWindow", "Add MS/MS target(s)", None))
        self.deleteMSMSTarget.setText(_translate("MainWindow", "Delete current MS/MS target", None))
        self.label_5.setText(_translate("MainWindow", "<html><head/><body><p align=\"right\">Define/Edit group</p></body></html>", None))
        self.configAdductsButton.setText(_translate("MainWindow", "Configure adducts", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.inputTab), _translate("MainWindow", "Input", None))
        self.label_6.setText(_translate("MainWindow", "<html><head/><body><p align=\"right\">Processing settings</p></body></html>", None))
        self.label.setText(_translate("MainWindow", "<html><head/><body><p>Please specify the settings for processing the defined MS/MS targets. First, the two MS/MS spectra will be insepcted for corresponding native and U-<span style=\" vertical-align:super;\">13</span>C-labeled peaks. These peaks will be annotated with putative sum formulas using the determined number of carbon atoms. </p></body></html>", None))
        self.label_7.setText(_translate("MainWindow", "<html><head/><body><p align=\"right\">Run tasks</p></body></html>", None))
        self.label_8.setText(_translate("MainWindow", "CPU cores", None))
        self.keepCPUCoreUnusedCheckbox.setText(_translate("MainWindow", "Keep one core unused", None))
        self.startProcessingButton.setText(_translate("MainWindow", "Start", None))
        self.groupBox.setTitle(_translate("MainWindow", "Scan selection", None))
        self.label_9.setText(_translate("MainWindow", "Full scan EIC ppm", None))
        self.label_10.setText(_translate("MainWindow", "Intensity threshold", None))
        self.groupBox_3.setTitle(_translate("MainWindow", "Peak matching", None))
        self.label_12.setText(_translate("MainWindow", "Matching max. PPM error", None))
        self.label_13.setText(_translate("MainWindow", "Max. relative intensity error", None))
        self.groupBox_2.setTitle(_translate("MainWindow", "Scan pre-processing", None))
        self.label_11.setText(_translate("MainWindow", "Min. scaled peak intensity [%]", None))
        self.scalePrecursorMZ.setText(_translate("MainWindow", "Scale to precursor mz", None))
        self.groupBox_4.setTitle(_translate("MainWindow", "Fragment annotation", None))
        self.fragmentAnnotationButton.setText(_translate("MainWindow", "Used elements", None))
        self.label_14.setText(_translate("MainWindow", "Max. ppm error", None))
        self.useTracExtractAnnotation.setText(_translate("MainWindow", "Allow additional Cn (e.g. biotransformation products)", None))
        self.label_15.setText(_translate("MainWindow", "Min. number of labeling atoms", None))
        self.useZeroLabelingAtoms.setText(_translate("MainWindow", "Use zero labeling atoms (for biotransformation products)", None))
        self.groupBox_5.setTitle(_translate("MainWindow", "Parent annotation", None))
        self.applyParentFragmentConsistencyRuleCheckBox.setText(_translate("MainWindow", "Apply parent-fragment consistency rule", None))
        self.groupBox_6.setTitle(_translate("MainWindow", "Save results", None))
        self.saveResultsAsTSVCheckBox.setText(_translate("MainWindow", "Save as TSV", None))
        self.saveResultsAsPDFCheckBox.setText(_translate("MainWindow", "Save as PDF", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.processTab), _translate("MainWindow", "Process", None))
        self.label_2.setText(_translate("MainWindow", "Processed file", None))
        self.openExternallyButton.setText(_translate("MainWindow", "Open externally", None))
        self.resultsTreeWidget.headerItem().setText(0, _translate("MainWindow", "Target name", None))
        self.resultsTreeWidget.headerItem().setText(2, _translate("MainWindow", "Cn", None))
        self.resultsTreeWidget.headerItem().setText(3, _translate("MainWindow", "Sum formula", None))
        self.resultsTreeWidget.headerItem().setText(4, _translate("MainWindow", "Charge", None))
        self.resultsTreeWidget.headerItem().setText(5, _translate("MainWindow", "Scan num native isotopolog", None))
        self.resultsTreeWidget.headerItem().setText(6, _translate("MainWindow", "Scan num labeled isotopolog", None))
        self.resultsTreeWidget.headerItem().setText(7, _translate("MainWindow", "Full scan event", None))
        self.resultsTreeWidget.headerItem().setText(8, _translate("MainWindow", "MS2 scan event native", None))
        self.resultsTreeWidget.headerItem().setText(9, _translate("MainWindow", "MS2 scan event labeled isotopolog", None))
        self.plotLabelledScanCheckBox.setText(_translate("MainWindow", "Plot labeled MS2-scan", None))
        self.plotCleanedScanCheckBox.setText(_translate("MainWindow", "Plot cleaned scan(s)", None))
        self.plotNativeScanCheckbox.setText(_translate("MainWindow", "Plot native MS2-scan", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.resultsTab), _translate("MainWindow", "Sample results", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.menuHelp.setTitle(_translate("MainWindow", "Help", None))
        self.actionLoad_Settings.setText(_translate("MainWindow", "Load Settings", None))
        self.actionSave_Settings.setText(_translate("MainWindow", "Save Settings", None))
        self.actionExit.setText(_translate("MainWindow", "Exit", None))
        self.actionHelp.setText(_translate("MainWindow", "Help", None))
        self.actionAbout.setText(_translate("MainWindow", "About", None))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

