# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\FTICRwindow.ui'
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


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1369, 907)
        MainWindow.setAcceptDrops(True)
        MainWindow.setStyleSheet(_fromUtf8(""))
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout_2 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.stackedWidget = QtGui.QStackedWidget(self.centralwidget)
        self.stackedWidget.setObjectName(_fromUtf8("stackedWidget"))
        self.loadSamplesPage = QtWidgets.QWidget()
        self.loadSamplesPage.setObjectName(_fromUtf8("loadSamplesPage"))
        self.gridLayout_8 = QtWidgets.QGridLayout(self.loadSamplesPage)
        self.gridLayout_8.setObjectName(_fromUtf8("gridLayout_8"))
        self.gridLayout_7 = QtWidgets.QGridLayout()
        self.gridLayout_7.setObjectName(_fromUtf8("gridLayout_7"))
        spacerItem = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.gridLayout_7.addItem(spacerItem, 2, 2, 1, 1)
        self.loadSamplesTSVButton = QtWidgets.QPushButton(self.loadSamplesPage)
        font = QtGui.QFont()
        font.setPointSize(32)
        self.loadSamplesTSVButton.setFont(font)
        self.loadSamplesTSVButton.setObjectName(_fromUtf8("loadSamplesTSVButton"))
        self.gridLayout_7.addWidget(self.loadSamplesTSVButton, 2, 1, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.gridLayout_7.addItem(spacerItem1, 2, 0, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_7.addItem(spacerItem2, 1, 1, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.loadSamplesPage)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout_7.addWidget(self.label_2, 0, 0, 1, 3)
        self.loadSamplesMZXMLButton = QtWidgets.QPushButton(self.loadSamplesPage)
        font = QtGui.QFont()
        font.setPointSize(32)
        self.loadSamplesMZXMLButton.setFont(font)
        self.loadSamplesMZXMLButton.setObjectName(_fromUtf8("loadSamplesMZXMLButton"))
        self.gridLayout_7.addWidget(self.loadSamplesMZXMLButton, 3, 1, 1, 1)
        spacerItem3 = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_7.addItem(spacerItem3, 4, 1, 1, 1)
        self.gridLayout_8.addLayout(self.gridLayout_7, 0, 0, 1, 1)
        self.stackedWidget.addWidget(self.loadSamplesPage)
        self.processSamplesPage = QtWidgets.QWidget()
        self.processSamplesPage.setObjectName(_fromUtf8("processSamplesPage"))
        self.gridLayout_10 = QtWidgets.QGridLayout(self.processSamplesPage)
        self.gridLayout_10.setObjectName(_fromUtf8("gridLayout_10"))
        self.gridLayout_9 = QtWidgets.QGridLayout()
        self.gridLayout_9.setObjectName(_fromUtf8("gridLayout_9"))
        spacerItem4 = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_9.addItem(spacerItem4, 3, 1, 1, 1)
        spacerItem5 = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_9.addItem(spacerItem5, 1, 1, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.processSamplesPage)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout_9.addWidget(self.label_7, 0, 0, 1, 4)
        spacerItem6 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.gridLayout_9.addItem(spacerItem6, 2, 0, 1, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(self.processSamplesPage)
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.gridLayout_4 = QtWidgets.QGridLayout(self.groupBox_2)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.process = QtWidgets.QPushButton(self.groupBox_2)
        self.process.setObjectName(_fromUtf8("process"))
        self.gridLayout_4.addWidget(self.process, 9, 2, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_14.sizePolicy().hasHeightForWidth())
        self.label_14.setSizePolicy(sizePolicy)
        self.label_14.setObjectName(_fromUtf8("label_14"))
        self.gridLayout_4.addWidget(self.label_14, 2, 0, 1, 1)
        self.label_9 = QtWidgets.QLabel(self.groupBox_2)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.gridLayout_4.addWidget(self.label_9, 3, 0, 1, 1)
        self.intensityThreshold = QtWidgets.QDoubleSpinBox(self.groupBox_2)
        self.intensityThreshold.setMaximum(1e11)
        self.intensityThreshold.setProperty("value", 1000000.0)
        self.intensityThreshold.setObjectName(_fromUtf8("intensityThreshold"))
        self.gridLayout_4.addWidget(self.intensityThreshold, 2, 1, 1, 1)
        self.bracketingPPM = QtWidgets.QDoubleSpinBox(self.groupBox_2)
        self.bracketingPPM.setProperty("value", 2.0)
        self.bracketingPPM.setObjectName(_fromUtf8("bracketingPPM"))
        self.gridLayout_4.addWidget(self.bracketingPPM, 4, 1, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.groupBox_2)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.gridLayout_4.addWidget(self.label_10, 4, 0, 1, 1)
        self.checkIsotopologPattern = QtWidgets.QCheckBox(self.groupBox_2)
        self.checkIsotopologPattern.setChecked(True)
        self.checkIsotopologPattern.setObjectName(_fromUtf8("checkIsotopologPattern"))
        self.gridLayout_4.addWidget(self.checkIsotopologPattern, 5, 0, 1, 2)
        self.calibrate = QtWidgets.QCheckBox(self.groupBox_2)
        self.calibrate.setChecked(True)
        self.calibrate.setObjectName(_fromUtf8("calibrate"))
        self.gridLayout_4.addWidget(self.calibrate, 6, 0, 1, 2)
        self.matchPPM = QtWidgets.QDoubleSpinBox(self.groupBox_2)
        self.matchPPM.setProperty("value", 0.5)
        self.matchPPM.setObjectName(_fromUtf8("matchPPM"))
        self.gridLayout_4.addWidget(self.matchPPM, 3, 1, 1, 1)
        self.groupBox_7 = QtWidgets.QGroupBox(self.groupBox_2)
        self.groupBox_7.setObjectName(_fromUtf8("groupBox_7"))
        self.gridLayout_14 = QtWidgets.QGridLayout(self.groupBox_7)
        self.gridLayout_14.setObjectName(_fromUtf8("gridLayout_14"))
        self.label_12 = QtWidgets.QLabel(self.groupBox_7)
        self.label_12.setObjectName(_fromUtf8("label_12"))
        self.gridLayout_14.addWidget(self.label_12, 1, 0, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.groupBox_7)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.gridLayout_14.addWidget(self.label_11, 0, 0, 1, 1)
        self.enr12C = QtWidgets.QDoubleSpinBox(self.groupBox_7)
        self.enr12C.setDecimals(4)
        self.enr12C.setSingleStep(0.001)
        self.enr12C.setProperty("value", 0.9893)
        self.enr12C.setObjectName(_fromUtf8("enr12C"))
        self.gridLayout_14.addWidget(self.enr12C, 0, 1, 1, 1)
        self.enr13C = QtWidgets.QDoubleSpinBox(self.groupBox_7)
        self.enr13C.setDecimals(4)
        self.enr13C.setSingleStep(0.001)
        self.enr13C.setProperty("value", 0.99)
        self.enr13C.setObjectName(_fromUtf8("enr13C"))
        self.gridLayout_14.addWidget(self.enr13C, 1, 1, 1, 1)
        self.label_13 = QtWidgets.QLabel(self.groupBox_7)
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.gridLayout_14.addWidget(self.label_13, 2, 0, 1, 1)
        self.maxEnrDeviation = QtWidgets.QDoubleSpinBox(self.groupBox_7)
        self.maxEnrDeviation.setProperty("value", 0.15)
        self.maxEnrDeviation.setObjectName(_fromUtf8("maxEnrDeviation"))
        self.gridLayout_14.addWidget(self.maxEnrDeviation, 2, 1, 1, 1)
        self.gridLayout_4.addWidget(self.groupBox_7, 0, 0, 1, 2)
        self.groupBox_4 = QtWidgets.QGroupBox(self.groupBox_2)
        self.groupBox_4.setObjectName(_fromUtf8("groupBox_4"))
        self.gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_4)
        self.gridLayout_5.setObjectName(_fromUtf8("gridLayout_5"))
        self.groupBox_5 = QtWidgets.QGroupBox(self.groupBox_4)
        self.groupBox_5.setObjectName(_fromUtf8("groupBox_5"))
        self.gridLayout_6 = QtWidgets.QGridLayout(self.groupBox_5)
        self.gridLayout_6.setObjectName(_fromUtf8("gridLayout_6"))
        self.label_3 = QtWidgets.QLabel(self.groupBox_5)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_6.addWidget(self.label_3, 0, 0, 1, 1)
        self.SGRElements_min = QtWidgets.QLineEdit(self.groupBox_5)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.SGRElements_min.sizePolicy().hasHeightForWidth()
        )
        self.SGRElements_min.setSizePolicy(sizePolicy)
        self.SGRElements_min.setObjectName(_fromUtf8("SGRElements_min"))
        self.gridLayout_6.addWidget(self.SGRElements_min, 0, 1, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.groupBox_5)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_6.sizePolicy().hasHeightForWidth())
        self.label_6.setSizePolicy(sizePolicy)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout_6.addWidget(self.label_6, 1, 0, 1, 1)
        self.SGRElements_max = QtWidgets.QLineEdit(self.groupBox_5)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.SGRElements_max.sizePolicy().hasHeightForWidth()
        )
        self.SGRElements_max.setSizePolicy(sizePolicy)
        self.SGRElements_max.setObjectName(_fromUtf8("SGRElements_max"))
        self.gridLayout_6.addWidget(self.SGRElements_max, 1, 1, 1, 1)
        self.gridLayout_5.addWidget(self.groupBox_5, 2, 0, 1, 3)
        self.groupBox_6 = QtWidgets.QGroupBox(self.groupBox_4)
        self.groupBox_6.setObjectName(_fromUtf8("groupBox_6"))
        self.gridLayout_13 = QtWidgets.QGridLayout(self.groupBox_6)
        self.gridLayout_13.setObjectName(_fromUtf8("gridLayout_13"))
        self.adduct_MKm2Hm = QtWidgets.QCheckBox(self.groupBox_6)
        self.adduct_MKm2Hm.setObjectName(_fromUtf8("adduct_MKm2Hm"))
        self.gridLayout_13.addWidget(self.adduct_MKm2Hm, 1, 0, 1, 1)
        self.adduct_MNam2Hm = QtWidgets.QCheckBox(self.groupBox_6)
        self.adduct_MNam2Hm.setObjectName(_fromUtf8("adduct_MNam2Hm"))
        self.gridLayout_13.addWidget(self.adduct_MNam2Hm, 1, 1, 1, 1)
        self.adduct_MpClm = QtWidgets.QCheckBox(self.groupBox_6)
        self.adduct_MpClm.setChecked(True)
        self.adduct_MpClm.setObjectName(_fromUtf8("adduct_MpClm"))
        self.gridLayout_13.addWidget(self.adduct_MpClm, 0, 1, 1, 1)
        self.adduct_MmHm = QtWidgets.QCheckBox(self.groupBox_6)
        self.adduct_MmHm.setChecked(True)
        self.adduct_MmHm.setObjectName(_fromUtf8("adduct_MmHm"))
        self.gridLayout_13.addWidget(self.adduct_MmHm, 0, 0, 1, 1)
        self.adduct_MpFAmHm = QtWidgets.QCheckBox(self.groupBox_6)
        self.adduct_MpFAmHm.setObjectName(_fromUtf8("adduct_MpFAmHm"))
        self.gridLayout_13.addWidget(self.adduct_MpFAmHm, 2, 0, 1, 1)
        self.adduct_MpBrm = QtWidgets.QCheckBox(self.groupBox_6)
        self.adduct_MpBrm.setObjectName(_fromUtf8("adduct_MpBrm"))
        self.gridLayout_13.addWidget(self.adduct_MpBrm, 2, 1, 1, 1)
        self.gridLayout_5.addWidget(self.groupBox_6, 0, 0, 1, 3)
        self.useSevenGoldenRules = QtWidgets.QCheckBox(self.groupBox_4)
        self.useSevenGoldenRules.setChecked(True)
        self.useSevenGoldenRules.setObjectName(_fromUtf8("useSevenGoldenRules"))
        self.gridLayout_5.addWidget(self.useSevenGoldenRules, 5, 0, 1, 3)
        self.label = QtWidgets.QLabel(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout_5.addWidget(self.label, 1, 0, 1, 2)
        self.annotationPPM = QtWidgets.QDoubleSpinBox(self.groupBox_4)
        self.annotationPPM.setProperty("value", 0.5)
        self.annotationPPM.setObjectName(_fromUtf8("annotationPPM"))
        self.gridLayout_5.addWidget(self.annotationPPM, 1, 2, 1, 1)
        self.groupBox_8 = QtWidgets.QGroupBox(self.groupBox_4)
        self.groupBox_8.setObjectName(_fromUtf8("groupBox_8"))
        self.gridLayout_15 = QtWidgets.QGridLayout(self.groupBox_8)
        self.gridLayout_15.setObjectName(_fromUtf8("gridLayout_15"))
        self.DBs = QtGui.QListView(self.groupBox_8)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Ignored
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.DBs.sizePolicy().hasHeightForWidth())
        self.DBs.setSizePolicy(sizePolicy)
        self.DBs.setMinimumSize(QtCore.QSize(0, 40))
        self.DBs.setObjectName(_fromUtf8("DBs"))
        self.gridLayout_15.addWidget(self.DBs, 0, 0, 1, 1)
        self.addDB = QtWidgets.QPushButton(self.groupBox_8)
        self.addDB.setObjectName(_fromUtf8("addDB"))
        self.gridLayout_15.addWidget(self.addDB, 1, 0, 1, 1)
        self.gridLayout_5.addWidget(self.groupBox_8, 6, 0, 1, 3)
        self.gridLayout_4.addWidget(self.groupBox_4, 0, 2, 9, 1)
        self.label_15 = QtWidgets.QLabel(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_15.sizePolicy().hasHeightForWidth())
        self.label_15.setSizePolicy(sizePolicy)
        self.label_15.setObjectName(_fromUtf8("label_15"))
        self.gridLayout_4.addWidget(self.label_15, 1, 0, 1, 1)
        self.cAtomsCounts = QtWidgets.QLineEdit(self.groupBox_2)
        self.cAtomsCounts.setObjectName(_fromUtf8("cAtomsCounts"))
        self.gridLayout_4.addWidget(self.cAtomsCounts, 1, 1, 1, 1)
        self.gridLayout_9.addWidget(self.groupBox_2, 2, 1, 1, 1)
        spacerItem7 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.gridLayout_9.addItem(spacerItem7, 2, 3, 1, 1)
        self.gridLayout_10.addLayout(self.gridLayout_9, 0, 0, 1, 1)
        self.stackedWidget.addWidget(self.processSamplesPage)
        self.showResultsPage = QtWidgets.QWidget()
        self.showResultsPage.setObjectName(_fromUtf8("showResultsPage"))
        self.gridLayout_12 = QtWidgets.QGridLayout(self.showResultsPage)
        self.gridLayout_12.setObjectName(_fromUtf8("gridLayout_12"))
        self.groupBox_3 = QtWidgets.QGroupBox(self.showResultsPage)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_3.sizePolicy().hasHeightForWidth())
        self.groupBox_3.setSizePolicy(sizePolicy)
        self.groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox_3)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.result_MassSpectrum = QtWidgets.QWidget(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.result_MassSpectrum.sizePolicy().hasHeightForWidth()
        )
        self.result_MassSpectrum.setSizePolicy(sizePolicy)
        self.result_MassSpectrum.setObjectName(_fromUtf8("result_MassSpectrum"))
        self.gridLayout_3.addWidget(self.result_MassSpectrum, 1, 0, 1, 2)
        self.gridLayout_12.addWidget(self.groupBox_3, 1, 1, 1, 1)
        self.gridLayout_11 = QtWidgets.QGridLayout()
        self.gridLayout_11.setObjectName(_fromUtf8("gridLayout_11"))
        self.groupBox = QtWidgets.QGroupBox(self.showResultsPage)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.loadedSamples = QtGui.QListView(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.loadedSamples.sizePolicy().hasHeightForWidth()
        )
        self.loadedSamples.setSizePolicy(sizePolicy)
        self.loadedSamples.setMinimumSize(QtCore.QSize(0, 0))
        self.loadedSamples.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.loadedSamples.setObjectName(_fromUtf8("loadedSamples"))
        self.gridLayout.addWidget(self.loadedSamples, 1, 0, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_4.sizePolicy().hasHeightForWidth())
        self.label_4.setSizePolicy(sizePolicy)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout.addWidget(self.label_4, 0, 0, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout.addWidget(self.label_5, 3, 0, 1, 1)
        self.results = QtWidgets.QTreeWidget(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.results.sizePolicy().hasHeightForWidth())
        self.results.setSizePolicy(sizePolicy)
        self.results.setMinimumSize(QtCore.QSize(350, 0))
        self.results.setDragDropMode(QtGui.QAbstractItemView.DropOnly)
        self.results.setObjectName(_fromUtf8("results"))
        self.gridLayout.addWidget(self.results, 4, 0, 1, 1)
        self.gridLayout_11.addWidget(self.groupBox, 0, 0, 1, 1)
        self.gridLayout_12.addLayout(self.gridLayout_11, 1, 0, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.showResultsPage)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridLayout_12.addWidget(self.label_8, 0, 0, 1, 2)
        self.stackedWidget.addWidget(self.showResultsPage)
        self.gridLayout_2.addWidget(self.stackedWidget, 2, 0, 1, 2)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1369, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionLoad_samples = QtGui.QAction(MainWindow)
        self.actionLoad_samples.setObjectName(_fromUtf8("actionLoad_samples"))
        self.actionSave_results = QtGui.QAction(MainWindow)
        self.actionSave_results.setObjectName(_fromUtf8("actionSave_results"))
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionLoad_results = QtGui.QAction(MainWindow)
        self.actionLoad_results.setObjectName(_fromUtf8("actionLoad_results"))
        self.actionNew_processing = QtGui.QAction(MainWindow)
        self.actionNew_processing.setObjectName(_fromUtf8("actionNew_processing"))
        self.menuFile.addAction(self.actionNew_processing)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menubar.addAction(self.menuFile.menuAction())

        self.retranslateUi(MainWindow)
        self.stackedWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(
            _translate("MainWindow", "MetExtract II - FT-ICR module", None)
        )
        self.loadSamplesTSVButton.setText(
            _translate("MainWindow", "Load samples (.tsv format)", None)
        )
        self.label_2.setText(
            _translate(
                "MainWindow",
                '<html><head/><body><p align="right">\n'
                "\n"
                '<span style=" font-size:18pt;">Load sample files</span></p>\n'
                "\n"
                '<p align="right">First load the sample files in the tab-separated-values (.tsv) format.<br>\n'
                "Two columns, one specifying the m/z value and one the observed intensity, must be present.<br>\n"
                "The respective columns can be specified later in a separate dialog.</p>\n"
                "\n"
                "</body></html>",
                None,
            )
        )
        self.loadSamplesMZXMLButton.setText(
            _translate("MainWindow", "Load samples (.mzXML format)", None)
        )
        self.label_7.setText(
            _translate(
                "MainWindow",
                '<html><head/><body><p align="right">\n'
                "\n"
                '<span style=" font-size:18pt;">Processing parameters</span></p>\n'
                "\n"
                '<p align="right">Please specify all necessary processing parameters for the loaded FT-ICR MS scans.</p>\n'
                "\n"
                "</body></html>",
                None,
            )
        )
        self.groupBox_2.setTitle(_translate("MainWindow", "Process samples", None))
        self.process.setText(_translate("MainWindow", "Process", None))
        self.label_14.setText(_translate("MainWindow", "Intensity threshold", None))
        self.label_9.setText(_translate("MainWindow", "Match ppm", None))
        self.label_10.setText(_translate("MainWindow", "Bracketing ppm", None))
        self.checkIsotopologPattern.setText(
            _translate("MainWindow", "M+1 and M'-1 must be present", None)
        )
        self.calibrate.setText(
            _translate(
                "MainWindow",
                "Calibrate with average MZ error based\n"
                "on generated, unique sum formulas",
                None,
            )
        )
        self.groupBox_7.setTitle(_translate("MainWindow", "Labeling enrichment", None))
        self.label_12.setText(_translate("MainWindow", "<sup>13</sup>C-labeled", None))
        self.label_11.setText(_translate("MainWindow", "Native", None))
        self.label_13.setText(_translate("MainWindow", "Maximum deviation +-", None))
        self.groupBox_4.setTitle(
            _translate("MainWindow", "Sum formula generation", None)
        )
        self.groupBox_5.setTitle(
            _translate("MainWindow", "Elements for sum formula generation", None)
        )
        self.label_3.setText(_translate("MainWindow", "Min", None))
        self.SGRElements_min.setText(_translate("MainWindow", "CH", None))
        self.label_6.setText(_translate("MainWindow", "Max", None))
        self.SGRElements_max.setText(
            _translate("MainWindow", "CH1000N10O100S10P10", None)
        )
        self.groupBox_6.setTitle(_translate("MainWindow", "Used adducts", None))
        self.adduct_MKm2Hm.setText(_translate("MainWindow", "[M+K-2H]-", None))
        self.adduct_MNam2Hm.setText(_translate("MainWindow", "[M+Na-2H]-", None))
        self.adduct_MpClm.setText(_translate("MainWindow", "[M+Cl]-", None))
        self.adduct_MmHm.setText(_translate("MainWindow", "[M-H]-", None))
        self.adduct_MpFAmHm.setText(_translate("MainWindow", "[M+FA-H]-", None))
        self.adduct_MpBrm.setText(_translate("MainWindow", "[M+Br]-", None))
        self.useSevenGoldenRules.setText(
            _translate("MainWindow", "Use Seven Golden Rules (Kind et al. 2007)", None)
        )
        self.label.setText(_translate("MainWindow", "Annotation ppm", None))
        self.groupBox_8.setTitle(_translate("MainWindow", "Sumformula databases", None))
        self.addDB.setText(_translate("MainWindow", "Add database", None))
        self.label_15.setText(
            _translate("MainWindow", "Carbon atoms to search for", None)
        )
        self.cAtomsCounts.setText(_translate("MainWindow", "3-60", None))
        self.groupBox_3.setTitle(_translate("MainWindow", "Results", None))
        self.groupBox.setTitle(
            _translate("MainWindow", "Loaded samples and results", None)
        )
        self.label_4.setText(_translate("MainWindow", "Samples", None))
        self.label_5.setText(_translate("MainWindow", "Results", None))
        self.results.headerItem().setText(0, _translate("MainWindow", "MZ", None))
        self.results.headerItem().setText(1, _translate("MainWindow", "Cn", None))
        self.results.headerItem().setText(
            2, _translate("MainWindow", "Sumformula", None)
        )
        self.label_8.setText(
            _translate(
                "MainWindow",
                '<html><head/><body><p align="right">\n'
                "\n"
                '<span style=" font-size:18pt;">Detected signal pairs</span></p>\n'
                "\n"
                '<p align="right">The data processing results can be quickly visualized here.</p>\n'
                "\n"
                "</body></html>",
                None,
            )
        )
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.actionLoad_samples.setText(_translate("MainWindow", "Load samples", None))
        self.actionSave_results.setText(_translate("MainWindow", "Save results", None))
        self.actionExit.setText(_translate("MainWindow", "Exit", None))
        self.actionLoad_results.setText(_translate("MainWindow", "Load results", None))
        self.actionNew_processing.setText(
            _translate("MainWindow", "New processing", None)
        )


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())
