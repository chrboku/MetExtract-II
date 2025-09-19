# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\mePyGuis\guis\ModuleSelectionWindow.ui'
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
        MainWindow.resize(638, 698)
        icon = QtGui.QIcon()
        icon.addPixmap(
            QtGui.QPixmap(_fromUtf8(":/MEIcon/resources/MEIcon.ico")),
            QtGui.QIcon.Normal,
            QtGui.QIcon.Off,
        )
        MainWindow.setWindowIcon(icon)
        MainWindow.setStyleSheet(_fromUtf8("background-color: rgb(255, 255, 255);"))
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout_2 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 6, 1, 1, 1)
        self.tracExtractIcon = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(90)
        sizePolicy.setVerticalStretch(90)
        sizePolicy.setHeightForWidth(
            self.tracExtractIcon.sizePolicy().hasHeightForWidth()
        )
        self.tracExtractIcon.setSizePolicy(sizePolicy)
        self.tracExtractIcon.setMinimumSize(QtCore.QSize(90, 90))
        self.tracExtractIcon.setStyleSheet(
            _fromUtf8("background-image: url(:/TracExtract/resources/TracExtract.png);")
        )
        self.tracExtractIcon.setText(_fromUtf8(""))
        self.tracExtractIcon.setFlat(True)
        self.tracExtractIcon.setObjectName(_fromUtf8("tracExtractIcon"))
        self.gridLayout.addWidget(self.tracExtractIcon, 2, 0, 1, 1)
        self.fragExtractIcon = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(90)
        sizePolicy.setVerticalStretch(90)
        sizePolicy.setHeightForWidth(
            self.fragExtractIcon.sizePolicy().hasHeightForWidth()
        )
        self.fragExtractIcon.setSizePolicy(sizePolicy)
        self.fragExtractIcon.setMinimumSize(QtCore.QSize(90, 90))
        self.fragExtractIcon.setStyleSheet(
            _fromUtf8("background-image: url(:/FragExtract/resources/FragExtract.png);")
        )
        self.fragExtractIcon.setText(_fromUtf8(""))
        self.fragExtractIcon.setFlat(True)
        self.fragExtractIcon.setObjectName(_fromUtf8("fragExtractIcon"))
        self.gridLayout.addWidget(self.fragExtractIcon, 4, 0, 1, 1)
        self.fragExtractLabel = QtWidgets.QLabel(self.centralwidget)
        self.fragExtractLabel.setObjectName(_fromUtf8("fragExtractLabel"))
        self.gridLayout.addWidget(self.fragExtractLabel, 4, 1, 1, 1)
        self.documentationLabel = QtWidgets.QLabel(self.centralwidget)
        self.documentationLabel.setObjectName(_fromUtf8("documentationLabel"))
        self.gridLayout.addWidget(self.documentationLabel, 11, 1, 1, 1)
        self.line_3 = QtWidgets.QFrame(self.centralwidget)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setObjectName(_fromUtf8("line_3"))
        self.gridLayout.addWidget(self.line_3, 9, 0, 1, 2)
        self.allExtractLabel = QtWidgets.QLabel(self.centralwidget)
        self.allExtractLabel.setStyleSheet(_fromUtf8(""))
        self.allExtractLabel.setObjectName(_fromUtf8("allExtractLabel"))
        self.gridLayout.addWidget(self.allExtractLabel, 0, 1, 1, 1)
        self.line = QtWidgets.QFrame(self.centralwidget)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout.addWidget(self.line, 1, 0, 1, 2)
        self.documentationIcon = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(90)
        sizePolicy.setVerticalStretch(90)
        sizePolicy.setHeightForWidth(
            self.documentationIcon.sizePolicy().hasHeightForWidth()
        )
        self.documentationIcon.setSizePolicy(sizePolicy)
        self.documentationIcon.setMinimumSize(QtCore.QSize(90, 90))
        self.documentationIcon.setStyleSheet(
            _fromUtf8(
                "background-image: url(:/Documentation/resources/Documentation.png);"
            )
        )
        self.documentationIcon.setText(_fromUtf8(""))
        self.documentationIcon.setFlat(True)
        self.documentationIcon.setObjectName(_fromUtf8("documentationIcon"))
        self.gridLayout.addWidget(self.documentationIcon, 11, 0, 1, 1)
        self.tracExtractLabel = QtWidgets.QLabel(self.centralwidget)
        self.tracExtractLabel.setObjectName(_fromUtf8("tracExtractLabel"))
        self.gridLayout.addWidget(self.tracExtractLabel, 2, 1, 1, 1)
        self.allExtractIcon = QtWidgets.QPushButton(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(90)
        sizePolicy.setVerticalStretch(90)
        sizePolicy.setHeightForWidth(
            self.allExtractIcon.sizePolicy().hasHeightForWidth()
        )
        self.allExtractIcon.setSizePolicy(sizePolicy)
        self.allExtractIcon.setMinimumSize(QtCore.QSize(90, 90))
        self.allExtractIcon.setStyleSheet(
            _fromUtf8("background-image: url(:/AllExtract/resources/AllExtract.png);")
        )
        self.allExtractIcon.setText(_fromUtf8(""))
        self.allExtractIcon.setFlat(True)
        self.allExtractIcon.setObjectName(_fromUtf8("allExtractIcon"))
        self.gridLayout.addWidget(self.allExtractIcon, 0, 0, 1, 1)
        self.line_4 = QtWidgets.QFrame(self.centralwidget)
        self.line_4.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName(_fromUtf8("line_4"))
        self.gridLayout.addWidget(self.line_4, 10, 0, 1, 2)
        self.combineResultsButton = QtWidgets.QPushButton(self.centralwidget)
        self.combineResultsButton.setMinimumSize(QtCore.QSize(90, 90))
        self.combineResultsButton.setStyleSheet(
            _fromUtf8(
                "background-image: url(:/combineResults/resources/combineResults.png);"
            )
        )
        self.combineResultsButton.setText(_fromUtf8(""))
        self.combineResultsButton.setFlat(True)
        self.combineResultsButton.setObjectName(_fromUtf8("combineResultsButton"))
        self.gridLayout.addWidget(self.combineResultsButton, 6, 0, 1, 1)
        self.line_5 = QtWidgets.QFrame(self.centralwidget)
        self.line_5.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_5.setObjectName(_fromUtf8("line_5"))
        self.gridLayout.addWidget(self.line_5, 5, 0, 1, 2)
        self.line_2 = QtWidgets.QFrame(self.centralwidget)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setObjectName(_fromUtf8("line_2"))
        self.gridLayout.addWidget(self.line_2, 3, 0, 1, 2)
        self.fticrExtractIcon = QtWidgets.QPushButton(self.centralwidget)
        self.fticrExtractIcon.setMinimumSize(QtCore.QSize(90, 90))
        self.fticrExtractIcon.setStyleSheet(
            _fromUtf8(
                "background-image: url(:/FTICRExtract/resources/FTICRExtract.png);"
            )
        )
        self.fticrExtractIcon.setText(_fromUtf8(""))
        self.fticrExtractIcon.setFlat(True)
        self.fticrExtractIcon.setObjectName(_fromUtf8("fticrExtractIcon"))
        self.gridLayout.addWidget(self.fticrExtractIcon, 8, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 8, 1, 1, 1)
        self.line_6 = QtWidgets.QFrame(self.centralwidget)
        self.line_6.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_6.setObjectName(_fromUtf8("line_6"))
        self.gridLayout.addWidget(self.line_6, 7, 0, 1, 2)
        self.gridLayout_2.addLayout(self.gridLayout, 1, 1, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_2.addItem(spacerItem, 0, 1, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.gridLayout_2.addItem(spacerItem1, 2, 1, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(
            40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.gridLayout_2.addItem(spacerItem2, 1, 0, 1, 1)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        spacerItem3 = QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.verticalLayout.addItem(spacerItem3)
        self.version = QtWidgets.QLabel(self.centralwidget)
        self.version.setStyleSheet(_fromUtf8("color: slategrey;\n"))
        self.version.setText(_fromUtf8(""))
        self.version.setAlignment(
            QtCore.Qt.AlignRight | QtCore.Qt.AlignTrailing | QtCore.Qt.AlignVCenter
        )
        self.version.setObjectName(_fromUtf8("version"))
        self.verticalLayout.addWidget(self.version)
        self.gridLayout_2.addLayout(self.verticalLayout, 2, 4, 1, 1)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem4 = QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.horizontalLayout.addItem(spacerItem4)
        self.label = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setMinimumSize(QtCore.QSize(80, 80))
        self.label.setStyleSheet(
            _fromUtf8("image: url(:/MEIcon_Large/resources/MEIcon_Large.png);")
        )
        self.label.setText(_fromUtf8(""))
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout.addWidget(self.label)
        self.verticalLayout_3.addLayout(self.horizontalLayout)
        spacerItem5 = QtWidgets.QSpacerItem(
            20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
        )
        self.verticalLayout_3.addItem(spacerItem5)
        spacerItem6 = QtWidgets.QSpacerItem(
            0, 0, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
        )
        self.verticalLayout_3.addItem(spacerItem6)
        self.gridLayout_2.addLayout(self.verticalLayout_3, 0, 4, 2, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menuBar = QtWidgets.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 638, 21))
        self.menuBar.setObjectName(_fromUtf8("menuBar"))
        self.menuTools = QtWidgets.QMenu(self.menuBar)
        self.menuTools.setObjectName(_fromUtf8("menuTools"))
        MainWindow.setMenuBar(self.menuBar)
        self.actionCalculate_isotopic_enrichment = QtGui.QAction(MainWindow)
        self.actionCalculate_isotopic_enrichment.setObjectName(
            _fromUtf8("actionCalculate_isotopic_enrichment")
        )
        self.menuTools.addAction(self.actionCalculate_isotopic_enrichment)
        self.menuBar.addAction(self.menuTools.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.label_2.setText(
            _translate(
                "MainWindow", "Combine results\n(AllExtract and TracExtract)\n", None
            )
        )
        self.fragExtractLabel.setText(
            _translate("MainWindow", "FragExtract\n(LC-HRMS)\n", None)
        )
        self.documentationLabel.setText(_translate("MainWindow", "Documentation", None))
        self.allExtractLabel.setText(
            _translate("MainWindow", "AllExtract\n(LC-HRMS)", None)
        )
        self.tracExtractLabel.setText(
            _translate("MainWindow", "TracExtract\n(LC-HRMS)", None)
        )
        self.label_3.setText(
            _translate("MainWindow", "FTICRExtract\n(FT-ICR-MS)\n", None)
        )
        self.menuTools.setTitle(_translate("MainWindow", "Tools", None))
        self.actionCalculate_isotopic_enrichment.setText(
            _translate("MainWindow", "Calculate isotopic enrichment", None)
        )


import resources_rc

if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())
