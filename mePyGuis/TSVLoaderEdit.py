import csv
import os

from PyQt4 import QtGui, QtCore

from TSVLoaderEditor import Ui_Dialog


class TSVLoaderEdit(QtGui.QDialog, Ui_Dialog):
    def __init__(self, parent=None, initDir=None, mapping={}, order=[]):
        QtGui.QDialog.__init__(self, parent)
        self.setWindowTitle("Adducts editor")
        self.setupUi(self)

        self.discardButton.clicked.connect(self.dialogCan)
        self.acceptButton.clicked.connect(self.dialogFin)
        self.loadTable.clicked.connect(self.loadFiles)

        self.mapping = mapping
        self.order = order
        if len(self.order) == 0:
            self.order = self.mapping.keys()
        else:
            assert len(self.order) == len(self.mapping)
            for ord in self.order:
                assert ord in self.mapping

        self.lastOpenDir = "."
        if os.environ.has_key('USERPROFILE'):
            self.lastOpenDir = os.getenv('USERPROFILE')
        elif os.environ.has_key('HOME'):
            self.lastOpenDir = os.getenv('HOME')
        if initDir is not None:
            self.lastOpenDir = initDir

        self.selectedFiles = None
        self.mappingComboBoxes = {}

    def loadFiles(self):

        scrollAreaWidgetContents = QtGui.QWidget()
        scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 291, 370))
        gridLayout_3 = QtGui.QGridLayout(scrollAreaWidgetContents)
        verticalLayout_3 = QtGui.QVBoxLayout()
        spacerItem6 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        verticalLayout_3.addItem(spacerItem6)

        openFiles = QtGui.QFileDialog.getOpenFileNames(caption="Select TSV/CSV file", directory=self.lastOpenDir,
                                                     filter="TSV file (*.tsv);;CSV file (*.csv);;All files(*.*)")

        if len(openFiles) > 0:
            self.selectedFiles = [str(f) for f in openFiles]
            self.mappingComboBoxes = {}

            horizontalLayout_5 = QtGui.QHBoxLayout()
            l1 = QtGui.QLabel(self.scrollAreaWidgetContents)
            sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(l1.sizePolicy().hasHeightForWidth())
            l1.setSizePolicy(sizePolicy)

            mFile=""
            if len(self.selectedFiles) == 1:
                mFile = self.selectedFiles[0]
                if len(mFile) > 26:
                    mFile = mFile[:3] + "..." + mFile[(len(mFile) - 23):]
            else:
                mFile = "Selected %d files"%len(self.selectedFiles)
            l1.setText("<html><head/><body><p>Selected file: %s</p></body></html>" % mFile)
            horizontalLayout_5.addWidget(l1)
            verticalLayout_3.addLayout(horizontalLayout_5)

            splitChar = '\t'
            if self.selectedFiles[0].lower().endswith(".csv"):
                splitChar = ';'

            headers = None
            for opFile in self.selectedFiles:
                print opFile
                with open(opFile, 'rb') as ofin:
                    tsvin = csv.reader(ofin, delimiter=splitChar)

                    i = 0
                    for row in tsvin:
                        if i == 0:
                            if headers == None:
                                headers = row
                            else:
                                headers = list(set(headers) & set(row))
                        else:
                            break
                        i = i + 1

            for map in self.order:

                horizontalLayout_5 = QtGui.QHBoxLayout()

                l1 = QtGui.QLabel(self.scrollAreaWidgetContents)
                sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
                sizePolicy.setHorizontalStretch(0)
                sizePolicy.setVerticalStretch(0)
                sizePolicy.setHeightForWidth(l1.sizePolicy().hasHeightForWidth())
                l1.setSizePolicy(sizePolicy)
                l1.setText(
                    "<html><head/><body><p><span style=\" text-decoration: underline;\">%s</span>:</p></body></html>" % map)

                horizontalLayout_5.addWidget(l1)
                c1 = QtGui.QComboBox(self.scrollAreaWidgetContents)
                self.mappingComboBoxes[map] = c1

                for header in headers:
                    c1.addItem(header)
                if self.mapping[map] in headers:
                    c1.setCurrentIndex(self.findMapInHeaders(self.mapping[map], headers))

                horizontalLayout_5.addWidget(c1)
                verticalLayout_3.addLayout(horizontalLayout_5)
        else:
            self.selectedFiles = []
            self.mappingComboBoxes = {}

            horizontalLayout_5 = QtGui.QHBoxLayout()

            l1 = QtGui.QLabel(self.scrollAreaWidgetContents)
            sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(l1.sizePolicy().hasHeightForWidth())
            l1.setSizePolicy(sizePolicy)
            l1.setText(
                "<html><head/><body><p>Please open TSV or CSV file using the 'Load table' button</p></body></html>")

            horizontalLayout_5.addWidget(l1)
            verticalLayout_3.addLayout(horizontalLayout_5)

        spacerItem7 = QtGui.QSpacerItem(0, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        verticalLayout_3.addItem(spacerItem7)
        gridLayout_3.addLayout(verticalLayout_3, 0, 0, 1, 1)
        self.scrollAreaWidgetContents = scrollAreaWidgetContents
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)

    def findMapInHeaders(self, mapping, headers):
        pos = -1
        if mapping in headers:
            i = 0
            for header in headers:
                if mapping == header:
                    pos = i
                i = i + 1
        return pos

    def setTitle(self, title):
        self.titleText.setText(title)

    def setDescription(self, desc):
        self.descriptionText.setText(desc)

    def getUserSelectedMapping(self):
        if len(self.mappingComboBoxes) == 0:
            raise Exception("No valid mapping was defined by the user")
        else:
            mapped = {}
            for map in self.mappingComboBoxes:
                mapped[map] = str(self.mappingComboBoxes[map].currentText())
        return mapped

    def getUserSelectedFile(self):
        return self.selectedFiles

    def dialogCan(self):
        self.reject()

    def dialogFin(self):
        if len(self.selectedFiles)==0:
            self.reject()
        else:
            self.accept()

    def executeDialog(self):
        self.loadFiles()
        x = self.exec_()
        return x


if __name__ == "__main__":
    import sys

    app = QtGui.QApplication(sys.argv)

    dialog = TSVLoaderEdit(
        mapping={"ID": "Id", "mz": "mz", "Cn": "C_Count", "Z": "Z", "RT": "RT_min", "MS scan event": "ScanEvent"},
        order=["ID", "mz", "Cn", "Z", "RT", "MS scan event"])

    dialog.setTitle("Select results file")
    dialog.setDescription("Select results file for some example")
    dialog.setWindowTitle("Select results file")
    dialog.resize(560, 80)

    if dialog.executeDialog() == QtGui.QDialog.Accepted:
        print dialog.getUserSelectedFile(), "\n", dialog.getUserSelectedMapping()

    sys.exit(app.exec_())

