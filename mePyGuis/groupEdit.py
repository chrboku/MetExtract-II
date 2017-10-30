import os

from PyQt4 import QtGui, QtCore

from mePyGuis.groupEditor import Ui_GroupEditor


class groupEdit(QtGui.QDialog, Ui_GroupEditor):
    def __init__(self, parent=None, initDir=None, colors=["Red", "Blue", "Green"], activeColor=0):
        self.groupfiles = []
        if initDir is not None:
            self.initDir = initDir
        else:
            self.initDir = "."
            if os.environ.has_key('USERPROFILE'):
                self.initDir = os.getenv('USERPROFILE')
            elif os.environ.has_key('HOME'):
                self.initDir = os.getenv('HOME')

        QtGui.QDialog.__init__(self, parent)
        self.setupUi(self)
        self.setWindowTitle("Group editor")

        self.dialogFinished.clicked.connect(self.dialogFin)
        self.dialogCanceled.clicked.connect(self.dialogCan)
        self.addFiles.clicked.connect(self.selectFiles)

        self.addFolder.setVisible(False)
        self.addFolder.clicked.connect(self.selectFolder)
        self.removeSelected.clicked.connect(self.removeSel)

        self.colors=QtCore.QStringList()
        for i in colors:
            self.colors.append(i)
        self.colorsComboBox.addItems(self.colors)
        self.colorsComboBox.setCurrentIndex(activeColor)
        self.col = self.colors[self.colorsComboBox.currentIndex()]

        self.dialogFinished.setFocus(True)

        self.groupName.setValidator(QtGui.QRegExpValidator(QtCore.QRegExp("[0-9a-zA-Z_]*")))
        self.setAcceptDrops(True)

    def removeSel(self):
        todel = []
        for selectedIndex in self.groupFiles.selectedIndexes():
            todel.append(selectedIndex.row())
        todel.sort()
        todel.reverse()

        for ind in todel:
            self.groupFiles.takeItem(ind)
            self.groupfiles.pop(ind)

    def dialogCan(self):
        self.reject()

    def dialogFin(self):
        if len(self.getGroupName()) < 1:
            QtGui.QMessageBox.information(None, "Group name", "Please specify a group name", QtGui.QMessageBox.Ok)
            return
        if len(self.getGroupFiles()) < 1:
            QtGui.QMessageBox.information(None, "Group files", "Please load one or more files into this group",
                                          QtGui.QMessageBox.Ok)
            return

        self.accept()

    def getGroupName(self):
        return str(self.groupName.text())

    def getGroupFiles(self):
        return [str(t) for t in self.groupfiles]

    def getMinimumGroupFound(self):
        return int(self.groupMinimumFound.value())

    def getOmitFeatures(self):
        return bool(self.omitFeatures.checkState() == QtCore.Qt.Checked)

    def getUseForMetaboliteGrouping(self):
        return bool(self.useForMetaboliteGrouping.checkState() == QtCore.Qt.Checked)

    def getRemoveAsFalsePositive(self):
        return bool(self.removeAsFalsePositive.checkState() == QtCore.Qt.Checked)

    def getOpenDir(self):
        return self.initDir

    def getGroupColor(self):
        return self.colors[self.colorsComboBox.currentIndex()]

    def selectFiles(self):
        filenames = QtGui.QFileDialog.getOpenFileNames(self, caption="Select group file(s)", directory=self.initDir,
                                                       filter="mzXML (*.mzxml);;mzML (*.mzml);;group file (*.grp);;All files (*.*)")
        for filename in filenames:
            self.groupFiles.addItem(filename.replace("\\", "/"))
            self.groupfiles.append(filename)
        if len(filenames) > 0:
            self.initDir = str(filenames[0]).replace("\\", "/")
            self.initDir = self.initDir[:self.initDir.rfind("/")]

    def selectFolder(self):
        foldername = str(QtGui.QFileDialog.getExistingDirectory(self, caption="Select folder containing mzxml files"),
                         directory=self.initDir)
        if len(foldername) > 0:
            self.initDir = foldername.replace("\\", "/")
            for root, dirs, files in os.walk(foldername):
                for file in files:
                    if file.lower().endswith(".mzxml") or file.lower().endswith(".mzml"):
                        self.groupFiles.addItem(root + "/" + file)
                        self.groupfiles.append(root + "/" + file)


    def executeDialog(self, groupName="", groupfiles=[], minimumGroupFound=1, omitFeatures=True, useForMetaboliteGrouping=True, removeAsFalsePositive=False, activeColor="Red"):

        for file in groupfiles:
            self.groupFiles.addItem(file)
            self.groupfiles.append(file)

            self.initDir = str(file).replace("\\", "/")
            self.initDir = self.initDir[:self.initDir.rfind("/")]

            colInd=-1
            for i, col in enumerate(self.colors):
                if col==activeColor:
                    colInd=i
            if colInd==-1:
                self.colors.append(activeColor)
                self.colorsComboBox.clear()
                self.colorsComboBox.addItems(self.colors)
                colInd=len(self.colors)-1
            self.colorsComboBox.setCurrentIndex(colInd)

        self.groupName.setText(groupName)
        self.groupMinimumFound.setValue(minimumGroupFound)
        if omitFeatures:
            self.omitFeatures.setCheckState(QtCore.Qt.Checked)
        else:
            self.omitFeatures.setCheckState(QtCore.Qt.Unchecked)
        if useForMetaboliteGrouping:
            self.useForMetaboliteGrouping.setCheckState(QtCore.Qt.Checked)
        else:
            self.useForMetaboliteGrouping.setCheckState(QtCore.Qt.Unchecked)
        if removeAsFalsePositive:
            self.removeAsFalsePositive.setCheckState(QtCore.Qt.Checked)
        else:
            self.removeAsFalsePositive.setCheckState(QtCore.Qt.Unchecked)

        self.groupName.setFocus()
        return self.exec_()

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        if event.mimeData().hasUrls:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()

            links = []
            for url in event.mimeData().urls():
                links.append(str(url.toLocalFile()).replace("\\", "/"))

            for link in links:
                if link.lower().endswith(".mzxml") or link.lower().endswith(".mzml"):
                    self.groupFiles.addItem(link)
                    self.groupfiles.append(link)

        else:
            event.ignore()


if __name__ == "__main__":
    import sys

    app = QtGui.QApplication(sys.argv)
    Dialog = groupEdit()

    Dialog.executeDialog()
    x = app.exec_()
    sys.exit(x)
