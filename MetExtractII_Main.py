#    MetExtract II
#    Copyright (C) 2015
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

# current MetExtract version
MetExtractVersion = "v2.9.0"


from mePyGuis.ModuleSelectionWindow import Ui_MainWindow
from PyQt4 import QtGui, QtCore
import subprocess
import platform

from mePyGuis import calcIsoEnrichmentDialog



class ModuleSelection(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, initDir=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.setWindowTitle("MetExtract II: Module selection")
        self.version.setText("MetExtract II (%s)"%MetExtractVersion)

        self.allExtractIcon.clicked.connect(self.openAllExtract)
        self.tracExtractIcon.clicked.connect(self.openTracExtract)
        self.fragExtractIcon.clicked.connect(self.openFragExtract)
        self.combineResultsButton.clicked.connect(self.openCombineResults)
        self.documentationIcon.clicked.connect(self.openDocumentation)

        self.actionCalculate_isotopic_enrichment.triggered.connect(self.openCalcIsotopologEnrichment)


    def openMetExtractModule(self, module):
        try:
            if module is "AllExtract" or module is "TracExtract":
                modStart="MExtract"
                if platform.system() == "Darwin":  #MAC
                    subprocess.Popen("%s -m %s"%(modStart, module))
                if platform.system() == "Windows":  #Windows
                    subprocess.Popen("%s.exe -m %s"%(modStart, module))
                if platform.system() == "Linux":  #Linux
                    subprocess.Popen("%s -m %s"%(modStart, module))
            elif module is "FragExtract":
                modStart="FragExtract"
                if platform.system() == "Darwin":  #MAC
                    subprocess.Popen("%s"%(modStart))
                if platform.system() == "Windows":  #Windows
                    subprocess.Popen("%s.exe"%(modStart))
                if platform.system() == "Linux":  #Linux
                    subprocess.Popen("%s"%(modStart))
            elif module is "combineResults":
                if platform.system() == "Darwin":  #MAC
                    subprocess.Popen("combineResults")
                if platform.system() == "Windows":  #Windows
                    subprocess.Popen("combineResults.exe")
                if platform.system() == "Linux":  #Linux
                    subprocess.Popen("combineResults")
            else:
                QtGui.QMessageBox.warning(self, "MetExtract", "Unknown module", QtGui.QMessageBox.Ok)
        except:
            QtGui.QMessageBox.warning(self, "MetExtract", "Requested module cannot be found. \nPlease try re-installing the software", QtGui.QMessageBox.Ok)


    def openCalcIsotopologEnrichment(self):
        diag = calcIsoEnrichmentDialog.calcIsoEnrichmentDialog()
        diag.executeDialog()

    def openAllExtract(self):
        print "starting AllExtract"
        self.openMetExtractModule(module="AllExtract")

    def openTracExtract(self):
        print "starting TracExtract"
        self.openMetExtractModule(module="TracExtract")

    def openFragExtract(self):
        print "starting FragExtract"
        self.openMetExtractModule(module="FragExtract")

    def openCombineResults(self):
        print "starting combineResults"
        self.openMetExtractModule(module="combineResults")

    def openDocumentation(self):
        import subprocess
        import webbrowser
        import sys
        from utils import get_main_dir

        url=get_main_dir()+"/documentation/index.html"
        if sys.platform == "darwin":    # in case of OS X
            subprocess.Popen(['open', url])
        else:
            webbrowser.open_new_tab(url)

if __name__ == "__main__":
    import sys

    app = QtGui.QApplication(sys.argv)
    Dialog = ModuleSelection()

    Dialog.show()
    x = app.exec_()

    sys.exit(x)
