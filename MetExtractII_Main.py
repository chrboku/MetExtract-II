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
MetExtractVersion = "v2.12.1"


from mePyGuis.ModuleSelectionWindow import Ui_MainWindow
from PySide6 import QtGui, QtCore, QtWidgets
import subprocess
import platform

from mePyGuis import calcIsoEnrichmentDialog



class ModuleSelection(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, initDir=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.setWindowTitle("MetExtract II: Module selection")
        self.version.setText("MetExtract II (%s)"%MetExtractVersion)

        self.allExtractIcon.clicked.connect(self.openAllExtract)
        self.tracExtractIcon.clicked.connect(self.openTracExtract)
        self.fragExtractIcon.clicked.connect(self.openFragExtract)
        self.combineResultsButton.clicked.connect(self.openCombineResults)
        self.fticrExtractIcon.clicked.connect(self.openFTICR)
        self.documentationIcon.clicked.connect(self.openDocumentation)

        self.actionCalculate_isotopic_enrichment.triggered.connect(self.openCalcIsotopologEnrichment)


    def openMetExtractModule(self, module):
        try:
            import sys
            import os
            
            # Get the Python executable from the current environment
            python_exe = sys.executable
            current_dir = os.path.dirname(os.path.abspath(__file__))
            
            if module == "AllExtract" or module == "TracExtract":
                script_path = os.path.join(current_dir, "MExtract.py")
                cmd = [python_exe, script_path, "-m", module]
                subprocess.Popen(cmd)
            elif module == "FragExtract":
                script_path = os.path.join(current_dir, "FragExtract.py")
                cmd = [python_exe, script_path]
                subprocess.Popen(cmd)
            elif module == "combineResults":
                script_path = os.path.join(current_dir, "resultsPostProcessing", "combineResults.py")
                if os.path.exists(script_path):
                    cmd = [python_exe, script_path]
                    subprocess.Popen(cmd)
                else:
                    QtWidgets.QMessageBox.warning(self, "MetExtract", "combineResults module not found", QtWidgets.QMessageBox.Ok)
            elif module == "FTICRExtract":
                script_path = os.path.join(current_dir, "FTICRModule.py")
                cmd = [python_exe, script_path]
                subprocess.Popen(cmd)
            else:
                QtWidgets.QMessageBox.warning(self, "MetExtract", "Unknown module", QtWidgets.QMessageBox.Ok)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, "MetExtract", f"Requested module cannot be found: {str(e)}\nPlease try re-installing the software", QtWidgets.QMessageBox.Ok)


    def openCalcIsotopologEnrichment(self):
        diag = calcIsoEnrichmentDialog.calcIsoEnrichmentDialog()
        diag.executeDialog()

    def openAllExtract(self):
        print("starting AllExtract")
        self.openMetExtractModule(module="AllExtract")

    def openTracExtract(self):
        print("starting TracExtract")
        self.openMetExtractModule(module="TracExtract")

    def openFragExtract(self):
        print("starting FragExtract")
        self.openMetExtractModule(module="FragExtract")

    def openCombineResults(self):
        print("starting combineResults")
        self.openMetExtractModule(module="combineResults")

    def openFTICR(self):
        print("starting FTICRExtract")
        self.openMetExtractModule(module="FTICRExtract")

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

def main():
    import sys

    app = QtWidgets.QApplication(sys.argv)
    Dialog = ModuleSelection()

    Dialog.show()
    x = app.exec()

    sys.exit(x)

if __name__ == "__main__":
    main()
