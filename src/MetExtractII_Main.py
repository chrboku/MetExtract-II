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
MetExtractVersion = "v3.1.0"


def get_version_from_toml(file_path):
    import tomllib

    try:
        with open(file_path, "rb") as f:
            data = tomllib.load(f)
            return data.get("project", {}).get("version", "unknown")
    except Exception as e:
        print(f"Error reading version from TOML file: {e}")
        return "unknown"


# Update MetExtractVersion from the TOML file
MetExtractVersion = get_version_from_toml("pyproject.toml")

# Handle both relative imports (when used as module) and absolute imports (when run directly)
try:
    from .mePyGuis.ModuleSelectionWindow import Ui_MainWindow
    from .mePyGuis import calcIsoEnrichmentDialog
except ImportError:
    # When run directly, use absolute imports
    import sys
    import os

    # Handle case where __file__ might not be defined (e.g., when using exec())
    if "__file__" in globals():
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    else:
        sys.path.insert(0, os.path.dirname(os.path.abspath(os.getcwd() + "/src")))
    from mePyGuis.ModuleSelectionWindow import Ui_MainWindow
    from mePyGuis import calcIsoEnrichmentDialog

from PySide6 import QtGui, QtCore, QtWidgets
import subprocess
import platform


class ModuleSelection(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, initDir=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.setWindowTitle("MetExtract II: Module selection")
        self.version.setText("MetExtract II (%s)" % MetExtractVersion)

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
            args = []
            if module == "AllExtract" or module == "TracExtract":
                args = ["uv", "run", "python", "-m", "src.MExtract", "-m", module]
            elif module == "FragExtract":
                args = ["uv", "run", "python", "-m", "src.FragExtract"]
            elif module == "combineResults":
                args = ["uv", "run", "python", "-m", "src.combineResults"]
            elif module == "FTICRExtract":
                args = ["uv", "run", "python", "-m", "src.FTICRModule"]
            else:
                QtWidgets.QMessageBox.warning(self, "MetExtract", "Unknown module", QtWidgets.QMessageBox.Ok)
            subprocess.Popen(args)
        except Exception as e:
            QtWidgets.QMessageBox.warning(
                self,
                "MetExtract",
                f"Requested module cannot be found: {str(e)}\nPlease try re-installing the software",
                QtWidgets.QMessageBox.Ok,
            )

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
        from .utils import get_main_dir

        url = get_main_dir() + "/documentation/index.html"
        if sys.platform == "darwin":  # in case of OS X
            subprocess.Popen(["open", url])
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
