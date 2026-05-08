"""Custom about dialog for MetExtract applications."""

import os
from PySide6 import QtWidgets, QtGui, QtCore


class AboutDialog(QtWidgets.QDialog):
    """Custom about dialog with logo, institute info, and disclaimer."""

    def __init__(self, parent=None, app_name="MetExtract", version="", institute_name=""):
        """Initialize the about dialog.

        Parameters
        ----------
        parent : QWidget, optional
            Parent widget
        app_name : str
            Application name
        version : str
            Application version
        institute_name : str
            Institute name for the about dialog
        """
        super().__init__(parent)
        self.setWindowTitle(app_name)
        self.setMinimumWidth(500)

        # Set white background
        self.setStyleSheet("background-color: white;")

        layout = QtWidgets.QVBoxLayout()
        layout.setContentsMargins(20, 16, 20, 16)
        layout.setSpacing(10)
        layout.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)

        # Logo at the top
        logo_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
            "resources",
            "logo.png",
        )

        if os.path.exists(logo_path):
            logo_label = QtWidgets.QLabel()
            pixmap = QtGui.QPixmap(logo_path)
            scaled_pixmap = pixmap.scaledToWidth(200, QtCore.Qt.SmoothTransformation)
            logo_label.setPixmap(scaled_pixmap)
            logo_label.setAlignment(QtCore.Qt.AlignCenter)
            logo_label.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
            logo_label.mousePressEvent = lambda e: self._open_boku_website()
            layout.addWidget(logo_label)

        # Application title and version
        title_label = QtWidgets.QLabel(f"<b>{app_name} {version}</b>")
        title_label.setAlignment(QtCore.Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 21px;")
        layout.addWidget(title_label)

        # Institute info
        institute_label = QtWidgets.QLabel()
        institute_label.setOpenExternalLinks(True)
        institute_label.setAlignment(QtCore.Qt.AlignCenter)

        if institute_name == "iBAM":
            institute_text = '<div style="font-size: 18px;"><p>(c) Institute for Bioanalytics and Agrometabolomics (iBAM), IFA-Tulln<br/>University of Natural Resources and Life Sciences, Vienna</p><p><a href="https://boku.ac.at/en/agri/ibam/research-groups/plant-microbe-metabolomics">Plant Microbe Metabolomics Research Group</a></p></div>'
        elif institute_name == "CAC":
            institute_text = '<div style="font-size: 18px;"><p>(c) Centre for Analytical Chemistry, IFA Tulln<br/>University of Natural Resources and Life Sciences, Vienna</p><p><a href="https://boku.ac.at/en/agri/ibam/research-groups/plant-microbe-metabolomics">Plant Microbe Metabolomics Research Group</a></p></div>'
        else:
            institute_text = '<div style="font-size: 18px;"><p>(c) Institute for Bioanalytics and Agrometabolomics (iBAM), IFA-Tulln<br/>University of Natural Resources and Life Sciences, Vienna</p><p><a href="https://boku.ac.at/en/agri/ibam/research-groups/plant-microbe-metabolomics">Plant Microbe Metabolomics Research Group</a></p></div>'

        institute_label.setText(institute_text)
        layout.addWidget(institute_label)

        # vertical space
        layout.addSpacing(40)

        # Disclaimer
        disclaimer_text = (
            'THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, '
            "INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. "
            "IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, "
            "WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
        )
        disclaimer_label = QtWidgets.QLabel(disclaimer_text)
        disclaimer_label.setAlignment(QtCore.Qt.AlignCenter)
        disclaimer_label.setWordWrap(True)
        disclaimer_label.setStyleSheet("font-size: 11px; color: #666;")
        layout.addWidget(disclaimer_label)

        # OK button
        button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok)
        button_box.setStyleSheet("font-size: 14px;")
        button_box.accepted.connect(self.accept)
        layout.addWidget(button_box)

        self.setLayout(layout)

    def _open_boku_website(self):
        """Open the BOKU website when logo is clicked."""
        import webbrowser

        webbrowser.open_new_tab("https://boku.ac.at")
