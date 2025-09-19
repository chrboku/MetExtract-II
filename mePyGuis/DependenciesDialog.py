import sys
from PySide6 import QtGui, QtWidgets


def setLabelBackground(qlabel, colorName=None, r=255, g=0, b=0, alpha=255):
    qlabel.setAutoFillBackground(True)

    if colorName is not None:
        color = QtGui.QColor(colorName)
    else:
        color = QtGui.QColor(r, g, b)

    alpha = alpha
    values = "{r},{g},{b},{a}".format(
        r=color.red(), g=color.green(), b=color.blue(), a=alpha
    )
    qlabel.setStyleSheet(
        "QLabel { background-color: rgba("
        + values
        + "); border-width: 1px;border-style: solid;border-color: rgb(170, 170, 170);}"
    )


# a dialog that shows the results of several tasks (e.g. check for dependencies)
class DependenciesDialog:
    def __init__(
        self,
        dependencies={},
        statuss={0: "olivedrab", 1: "orange", 2: "firebrick"},
        dependencyOrder=None,
        parent=None,
    ):
        if dependencyOrder is None:
            dependencyOrder = dependencies.keys()
        self.pd = QtWidgets.QDialog(parent)
        self.pd.setModal(True)
        self.pd.setWindowTitle("Progress Wrapper")

        self.dependencies = {}
        self.dependenciesTexts = {}
        self.dependenciesIndicators = {}
        self.statuss = statuss

        l = QtWidgets.QGridLayout(self.pd)

        i = 0
        for dependency in dependencyOrder:
            assert dependency in dependencies.keys()
            status = dependencies[dependency]

            self.dependencies[dependency] = status

            if i > 0:
                line = QtWidgets.QFrame()
                line.setFrameShape(QtWidgets.QFrame.HLine)
                l.addWidget(line)
            i += 1

            k = QtWidgets.QWidget()
            l.addWidget(k)
            o = QtWidgets.QGridLayout(k)
            o.setContentsMargins(0, 0, 0, 0)

            text = QtWidgets.QLabel("")
            text.setFixedSize(15, 15)
            text.setToolTip("")
            setLabelBackground(text, colorName=statuss[status])
            self.dependenciesIndicators[dependency] = text
            o.addWidget(text, 0, 0)

            text = QtWidgets.QLabel(dependency)
            self.dependenciesTexts[dependency] = text
            o.addWidget(text, 0, 1)

    def show(self):
        self.pd.show()

    def hide(self):
        self.pd.hide()
        self.pd.close()

    def setDependencyText(self, dependency, newText):
        self.dependenciesTexts[dependency].setText(newText)

    def setDependencyStatus(self, dependency, newStatus):
        setLabelBackground(
            self.dependenciesIndicators[dependency], colorName=self.statuss[newStatus]
        )


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    # create a DependenciesDialog with 3 tasks and different statuses
    pw = DependenciesDialog(
        dependencies={"first": 1, "second": 2, "third": 0},
        dependencyOrder=["first", "second", "third"],
    )
    pw.show()

    sys.exit(app.exec())
