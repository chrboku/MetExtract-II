import sys
from PySide6 import QtGui, QtWidgets

from utils import natSort

import re


# a dialog that shows the results of several tasks (e.g. check for dependencies)
class RegExTestDialog(QtWidgets.QDialog):
    def __init__(self, strings=[], parent=None):
        super(RegExTestDialog, self).__init__()
        self.setModal(True)
        self.setWindowTitle("RegEx group import")
        self.strings = strings

        self.setMinimumHeight(700)
        self.setMinimumWidth(400)

        l = QtWidgets.QGridLayout(self)

        k = QtWidgets.QWidget()
        l.addWidget(k)
        o = QtWidgets.QGridLayout(k)
        o.setContentsMargins(0, 0, 0, 0)

        text = QtWidgets.QLabel("RegEx file path match: ")
        # text.setFixedSize(15, 15)
        text.setToolTip("")
        o.addWidget(text, 0, 0)

        self.regEx = QtWidgets.QLineEdit(".*/(.*?)_(.*).(mzXML|mzML)")
        self.regEx.textChanged.connect(self.updateReg)
        o.addWidget(self.regEx, 0, 1)

        text = QtWidgets.QLabel("RegEx group name: ")
        # text.setFixedSize(15, 15)
        text.setToolTip("")
        o.addWidget(text, 1, 0)

        self.regExRes = QtWidgets.QLineEdit("{0}")
        self.regExRes.textChanged.connect(self.updateReg)
        o.addWidget(self.regExRes, 1, 1)

        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidgetResizable(True)

        self.res = QtWidgets.QLabel("Parsed: \n\nEnter Regex")
        self.scroll.setWidget(self.res)
        o.addWidget(self.scroll, 2, 1)

        self.okBut = QtWidgets.QPushButton("OK")
        self.okBut.clicked.connect(self.hide)
        o.addWidget(self.okBut, 3, 1)

        self.updateReg()

    def updateReg(self):
        regEx = str(self.regEx.text())
        regExRes = str(self.regExRes.text())

        groups = {}

        try:
            resA = []
            for string in self.strings:
                try:
                    res = re.match(regEx, string)
                    finSt = regExRes.format(*res.groups())

                    if finSt not in groups.keys():
                        groups[finSt] = []
                    groups[finSt].append(string)

                    # resA.append("%s --> %s"%(string, finSt))
                except Exception as ex:
                    resA = ["Error in RegEx"]
                    # resA.append("%s --> Error in RegEx (%s)"%(string, ex.message))

            resA.append("\nExtracted groups:\n")
            for g in natSort(groups.keys()):
                resA.append("%s:" % (g))
                for f in groups[g]:
                    resA.append("   %s" % (f))
                resA.append("")

            self.res.setText("%s" % ("\n".join(resA)))

        except Exception as ex:
            self.res.setText("Error in RegEx/Result (%s)" % ex.message)

    def getRegEx(self):
        return str(self.regEx.text())

    def getRegExRes(self):
        return str(self.regExRes.text())


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    # create a DependenciesDialog with 3 tasks and different statuses
    pw = RegExTestDialog(
        strings=[
            "asldfkj/TA_234234_S1_0h.mzXML",
            "asldfkj/TA_34535_S2_0h.mzML",
            "asldfkj/TA_2342342_S1_4h.mzXML",
        ]
    )
    pw.show()

    sys.exit(app.exec())
