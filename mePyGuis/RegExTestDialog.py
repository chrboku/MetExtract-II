import sys
from PyQt4 import QtGui

from utils import natSort

import re


# a dialog that shows the results of several tasks (e.g. check for dependencies)
class RegExTestDialog(QtGui.QDialog):
    def __init__(self, strings=[], parent=None):
        super(RegExTestDialog, self).__init__()
        self.setModal(True)
        self.setWindowTitle("RegEx test dialog")
        self.strings=strings


        l = QtGui.QGridLayout(self)

        k = QtGui.QWidget()
        l.addWidget(k)
        o = QtGui.QGridLayout(k)
        o.setContentsMargins(0, 0, 0, 0)

        text = QtGui.QLabel("RegEx: ")
        #text.setFixedSize(15, 15)
        text.setToolTip("")
        o.addWidget(text, 0, 0)

        self.regEx = QtGui.QLineEdit(".*/(.*)_(.*).mzXML")
        self.regEx.textChanged.connect(self.updateReg)
        o.addWidget(self.regEx, 0, 1)

        text = QtGui.QLabel("RegEx: ")
        #text.setFixedSize(15, 15)
        text.setToolTip("")
        o.addWidget(text, 1, 0)

        self.regExRes = QtGui.QLineEdit("{0}")
        self.regExRes.textChanged.connect(self.updateReg)
        o.addWidget(self.regExRes, 1, 1)


        self.res=QtGui.QLabel("Parsed: \n\nEnter Regex")
        o.addWidget(self.res, 2, 1)

        self.okBut=QtGui.QPushButton("OK")
        self.okBut.clicked.connect(self.hide)
        o.addWidget(self.okBut, 3,1)


        self.updateReg()


    def updateReg(self):
        regEx=str(self.regEx.text())
        regExRes=str(self.regExRes.text())

        groups={}

        try:
            resA=[]
            for string in self.strings:
                try:
                    res=re.match(regEx, string)
                    finSt=regExRes.format(*res.groups())

                    if finSt not in groups.keys():
                        groups[finSt]=[]
                    groups[finSt].append(string)

                    resA.append("%s --> %s"%(string, finSt))
                except Exception as ex:
                    resA.append("%s --> Error in RegEx (%s)"%(string, ex.message))


            resA.append("")
            resA.append("Extracted groups:\n")
            for g in natSort(groups.keys()):
                resA.append("%s:"%(g))
                for f in groups[g]:
                    resA.append("   %s"%(f))
                resA.append("")


            self.res.setText("Parsed: \n\n%s"%("\n".join(resA)))

        except Exception as ex:
            self.res.setText("Error in RegEx/Result (%s)"%ex.message)

    def getRegEx(self):
        return str(self.regEx.text())
    def getRegExRes(self):
        return str(self.regExRes.text())


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    # create a DependenciesDialog with 3 tasks and different statuses
    pw = RegExTestDialog(strings=["TA_234234_S1_0h", "TA_34535_S2_0h", "TA_2342342_S1_4h"])
    pw.show()

    sys.exit(app.exec_())













