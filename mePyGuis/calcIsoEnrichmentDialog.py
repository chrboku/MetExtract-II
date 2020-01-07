from PyQt4 import QtGui, QtCore
from mePyGuis.calcIsotopeEnrichmentDialog import Ui_Dialog

from math import factorial

def choose(n, r):
    return factorial(n) // factorial(r) // factorial(n-r)


def calcIsoEnrichment(a,s,r):
    return choose(a,s)**(1./s) / (choose(a,s)**(1./s)+r**(1./s))

class calcIsoEnrichmentDialog(QtGui.QDialog, Ui_Dialog):
    def __init__(self, parent=None, initDir=None):
        QtGui.QDialog.__init__(self, parent)
        self.setupUi(self)
        self.setWindowTitle("Calculate isotopic enrichment")

        self.tableWidget.setColumnWidth(0, 120)
        self.tableWidget.setColumnWidth(1, 120)
        self.tableWidget.setColumnWidth(2, 120)
        self.tableWidget.setColumnWidth(3, 120)
        self.tableWidget.setColumnWidth(4, 60)
        self.tableWidget.setColumnWidth(5, 150)
        self.tableWidget.setColumnWidth(6, 150)
        self.tableWidget.cellChanged.connect(self.updateTable)

    def updateTable(self, x, y):
        if y<5:
            rowi=x
            try:
                m=float(self.tableWidget.item(rowi, 0).text())
                mpo=float(self.tableWidget.item(rowi, 1).text())
                cn=int(self.tableWidget.item(rowi, 4).text())

                self.tableWidget.setItem(rowi, 5, QtGui.QTableWidgetItem("%.4f%%"%(100.*calcIsoEnrichment(  cn, 1, mpo/m  ))))
            except Exception as ex:
                pass
            try:
                mpmo=float(self.tableWidget.item(rowi, 2).text())
                mp=float(self.tableWidget.item(rowi, 3).text())
                cn=int(self.tableWidget.item(rowi, 4).text())

                self.tableWidget.setItem(rowi, 6, QtGui.QTableWidgetItem("%.4f%%"%(100.*calcIsoEnrichment(  cn, 1, mpmo/mp  ))))
            except Exception as ex:
                pass


    def executeDialog(self):
        x = self.exec_()

        return x


if __name__ == "__main__":
    import sys

    app = QtGui.QApplication(sys.argv)
    Dialog = calcIsoEnrichmentDialog()

    Dialog.show()
    app.exec_()

    sys.exit()
