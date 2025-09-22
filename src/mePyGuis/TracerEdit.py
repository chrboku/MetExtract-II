import sys
import os
import pickle
import base64
from copy import copy, deepcopy

from PySide6 import QtCore, QtGui, QtWidgets

from .TracerEditor import Ui_Dialog
from ..utils import getRatio, getXCombinations
from ..formulaTools import formulaTools, getIsotopeMass


class ConfiguredTracer:
    def __init__(
        self,
        name="Deoxynivalenol",
        elementCount=15,
        isotopeA="12C",
        isotopeB="13C",
        enrichmentA=0.9893,
        enrichmentB=0.995,
        amountA=50.0,
        amountB=50.0,
        monoisotopicRatio=1,
        checkRatio=True,
        maxRelNegBias=70,
        maxRelPosBias=130,
        tracerType="user",
        id=-1,
        mzDelta=None,
    ):
        self.id = id

        self.name = name  # 0
        self.elementCount = elementCount  # 1
        self.isotopeA = isotopeA  # 2
        self.isotopeB = isotopeB  # 3
        self.enrichmentA = enrichmentA  # 4
        self.enrichmentB = enrichmentB  # 5
        self.amountA = amountA  # 6
        self.amountB = amountB  # 7
        self.monoisotopicRatio = monoisotopicRatio  # 8
        self.checkRatio = checkRatio  # 9
        self.maxRelNegBias = maxRelNegBias  # 10
        self.maxRelPosBias = maxRelPosBias  # 11
        self.tracerType = tracerType  # 12

        if mzDelta is None:
            self.mzDelta = getIsotopeMass(self.isotopeB)[0] - getIsotopeMass(self.isotopeA)[0]
        else:
            self.mzDelta = mzDelta

    def __str__(self):
        return "ConfiguredTracer: %s %s %s %d  enrichment: %.3f %.3f amount %.1f %.1f monoisotopicRatio %.3f (check ratio %s) bias: %.1f %.1f" % (
            self.name,
            self.isotopeA,
            self.isotopeB,
            self.elementCount,
            self.enrichmentA,
            self.enrichmentB,
            self.amountA,
            self.amountB,
            self.monoisotopicRatio,
            self.checkRatio,
            self.maxRelNegBias,
            self.maxRelPosBias,
        )


def getShortName(element):
    fT = formulaTools()
    for x in fT.elemDetails:
        d = fT.elemDetails[x]
        if d[0] == element:
            return x
    return ""


class tracerEdit(QtWidgets.QDialog, Ui_Dialog):
    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self, parent)
        self.setupUi(self)
        self.setWindowTitle("Tracer editor")

        self.configuredTracer = ConfiguredTracer()

        self.discardButton.clicked.connect(self.dialogCan)
        self.acceptButton.clicked.connect(self.dialogFin)

        self.acceptButton.setFocus()

        self.tName.setText(self.configuredTracer.name)
        self.tAtomCount.setValue(self.configuredTracer.elementCount)
        self.tNativeIsotope.setText(self.configuredTracer.isotopeA)
        self.tNativeEnrichment.setValue(self.configuredTracer.enrichmentA * 100.0)
        self.tLabeledIsotope.setText(self.configuredTracer.isotopeB)
        self.tLabeledEnrichment.setValue(self.configuredTracer.enrichmentB * 100.0)
        self.tAmountA.setValue(self.configuredTracer.amountA)
        self.tAmountB.setValue(self.configuredTracer.amountB)
        self.tCheckRatio.setChecked(True)
        self.tMinRatio.setValue(self.configuredTracer.maxRelNegBias)
        self.tMaxRatio.setValue(self.configuredTracer.maxRelPosBias)
        self.updateTracer()

        self.tAtomCount.valueChanged.connect(self.updateTracer)
        self.tNativeEnrichment.valueChanged.connect(self.updateTracer)
        self.tLabeledEnrichment.valueChanged.connect(self.updateTracer)
        self.tAmountA.valueChanged.connect(self.updateTracer)
        self.tAmountB.valueChanged.connect(self.updateTracer)
        self.tCheckRatio.stateChanged.connect(self.updateTracer)
        self.tMinRatio.valueChanged.connect(self.updateTracer)
        self.tMaxRatio.valueChanged.connect(self.updateTracer)

    def updateTracer(self):
        ratio = getRatio(self.tNativeEnrichment.value() / 100.0, self.tAtomCount.value(), 0) * self.tAmountA.value() / (getRatio(self.tLabeledEnrichment.value() / 100.0, self.tAtomCount.value(), 0) * self.tAmountB.value())
        self.tTheoreticalSignalRatio.setText("M:M' ~ %.3f" % ratio)

        self.tRatioGroupBox.setEnabled(self.tCheckRatio.isChecked())

        # self.tMinRatio.setPrefix("< %.4f * "%(ratio))
        # self.tMaxRatio.setPrefix("> %.4f * "%(ratio))
        self.tMaxRatio.setValue
        self.rRatioInfo.setText(
            "Ratio native monoisotopic / fully labeled (abundance / abundance): %.5f-%.5f"
            % (
                ratio * float(self.tMinRatio.value()) / 100.0,
                ratio * float(self.tMaxRatio.value()) / 100.0,
            )
        )

        return ratio

    def setTracer(self, tracer):
        if tracer is not None:
            self.configuredTracer = tracer
        else:
            self.configuredTracer = ConfiguredTracer()

        self.tName.setText(self.configuredTracer.name)
        self.tAtomCount.setValue(self.configuredTracer.elementCount)
        self.tNativeIsotope.setText(self.configuredTracer.isotopeA)
        self.tLabeledIsotope.setText(self.configuredTracer.isotopeB)
        self.tNativeEnrichment.setValue(self.configuredTracer.enrichmentA * 100.0)
        self.tLabeledEnrichment.setValue(self.configuredTracer.enrichmentB * 100.0)
        self.tAmountA.setValue(self.configuredTracer.amountA)
        self.tAmountB.setValue(self.configuredTracer.amountB)
        self.tCheckRatio.setChecked(self.configuredTracer.checkRatio)
        self.tMinRatio.setValue(self.configuredTracer.maxRelNegBias * 100.0)
        self.tMaxRatio.setValue(self.configuredTracer.maxRelPosBias * 100.0)
        self.mzDelta = getIsotopeMass(self.configuredTracer.isotopeB)[0] - getIsotopeMass(self.configuredTracer.isotopeA)[0]
        self.updateTracer()

    def getTracer(self):
        x = ConfiguredTracer(
            name=str(self.tName.text()),
            elementCount=int(self.tAtomCount.value()),
            isotopeA=str(self.tNativeIsotope.text()),
            isotopeB=str(self.tLabeledIsotope.text()),
            enrichmentA=float(self.tNativeEnrichment.value()) / 100.0,
            enrichmentB=float(self.tLabeledEnrichment.value()) / 100.0,
            amountA=float(self.tAmountA.value()),
            amountB=float(self.tAmountB.value()),
            monoisotopicRatio=self.updateTracer(),
            checkRatio=self.tCheckRatio.isChecked(),
            maxRelNegBias=float(self.tMinRatio.value()) / 100.0,
            maxRelPosBias=float(self.tMaxRatio.value()) / 100.0,
            tracerType="user",
        )
        self.configuredTracer = x
        return x

    def checkConfiguration(self):
        ok = True
        severe = False
        tracer = self.getTracer()

        ma, ea = getIsotopeMass(tracer.isotopeA)
        mb, eb = getIsotopeMass(tracer.isotopeB)
        if tracer.name == "":
            ok = False
            severe = True
            QtWidgets.QMessageBox.question(
                self,
                "MetExtract",
                "Error in tracer number %d \nYou did not specify a name. \nPlease provide a name",
                QtWidgets.QMessageBox.Ok,
            )

        if len(ea) == 0 or ea != eb:
            ok = False
            QtWidgets.QMessageBox.question(
                self,
                "MetExtract",
                "Error in tracer %s \nYou cannot use two different elements for the labelling process. \nPlease enter isotopes of the same element" % tracer.name,
                QtWidgets.QMessageBox.Ok,
            )

        return ok and not severe

    def dialogCan(self):
        self.reject()

    def dialogFin(self):
        if self.checkConfiguration():
            self.accept()

    def executeDialog(self):
        x = self.exec()

        return x


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    Dialog = tracerEdit()

    Dialog.show()
    x = app.exec()

    print("final tracer:", Dialog.getTracer())

    sys.exit(x)
