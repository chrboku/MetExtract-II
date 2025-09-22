########################################################################################################################
########################################################################################################################

import sys

# sys.path.append("C:/PyMetExtract/PyMetExtract")  # Removed hardcoded path

from ..TableUtils import TableUtils
from ..utils import Bunch

from ..mePyGuis import combineResultsDialog
from PySide6 import QtCore, QtGui, QtWidgets

import csv
import os


########################################################################################################################
########################################################################################################################
#############################  Function definitions


## read table into table-like object


def readDataMatrixToTable(matFile):
    ## import data matrix
    res = Bunch(columns=[], data=[], comments=[])
    with open(matFile, "rb") as inF:
        for iline, line in enumerate(inF):
            line = line.strip().replace("\r", "").replace("\n", "")
            if iline == 0:
                res.columns = line.split("\t")
            else:
                if line.startswith("#"):
                    res.comments.append(line)
                else:
                    d = {}
                    cells = line.split("\t")
                    for i, h in enumerate(res.columns):
                        if i < len(cells):
                            d[h] = cells[i]
                        else:
                            d[h] = ""
                    res.data.append(Bunch(**d))

    ## check column types and parse them
    for col in res.columns:
        type = "int"
        for row in res.data:
            cell = getattr(row, col)
            if type == "int":
                try:
                    int(cell)
                except:
                    type = "float"
            if type == "float":
                try:
                    float(cell)
                except:
                    type = "str"
        for irow, row in enumerate(res.data):
            if type == "int":
                setattr(row, col, int(getattr(row, col)))
            elif type == "float":
                setattr(row, col, float(getattr(row, col)))

    return res


## combine the results of two experiments using the m/z, rt, xn and loading properties of the detected feature pairs
def combineResults(
    experiments,
    experimentsOrder,
    saveToFile,
    maxPPMError=5.0,
    maxRTShift=0.15,
    checkXn=False,
    numCol="Num",
    mzCol="MZ",
    rtCol="RT",
    xnCol="Xn",
    chargeCol="Charge",
    ionModeCol="Ionisation_Mode",
    importantCols=[],
):
    results = {}

    ## read the results files
    experimentsData = {}
    for expi, expDesc in enumerate(experimentsOrder):
        expFile = experiments[expDesc]
        res = readDataMatrixToTable(expFile)
        experimentsData[expDesc] = res

        results[expDesc] = {}
        for row in res.data:
            num = getattr(row, numCol)
            setattr(row, numCol, expDesc + str(num))
            ogroup = getattr(row, "OGroup")
            setattr(row, "OGroup", expDesc + str(ogroup))

            results[expDesc][num] = row

    ## outer join of all results
    with open(saveToFile, "w") as tsvFile:
        resWriter = csv.writer(tsvFile, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL)

        ## meta-info columns and important columns
        rowVals = [
            "Num",
            "Comment",
            "OGroup",
            "Results",
            "MeanMZ",
            "MeanRT",
            "IonisationMode",
            "Charge",
            "PPMDev",
            "RTDev",
        ]
        for expDesc in experimentsOrder:
            rowVals.extend(
                [
                    expDesc + "__Num",
                    expDesc + "__MZ",
                    expDesc + "__RT",
                    expDesc + "__Xn",
                    expDesc + "__Charge",
                    expDesc + "__Ionisation_Mode",
                ]
            )
            rowVals.extend([expDesc + "__" + t for t in importantCols])
        ## all remaining columns from the different experiments
        for expDesc in experimentsOrder:
            for colName in experimentsData[expDesc].columns:
                rowVals.append(expDesc + "__" + colName)
        resWriter.writerow(rowVals)

        resI = 1
        run = True
        while run:
            result = None
            for expDesc in experimentsOrder:
                if len(results[expDesc]) > 0:
                    result = next(iter(results[expDesc].values()))
                    break
            run = result != None

            ## match all results for the current one (variable result)
            minRT = 1000000
            maxRT = 0
            minMZ = 10000000
            maxMZ = 0
            if run:
                ## find possible matching feature pairs
                possibleHits = {}
                mzs = []
                rts = []
                ionMode = ""
                charge = 0
                for expDesc in experimentsOrder:
                    possibleHits[expDesc] = []
                    for pHitNum, pHitRow in results[expDesc].items():
                        if (
                            abs(getattr(pHitRow, mzCol) - getattr(result, mzCol)) * 1e6 / getattr(result, mzCol) <= maxPPMError
                            and abs(getattr(pHitRow, rtCol) - getattr(result, rtCol)) <= maxRTShift
                            and getattr(pHitRow, chargeCol) == getattr(result, chargeCol)
                            and getattr(pHitRow, ionModeCol) == getattr(result, ionModeCol)
                            and ((not checkXn) or getattr(pHitRow, xnCol) == getattr(result, xnCol))
                        ):
                            possibleHits[expDesc].append(pHitNum)
                            mzs.append(getattr(pHitRow, mzCol))
                            rts.append(getattr(pHitRow, rtCol))
                            ionMode = getattr(pHitRow, ionModeCol)
                            charge = getattr(pHitRow, chargeCol)

                            minRT = min(minRT, getattr(pHitRow, rtCol))
                            maxRT = max(maxRT, getattr(pHitRow, rtCol))
                            minMZ = min(minMZ, getattr(pHitRow, mzCol))
                            maxMZ = max(maxMZ, getattr(pHitRow, mzCol))

                ## write meta-information of matched feature pairs to new TSV file
                rowVals = [
                    resI,
                    "",
                    "",
                    sum([len(possibleHits[j]) for j in possibleHits.keys()]),
                    sum(mzs) / len(mzs),
                    sum(rts) / len(rts),
                    ionMode,
                    charge,
                    "%.5f" % ((maxMZ - minMZ) * 1e6 / getattr(result, mzCol)),
                    "%.2f" % (maxRT - minRT),
                ]
                for expDesc in experimentsOrder:
                    if len(possibleHits[expDesc]) == 0:
                        rowVals.extend(["", "", "", "", "", ""])
                        rowVals.extend(["" for t in importantCols])
                    elif len(possibleHits[expDesc]) == 1:
                        rowVals.extend(
                            [
                                getattr(results[expDesc][possibleHits[expDesc][0]], numCol),
                                getattr(results[expDesc][possibleHits[expDesc][0]], mzCol),
                                getattr(results[expDesc][possibleHits[expDesc][0]], rtCol),
                                getattr(results[expDesc][possibleHits[expDesc][0]], xnCol),
                                getattr(
                                    results[expDesc][possibleHits[expDesc][0]],
                                    chargeCol,
                                ),
                                getattr(
                                    results[expDesc][possibleHits[expDesc][0]],
                                    ionModeCol,
                                ),
                            ]
                        )
                        for t in importantCols:
                            rowVals.append(getattr(results[expDesc][possibleHits[expDesc][0]], t))
                    else:
                        rowVals.extend(
                            [
                                ";".join([expDesc + str(t) for t in possibleHits[expDesc]]),
                                "%f" % (sum([getattr(results[expDesc][t], mzCol) for t in possibleHits[expDesc]]) / len(possibleHits[expDesc])),
                                ";".join(["%.2f" % getattr(results[expDesc][t], rtCol) for t in possibleHits[expDesc]]),
                                ";".join([str(getattr(results[expDesc][t], xnCol)) for t in possibleHits[expDesc]]),
                                getattr(
                                    results[expDesc][possibleHits[expDesc][0]],
                                    chargeCol,
                                ),
                                getattr(
                                    results[expDesc][possibleHits[expDesc][0]],
                                    ionModeCol,
                                ),
                            ]
                        )
                        rowVals.extend([";".join([str(getattr(results[expDesc][t], x)) for t in possibleHits[expDesc]]) for x in importantCols])
                ## write remaining columns of all experiments to the new TSV file
                for expDesc in experimentsOrder:
                    if len(possibleHits[expDesc]) == 0:
                        for colName in experimentsData[expDesc].columns:
                            rowVals.append("")
                    elif len(possibleHits[expDesc]) == 1:
                        for colName in experimentsData[expDesc].columns:
                            rowVals.append(getattr(results[expDesc][possibleHits[expDesc][0]], colName))
                    else:
                        for colName in experimentsData[expDesc].columns:
                            rowVals.append(
                                ";".join(
                                    [
                                        str(
                                            getattr(
                                                results[expDesc][possibleHits[expDesc][ci]],
                                                colName,
                                            )
                                        )
                                        for ci in range(len(possibleHits[expDesc]))
                                    ]
                                )
                            )
                            # rowVals.append("")

                resWriter.writerow(rowVals)

                ## remove all matched feature pairs from further data processing steps
                for expDesc in experimentsOrder:
                    for pH in possibleHits[expDesc]:
                        del results[expDesc][pH]

            resI += 1

        ## add processing information as comments to new file
        for expDesc in experimentsOrder:
            resWriter.writerow(["###### Experiment %s" % expDesc])
            for comment in experimentsData[expDesc].comments:
                resWriter.writerow(["## %s" % comment])

        resWriter.writerow(["###### Processing combination of the processing results"])
        b = Bunch()
        b.maxPPMError = maxPPMError
        b.maxRTShift = maxRTShift
        b.checkXn = checkXn
        for expi, expDesc in enumerate(experimentsOrder):
            expFile = experiments[expDesc]
            resWriter.writerow(["## Experiment %s (file: %s)" % (expDesc, expFile)])
        resWriter.writerow(["## Combination parameters %s" % (b.dumpAsJSon().replace('"', "'"))])
        resWriter.writerow(["## !!! IMPORTANT: Nums, OGroups and other ID values cannot be compared accross different processings!"])


class CombineDialog:
    def __init__(self, dialog):
        self.ui = combineResultsDialog.Ui_Dialog()
        self.ui.setupUi(dialog)

        self.ui.expALoad.clicked.connect(self.selectExpA)
        self.ui.expBLoad.clicked.connect(self.selectExpB)
        self.ui.expCLoad.clicked.connect(self.selectExpC)
        self.ui.expDLoad.clicked.connect(self.selectExpD)
        self.ui.expELoad.clicked.connect(self.selectExpE)
        self.ui.expFLoad.clicked.connect(self.selectExpF)

        self.ui.loadResults.clicked.connect(self.selectResults)

        self.ui.run.clicked.connect(self.start)

        self.lastOpenDir = "."

    def selectExpA(self):
        sel = str(
            QtWidgets.QFileDialog.getOpenFileName(
                caption="Select experiment A file",
                dir=self.lastOpenDir,
                filter="Data matrix (*.tsv);;All files(*.*)",
            )
        )
        if sel is not None:
            self.ui.expASave.setText(sel)
            self.lastOpenDir = os.path.dirname(sel)

    def selectExpB(self):
        sel = str(
            QtWidgets.QFileDialog.getOpenFileName(
                caption="Select experiment B file",
                dir=self.lastOpenDir,
                filter="Data matrix (*.tsv);;All files(*.*)",
            )
        )
        if sel is not None:
            self.ui.expBSave.setText(sel)
            self.lastOpenDir = os.path.dirname(sel)

    def selectExpC(self):
        sel = str(
            QtWidgets.QFileDialog.getOpenFileName(
                caption="Select experiment C file",
                dir=self.lastOpenDir,
                filter="Data matrix (*.tsv);;All files(*.*)",
            )
        )
        if sel is not None:
            self.ui.expCSave.setText(sel)
            self.lastOpenDir = os.path.dirname(sel)

    def selectExpD(self):
        sel = str(
            QtWidgets.QFileDialog.getOpenFileName(
                caption="Select experiment D file",
                dir=self.lastOpenDir,
                filter="Data matrix (*.tsv);;All files(*.*)",
            )
        )
        if sel is not None:
            self.ui.expDSave.setText(sel)
            self.lastOpenDir = os.path.dirname(sel)

    def selectExpE(self):
        sel = str(
            QtWidgets.QFileDialog.getOpenFileName(
                caption="Select experiment E file",
                dir=self.lastOpenDir,
                filter="Data matrix (*.tsv);;All files(*.*)",
            )
        )
        if sel is not None:
            self.ui.expESave.setText(sel)
            self.lastOpenDir = os.path.dirname(sel)

    def selectExpF(self):
        sel = str(
            QtWidgets.QFileDialog.getOpenFileName(
                caption="Select experiment F file",
                dir=self.lastOpenDir,
                filter="Data matrix (*.tsv);;All files(*.*)",
            )
        )
        if sel is not None:
            self.ui.expFSave.setText(sel)
            self.lastOpenDir = os.path.dirname(sel)

    def selectResults(self):
        sel = str(
            QtWidgets.QFileDialog.getSaveFileName(
                caption="Select results file",
                dir=self.lastOpenDir,
                filter="Data matrix (*.tsv);;All files(*.*)",
            )
        )
        if sel is not None:
            self.ui.saveResults.setText(sel)
            self.lastOpenDir = os.path.dirname(sel)

    def start(self):
        if str(self.ui.expASave.text()) == "" or str(self.ui.expAPrefix.text()) == "":
            QtWidgets.QMessageBox.warning(
                None,
                "MetExtract",
                "Error: You need to specify at least the two experiments A and B",
                QtWidgets.QMessageBox.Ok,
            )
            return
        if str(self.ui.expBSave.text()) == "" or str(self.ui.expBPrefix.text()) == "":
            QtWidgets.QMessageBox.warning(
                None,
                "MetExtract",
                "Error: You need to specify at least the two experiments A and B",
                QtWidgets.QMessageBox.Ok,
            )
            return

        if str(self.ui.expCSave.text()) == "" != str(self.ui.expCPrefix.text()) == "":
            QtWidgets.QMessageBox.warning(
                None,
                "MetExtract",
                "Error: You need to specify both the input file and the prefix of experiment C",
                QtWidgets.QMessageBox.Ok,
            )
            return
        if str(self.ui.expDSave.text()) == "" != str(self.ui.expDPrefix.text()) == "":
            QtWidgets.QMessageBox.warning(
                None,
                "MetExtract",
                "Error: You need to specify both the input file and the prefix of experiment D",
                QtWidgets.QMessageBox.Ok,
            )
            return
        if str(self.ui.expESave.text()) == "" != str(self.ui.expEPrefix.text()) == "":
            QtWidgets.QMessageBox.warning(
                None,
                "MetExtract",
                "Error: You need to specify both the input file and the prefix of experiment E",
                QtWidgets.QMessageBox.Ok,
            )
            return
        if str(self.ui.expFSave.text()) == "" != str(self.ui.expFPrefix.text()) == "":
            QtWidgets.QMessageBox.warning(
                None,
                "MetExtract",
                "Error: You need to specify both the input file and the prefix of experiment F",
                QtWidgets.QMessageBox.Ok,
            )
            return

        if str(self.ui.saveResults.text()) == "":
            QtWidgets.QMessageBox.warning(
                None,
                "MetExtract",
                "Error: You need to specify a file to save the results to",
                QtWidgets.QMessageBox.Ok,
            )
            return

        experiments = {}
        experimentsOrder = []

        ex = str(self.ui.expAPrefix.text())
        experiments[ex] = str(self.ui.expASave.text())
        experimentsOrder.append(ex)

        ex = str(self.ui.expBPrefix.text())
        experiments[ex] = str(self.ui.expBSave.text())
        experimentsOrder.append(ex)

        ex = str(self.ui.expCPrefix.text())
        if ex != "":
            experiments[ex] = str(self.ui.expCSave.text())
            experimentsOrder.append(ex)

        ex = str(self.ui.expDPrefix.text())
        if ex != "":
            experiments[ex] = str(self.ui.expDSave.text())
            experimentsOrder.append(ex)

        ex = str(self.ui.expEPrefix.text())
        if ex != "":
            experiments[ex] = str(self.ui.expESave.text())
            experimentsOrder.append(ex)

        ex = str(self.ui.expFPrefix.text())
        if ex != "":
            experiments[ex] = str(self.ui.expFSave.text())
            experimentsOrder.append(ex)

        print("Combining results of experiments...")
        combineResults(
            experiments,
            experimentsOrder,
            maxPPMError=self.ui.maxPPMDev.value(),
            maxRTShift=self.ui.maxRTDev.value(),
            checkXn=False,
            saveToFile=str(self.ui.saveResults.text()),
            importantCols=["OGroup", "Ion", "Loss", "M"],
        )

        QtWidgets.QMessageBox.information(None, "MetExtract", "Processing finished..")
        sys.exit(0)


app = QtWidgets.QApplication(sys.argv)
Dialog = QtWidgets.QDialog()
ui = CombineDialog(Dialog)
Dialog.show()
sys.exit(app.exec())
