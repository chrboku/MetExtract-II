########################################################################################################################
########################################################################################################################
#############################  Imports


from TableUtils import TableUtils
from utils import Bunch

import guidata
import guidata.dataset.datatypes as dt
import guidata.dataset.dataitems as di

import csv


########################################################################################################################
########################################################################################################################
#############################  Function definitions

## combine the results of two experiments using the m/z, rt, xn and loading properties of the detected feature pairs
def combineResults(resAFile, resBFile, resAName, resBName, saveToFile,
                   mergeName="", maxPPMError=5., maxRTShift=0.15, checkXn=False,
                   numCol="Num", mzCol="MZ", rtCol="RT", xnCol="Xn", chargeCol="Charge", ionModeCol="Ionisation_Mode",
                   importantCols=[]):
    ## read the results files
    resA = TableUtils.readFile(resAFile)
    resB = TableUtils.readFile(resBFile)

    ## fetch all feature pairs from the two results files
    resAFPs = {}
    resAFPslist = []
    for rowA in resA.getData(getResultsAsBunchObjects=True):
        num = getattr(rowA, numCol)
        resAFPs[num] = rowA
        resAFPslist.append(rowA)

    resBFPs = {}
    resBFPslist = []
    for rowB in resB.getData(getResultsAsBunchObjects=True):
        num = getattr(rowB, numCol)
        rowB.unused = True
        resBFPs[num] = rowB
        resBFPslist.append(rowB)

    ## perform outer join on mz, rt, xn, charge, and ionmode values
    matches = []
    for a in resAFPslist:
        anum = getattr(a, numCol)
        amz = getattr(a, mzCol)
        art = getattr(a, rtCol)
        axn = getattr(a, xnCol)
        az = getattr(a, chargeCol)
        aionMode = getattr(a, ionModeCol)

        wi = []
        for b in resBFPslist:
            bnum = getattr(b, numCol)
            bmz = getattr(b, mzCol)
            brt = getattr(b, rtCol)
            bxn = getattr(b, xnCol)
            bz = getattr(b, chargeCol)
            bionMode = getattr(b, ionModeCol)

            if az == bz and aionMode == bionMode and \
                                            abs(amz - bmz) * 1000000 / amz <= maxPPMError and abs(
                        art - brt) <= maxRTShift \
                    and (not checkXn or axn == bxn):
                wi.append(bnum)
                setattr(b, "unused", False)

        matches.append(Bunch(anum=anum, bnums=wi))

    for b in resBFPslist:
        if getattr(b, "unused"):
            matches.append(Bunch(anum=None, bnums=[getattr(b, numCol)]))

    ## fetch all columns of both result files
    firstCols = [numCol, mzCol, rtCol, xnCol, chargeCol, ionModeCol] + importantCols
    colsA = [x.getName() for x in resA.getColumns()]
    colsB = [x.getName() for x in resB.getColumns()]
    colsRes = ["Num"] + \
              [mzCol, rtCol, xnCol, chargeCol, ionModeCol] + \
              ["%s_%s" % (resAName, c) for c in firstCols] + \
              ["%s_%s" % (resBName, c) for c in firstCols] + \
              ["%s_%s" % (resAName, c) for c in colsA if c not in firstCols] + \
              ["%s_%s" % (resBName, c) for c in colsB if c not in firstCols]

    ## save results
    num = 1
    with open(saveToFile, "wb") as tsvFile:
        resWriter = csv.writer(tsvFile, delimiter="\t", quotechar="\"", quoting=csv.QUOTE_MINIMAL)

        resWriter.writerow(colsRes)

        for match in matches:
            row = [num]
            if match.anum != None:
                row.extend([getattr(resAFPs[match.anum], c) for c in [mzCol, rtCol, xnCol, chargeCol, ionModeCol]])
                row.extend([getattr(resAFPs[match.anum], c) for c in firstCols])
            elif len(match.bnums) == 1:
                row.extend([getattr(resBFPs[match.bnums[0]], c) for c in [mzCol, rtCol, xnCol, chargeCol, ionModeCol]])
                row.extend([""] * (len(firstCols)))
            else:
                row.extend(["Error"] + [""] * (len(firstCols) - 2))

            if len(match.bnums) == 0:
                row.extend([""] * len(firstCols))
            elif len(match.bnums) == 1:
                row.extend([getattr(resBFPs[match.bnums[0]], c) for c in firstCols])
            else:
                row.extend(["Several(%d)" % len(match.bnums)] + [""] * (len(firstCols) - 1))

            if match.anum != None:
                row.extend([getattr(resAFPs[match.anum], c) for c in colsA if c not in firstCols])
            else:
                row.extend(["" for c in colsA if c not in firstCols])

            if len(match.bnums) == 1:
                row.extend([getattr(resBFPs[match.bnums[0]], c) for c in colsB if c not in firstCols])
            else:
                row.extend(["" for c in colsA if c not in firstCols])

            resWriter.writerow(row)

            num = num + 1

        resWriter.writerow(["###### Processing %s" % resAName])
        for comment in resA.getComments():
            resWriter.writerow(["## %s" % comment])
        resWriter.writerow(["###### Processing %s" % resBName])
        for comment in resB.getComments():
            resWriter.writerow(["## Parameters %s" % comment])
        resWriter.writerow(["###### Processing combination of the two processing results"])

        b = Bunch()
        b.resAFile = resAFile
        b.resBFile = resBFile
        b.resAName = resAName
        b.resBName = resBName
        b.saveToFile = saveToFile
        b.mergeName = mergeName
        b.maxPPMError = maxPPMError
        b.maxRTShift = maxRTShift
        b.checkXn = checkXn
        resWriter.writerow(["## Combination parameters %s" % ((b.dumpAsJSon().replace("\"", "'")))])
        resWriter.writerow(
            ["## IMPORTANT: Nums, OGroups and other ID values cannot be compared accross different processings!"])


########################################################################################################################
########################################################################################################################
#############################  Execute script and show a graphical user interface to the user




## Combination for Benedikt. Comparison of different data evaluation parameters
if __name__ == "__main__":

    params = Bunch()
    params.expA = ""
    params.expADesc = ""
    params.expB = ""
    params.expBDesc = ""
    params.maxPPMError = 5
    params.maxRTError = 0.05
    params.checkXn = 1
    params.saveResTo = ""


    ## necessary, construct QT application
    _app = guidata.qapplication()  # not required if a QApplication has already been created


    ## annotate the input dialog
    class Processing(dt.DataSet):
        """Combine results of two experiments"""
        _bg1 = dt.BeginGroup("Input files")
        expA = di.FileOpenItem("Select results A", ("tsv", "txt"), params.expA)
        expADescString = di.StringItem("Experiment A prefix", params.expADesc)

        expB = di.FileOpenItem("Select results B", ("tsv", "txt"), params.expB)
        expBDescString = di.StringItem("Experiment B prefix", params.expBDesc)
        _eg1 = dt.EndGroup("Input files")

        _bg2 = dt.BeginGroup("Processing parameters")
        maxPPMError = di.FloatItem("Max. ppm error", min=0, max=1000, default=params.maxPPMError)
        maxRTError = di.FloatItem("Max. RT error [min]", min=0, max=10, default=params.maxRTError)
        checkXn = di.ChoiceItem("Require Xn identity", ["No", "Yes"], default=params.checkXn)
        _eg2 = dt.EndGroup("Processing parameters")

        saveResTo = di.FileSaveItem("Select new results file", ("tsv", "txt"), params.saveResTo)


    run = True

    while run:
        ## show the input dialog to the user
        dialog = Processing()
        if dialog.edit():

            params = Bunch()
            params.expA = str(dialog.expA)
            params.expADesc = str(dialog.expADescString)
            params.expB = str(dialog.expB)
            params.expBDesc = str(dialog.expBDescString)
            params.maxPPMError = float(dialog.maxPPMError)
            params.maxRTError = float(dialog.maxRTError)
            params.checkXn = [False, True][dialog.checkXn]
            params.saveResTo = str(dialog.saveResTo)

            ## combine the results from the two experiments
            print "Combining ... (%s and %s)" % (params.expADesc, params.expBDesc)
            combineResults(params.expA, params.expB, params.expADesc, params.expBDesc,
                           maxPPMError=params.maxPPMError, maxRTShift=params.maxRTError,
                           checkXn=params.checkXn, saveToFile=params.saveResTo,
                           importantCols=["OGroup", "Ion", "Loss", "M"])
        else:
            run = False
