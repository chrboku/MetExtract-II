########################################################################################################################
########################################################################################################################
#############################  Imports

import sys
sys.path.append("C:/PyMetExtract/PyMetExtract")

from TableUtils import TableUtils
from utils import Bunch

import guidata
import guidata.dataset.datatypes as dt
import guidata.dataset.dataitems as di

import csv


########################################################################################################################
########################################################################################################################
#############################  Function definitions


## read table into table-like object


def readDataMatrixToTable(matFile):

    ## import data matrix
    res = Bunch(columns=[], data=[], comments=[])
    with open(matFile, "rb") as inF:
        for iline, line in enumerate(inF):
            line=line.strip()
            if iline == 0:
                res.columns = line.split("\t")
            else:
                if line.startswith("#"):
                    res.comments.append(line)
                else:
                    d = {}
                    cells = line.split("\t")
                    for i, h in enumerate(res.columns):
                        d[h] = cells[i]
                    res.data.append(Bunch(**d))

    ## check column types and parse them
    for col in res.columns:
        type="int"
        for row in res.data:
            cell=getattr(row, col)
            if type=="int":
                try:
                   int(cell)
                except:
                    type="float"
            if type=="float":
                try:
                    float(cell)
                except:
                    type="str"
        for irow, row in enumerate(res.data):
            if type=="int":
                setattr(row, col, int(getattr(row, col)))
            elif type=="float":
                setattr(row, col, float(getattr(row, col)))

    return res


## combine the results of two experiments using the m/z, rt, xn and loading properties of the detected feature pairs
def combineResults(experiments, experimentsOrder, saveToFile,
                   maxPPMError=5., maxRTShift=0.15, checkXn=False,
                   numCol="Num", mzCol="MZ", rtCol="RT", xnCol="Xn", chargeCol="Charge", ionModeCol="Ionisation_Mode",
                   importantCols=[]):

    results = {}

    ## read the results files
    experimentsData={}
    for expi, expDesc in enumerate(experimentsOrder):
        expFile=experiments[expDesc]
        res = readDataMatrixToTable(expFile)
        experimentsData[expDesc] = res

        results[expDesc]={}
        for row in res.data:
            num = getattr(row, numCol)
            setattr(row, numCol, str(int(int(num)+(expi+1)*1E6)))
            ogroup = getattr(row, "OGroup")
            setattr(row, "OGroup", str(int(int(ogroup)+(expi+1)*1E6)))

            results[expDesc][num] = row


    ## outer join of all results
    with open(saveToFile, "wb") as tsvFile:
        resWriter = csv.writer(tsvFile, delimiter="\t", quotechar="\"", quoting=csv.QUOTE_MINIMAL)

        ## meta-info columns and important columns
        rowVals=["Num", "Comment", "OGroup", "Results", "MeanMZ", "MeanRT", "IonisationMode", "Charge", "PPMDev", "RTDev"]
        for expDesc in experimentsOrder:
            rowVals.extend([expDesc+"__Num", expDesc+"__MZ", expDesc+"__RT", expDesc+"__Xn", expDesc+"__Charge", expDesc+"__Ionisation_Mode"])
            rowVals.extend([expDesc+"__"+t for t in importantCols])
        ## all remaining columns from the different experiments
        for expDesc in experimentsOrder:
            for colName in experimentsData[expDesc].columns:
                rowVals.append(expDesc+"__"+colName)
        resWriter.writerow(rowVals)

        resI=1
        run = True
        while run:

            result=None
            for expDesc in experimentsOrder:
                if len(results[expDesc])>0:
                    result = results[expDesc].itervalues().next()
                    break
            run = result!=None

            ## match all results for the current one (variable result)
            minRT=1000000
            maxRT=0
            minMZ=10000000
            maxMZ=0
            if run:
                ## find possible matching feature pairs
                possibleHits={}
                mzs=[]
                rts=[]
                ionMode=""
                charge=0
                for expDesc in experimentsOrder:
                    possibleHits[expDesc]=[]
                    for pHitNum, pHitRow in results[expDesc].items():
                        if abs(getattr(pHitRow, mzCol)-getattr(result, mzCol))*1E6/getattr(result, mzCol) <= maxPPMError and \
                           abs(getattr(pHitRow, rtCol)-getattr(result, rtCol)) <= maxRTShift and \
                           getattr(pHitRow, chargeCol) == getattr(result, chargeCol) and \
                           getattr(pHitRow, ionModeCol) == getattr(result, ionModeCol) and \
                           ((not checkXn) or getattr(pHitRow, xnCol) == getattr(result, xnCol)):

                            possibleHits[expDesc].append(pHitNum)
                            mzs.append(getattr(pHitRow, mzCol))
                            rts.append(getattr(pHitRow, rtCol))
                            ionMode=getattr(pHitRow, ionModeCol)
                            charge=getattr(pHitRow, chargeCol)

                            minRT=min(minRT, getattr(pHitRow, rtCol))
                            maxRT=max(maxRT, getattr(pHitRow, rtCol))
                            minMZ=min(minMZ, getattr(pHitRow, mzCol))
                            maxMZ=max(maxMZ, getattr(pHitRow, mzCol))

                ## write meta-information of matched feature pairs to new TSV file
                rowVals = [resI, "", "", sum([len(possibleHits[j]) for j in possibleHits.keys()]), sum(mzs)/len(mzs), sum(rts)/len(rts), ionMode, charge, "%.5f"%((maxMZ-minMZ)*1E6/getattr(result, mzCol)), "%.2f"%(maxRT-minRT)]
                for expDesc in experimentsOrder:

                    if len(possibleHits[expDesc])==0:
                        rowVals.extend(["", "", "", "", "", ""])
                        rowVals.extend(["" for t in importantCols])
                    elif len(possibleHits[expDesc])==1:
                        rowVals.extend([getattr(results[expDesc][possibleHits[expDesc][0]], numCol),
                                        getattr(results[expDesc][possibleHits[expDesc][0]], mzCol),
                                        getattr(results[expDesc][possibleHits[expDesc][0]], rtCol),
                                        getattr(results[expDesc][possibleHits[expDesc][0]], xnCol),
                                        getattr(results[expDesc][possibleHits[expDesc][0]], chargeCol),
                                        getattr(results[expDesc][possibleHits[expDesc][0]], ionModeCol)])
                        for t in importantCols:
                            rowVals.append(getattr(results[expDesc][possibleHits[expDesc][0]], t))
                    else:
                        rowVals.extend([";".join([str(t) for t in possibleHits[expDesc]]),
                                        "%f"%(sum([getattr(results[expDesc][t], mzCol) for t in possibleHits[expDesc]])/len(possibleHits[expDesc])),
                                        ";".join(["%.2f"%getattr(results[expDesc][t], rtCol) for t in possibleHits[expDesc]]),
                                        ";".join([str(getattr(results[expDesc][t], xnCol)) for t in possibleHits[expDesc]]),
                                        getattr(results[expDesc][possibleHits[expDesc][0]], chargeCol),
                                        getattr(results[expDesc][possibleHits[expDesc][0]], ionModeCol)])
                        rowVals.extend([";".join([str(getattr(results[expDesc][t], x)) for t in possibleHits[expDesc]]) for x in importantCols])
                ## write remaining columns of all experiments to the new TSV file
                for expDesc in experimentsOrder:
                    if len(possibleHits[expDesc])==0:
                        for colName in experimentsData[expDesc].columns:
                            rowVals.append("")
                    elif len(possibleHits[expDesc])==1:
                        for colName in experimentsData[expDesc].columns:
                            rowVals.append(getattr(results[expDesc][possibleHits[expDesc][0]], colName))
                    else:
                        for colName in experimentsData[expDesc].columns:
                            rowVals.append("")

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
            resWriter.writerow(["## Experiment %s (file: %s)"%(expDesc, expFile)])
        resWriter.writerow(["## Combination parameters %s" % ((b.dumpAsJSon().replace("\"", "'")))])
        resWriter.writerow(["## !!! IMPORTANT: Nums, OGroups and other ID values cannot be compared accross different processings!"])









########################################################################################################################
########################################################################################################################
#############################  Execute script and show a graphical user interface to the user
if __name__ == "__main__":

    params = Bunch()
    params.expA = "D:/Maria_20171219/EVAL_171206_FML_Remus_untreated_FH/EVAL_171206_FML_Remus_untreated_FH/test.tsv"
    params.expADesc = "FML"
    params.expB = "D:/Maria_20171219/EVAL_171206_PHE_Remus_untreated_FH/EVAL_171206_PHE_Remus_untreated_FH/test.tsv"
    params.expBDesc = "PHE"
    params.expC = "D:/Maria_20171219/EVAL_171206_TRP_Remus_untreated_FH/test.tsv"
    params.expCDesc = "TRP"
    params.expD = ""
    params.expDDesc = ""
    params.expE = ""
    params.expEDesc = ""
    params.expF = ""
    params.expFDesc = ""

    params.maxPPMError = 5
    params.maxRTError = 0.15
    params.checkXn = 1
    params.saveResTo = "D:/Maria_20171219/results_FML_PHE_TRP.tsv"


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

        expC = di.FileOpenItem("Select results C", ("tsv", "txt"), params.expC, check=False)
        expCDescString = di.StringItem("Experiment C prefix", params.expCDesc)
        _eg1 = dt.EndGroup("Input files")

        expD = di.FileOpenItem("Select results D", ("tsv", "txt"), params.expD, check=False)
        expDDescString = di.StringItem("Experiment D prefix", params.expDDesc)
        _eg1 = dt.EndGroup("Input files")

        expE = di.FileOpenItem("Select results E", ("tsv", "txt"), params.expE, check=False)
        expEDescString = di.StringItem("Experiment E prefix", params.expEDesc)
        _eg1 = dt.EndGroup("Input files")

        expF = di.FileOpenItem("Select results F", ("tsv", "txt"), params.expF, check=False)
        expFDescString = di.StringItem("Experiment F prefix", params.expFDesc)
        _eg1 = dt.EndGroup("Input files")

        _bg2 = dt.BeginGroup("Processing parameters")
        maxPPMError = di.FloatItem("Max. ppm error", min=0, max=1000, default=params.maxPPMError)
        maxRTError = di.FloatItem("Max. RT error [min]", min=0, max=10, default=params.maxRTError)
        _eg2 = dt.EndGroup("Processing parameters")

        saveResTo = di.FileSaveItem("Select new results file", ("tsv", "txt"), params.saveResTo)


    ## show the input dialog to the user
    dialog = Processing()
    if dialog.edit():

        params = Bunch()
        params.expA = str(dialog.expA)
        params.expADesc = str(dialog.expADescString)
        params.expB = str(dialog.expB)
        params.expBDesc = str(dialog.expBDescString)
        params.expC = str(dialog.expC)
        params.expCDesc = str(dialog.expCDescString)
        params.expD = str(dialog.expD)
        params.expDDesc = str(dialog.expDDescString)
        params.expE = str(dialog.expE)
        params.expEDesc = str(dialog.expEDescString)
        params.expF = str(dialog.expF)
        params.expFDesc = str(dialog.expFDescString)

        params.maxPPMError = float(dialog.maxPPMError)
        params.maxRTError = float(dialog.maxRTError)
        params.saveResTo = str(dialog.saveResTo)

        experiments={}
        experimentsOrder=[]
        if params.expADesc != "":
            experimentsOrder.append(params.expADesc)
            experiments[params.expADesc] = params.expA
        if params.expBDesc != "":
            experimentsOrder.append(params.expBDesc)
            experiments[params.expBDesc] = params.expB
        if params.expCDesc != "":
            experimentsOrder.append(params.expCDesc)
            experiments[params.expCDesc] = params.expC
        if params.expDDesc != "":
            experimentsOrder.append(params.expDDesc)
            experiments[params.expDDesc] = params.expD
        if params.expEDesc != "":
            experimentsOrder.append(params.expEDesc)
            experiments[params.expEDesc] = params.expE
        if params.expFDesc != "":
            experimentsOrder.append(params.expFDesc)
            experiments[params.expFDesc] = params.expF

        ## combine the results from the two experiments
        print "Combining results of %d experiments..." % (len(experiments))
        combineResults(experiments, experimentsOrder,
                       maxPPMError=params.maxPPMError, maxRTShift=params.maxRTError,
                       checkXn=False, saveToFile=params.saveResTo,
                       importantCols=["OGroup", "Ion", "Loss", "M"])

