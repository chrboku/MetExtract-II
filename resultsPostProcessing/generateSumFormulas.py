import sys
sys.path.append("C:/PyMetExtract/PyMetExtract")

from SGR import SGRGenerator


import TableUtils

from copy import deepcopy



exID = "Num"
exMZ = "MZ"
exRT = "RT"
exAccMass = "M"
exXCount = "Xn"
exIonMode = "Ionisation_Mode"
exCharge = "Charge"


def processFile(file, columns, adducts, ppm=5., useAtoms=[], atomsRange=[], smCol="sumFormula_", toFile=None, useSevenGoldenRules=True, useCn=True):

    if toFile is None:
        toFile = file.replace(".tsv", ".SFs.tsv").replace(".txt", ".SFs.txt")

    if len(useAtoms) == 0:
        useAtoms = ["C", "H", "O", "N", "P"]
        atomsRange = [[-1, -1]]  #C
        atomsRange.append([0, 130])  #H
        atomsRange.append([0, 40])  #O
        atomsRange.append([0, 10])  #N
        atomsRange.append([0, 10])  #P

    table = TableUtils.TableUtils.readFile(file)

    if not smCol+"_CHO_count" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol+"_CHO_count", "INTEGER")
    if not smCol + "_CHON"+"_count" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHON"+"_count", "INTEGER")
    if not smCol + "_CHOP"+"_count" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHOP"+"_count", "INTEGER")
    if not smCol + "_CHOS"+"_count" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHOS"+"_count", "INTEGER")
    if not smCol + "_CHONP"+"_count" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHONP"+"_count", "INTEGER")
    if not smCol + "_CHONS"+"_count" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHONS"+"_count", "INTEGER")
    if not smCol + "_CHOPS"+"_count" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHOPS"+"_count", "INTEGER")
    if not smCol + "_CHONPS"+"_count" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHONPS"+"_count", "INTEGER")

    if not smCol+"_CHO" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol+"_CHO", "TEXT")
    if not smCol + "_CHON" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHON", "TEXT")
    if not smCol + "_CHOP" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHOP", "TEXT")
    if not smCol + "_CHOS" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHOS", "TEXT")
    if not smCol + "_CHONP" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHONP", "TEXT")
    if not smCol + "_CHONS" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHONS", "TEXT")
    if not smCol + "_CHOPS" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHOPS", "TEXT")
    if not smCol + "_CHONPS" in [col.name for col in table.getColumns()]:
        table.addColumn(smCol + "_CHONPS", "TEXT")

    if not "all_"+smCol+"_count" in [col.name for col in table.getColumns()]:
        table.addColumn("all_"+smCol+"_count", "INTEGER")

    class formulaSearch:
        def __init__(self, columns, adducts, ppm, useAtoms, atomsRange, useSevenGoldenRules, useCn):
            self.columns = columns
            self.adducts = adducts
            self.ppm = ppm
            self.useAtoms = useAtoms
            self.atomsRange = atomsRange
            self.useSevenGoldenRules=useSevenGoldenRules
            self.useCn=useCn

        def updateSumFormulaCol(self, x):

            sfg = SGRGenerator()
            mz = x[self.columns[exMZ]]
            if self.columns[exAccMass] in x.keys():
                accMass = x[self.columns[exAccMass]]
            else:
                accMass=""

            xCount = int(x[self.columns[exXCount]])
            atomsRange = deepcopy(self.atomsRange)
            charge = x[self.columns[exCharge]]
            ionMode = x[self.columns[exIonMode]]

            if self.useCn.lower()=="exact":
                atomsRange[0] = xCount
            elif self.useCn.lower()=="don't use":
                pass
            elif self.useCn.lower()=="min" or self.useCn.lower()=="minimum":
                atomsRange[0]=(xCount, atomsRange[0][1])
            elif self.useCn.lower().startswith("plusminus"):
                g=2
                if len("plusminus_")<len(self.useCn):
                    g=int(self.useCn[len("plusminus_"):])

                atomsRange[0]=(xCount-g, xCount+g)

            #print "--------"
            dbe = {}
            dbe[smCol + "_CHO"] = []
            dbe[smCol + "_CHOS"] = []
            dbe[smCol + "_CHOP"] = []
            dbe[smCol + "_CHON"] = []
            dbe[smCol + "_CHONP"] = []
            dbe[smCol + "_CHOPS"] = []
            dbe[smCol + "_CHONS"] = []
            dbe[smCol + "_CHONPS"] = []

            if accMass == "":
                for adduct in self.adducts:
                    if ionMode == adduct[2] and charge==adduct[3]:
                        ret = sfg.findFormulas((mz - adduct[1])*adduct[3]/adduct[4], self.ppm, useAtoms=useAtoms, atomsRange=atomsRange,
                                               fixed="C", useSevenGoldenRules=self.useSevenGoldenRules)
                        for e in ret:
                            #print e
                            ent = ["[M" + adduct[0] + "]: ", e]
                            if "P" in e and "N" in e and "S" in e:
                                dbe[smCol + "_CHONPS"].append(ent)
                            elif "P" in e and "S" in e:
                                dbe[smCol + "_CHOPS"].append(ent)
                            elif "S" in e and "N" in e:
                                dbe[smCol + "_CHONS"].append(ent)
                            elif "N" in e and "P" in e:
                                dbe[smCol + "_CHONP"].append(ent)
                            elif "P" in e:
                                dbe[smCol + "_CHOP"].append(ent)
                            elif "N" in e:
                                dbe[smCol + "_CHON"].append(ent)
                            elif "S" in e:
                                dbe[smCol + "_CHOS"].append(ent)
                            else:
                                dbe[smCol + "_CHO"].append(ent)
            else:
                for m in str(accMass).split(","):
                    m = m.replace("\"", "")
                    m = m.strip()
                    ret = sfg.findFormulas(float(m), self.ppm, useAtoms=useAtoms, atomsRange=atomsRange, fixed="C",
                                           useSevenGoldenRules=self.useSevenGoldenRules)
                    for e in ret:
                        #print e
                        ent = ["[M]: ", e]
                        if "P" in e and "N" in e and "S" in e:
                            dbe[smCol + "_CHONPS"].append(ent)
                        elif "P" in e and "S" in e:
                            dbe[smCol + "_CHOPS"].append(ent)
                        elif "S" in e and "N" in e:
                            dbe[smCol + "_CHONS"].append(ent)
                        elif "N" in e and "P" in e:
                            dbe[smCol + "_CHONP"].append(ent)
                        elif "P" in e:
                            dbe[smCol + "_CHOP"].append(ent)
                        elif "N" in e:
                            dbe[smCol + "_CHON"].append(ent)
                        elif "S" in e:
                            dbe[smCol + "_CHOS"].append(ent)
                        else:
                            dbe[smCol + "_CHO"].append(ent)

            if len(dbe[smCol + "_CHO"]) > 0:
                if len(dbe[smCol + "_CHO"])>250:
                    x[smCol + "_CHO"] = "Too many hits %d"%(len(dbe[smCol + "_CHO"]))
                else:
                    x[smCol + "_CHO"] = ", ".join([h[0] + str(h[1]) for h in dbe[smCol + "_CHO"]])
                x[smCol + "_CHO_count"] = len(dbe[smCol + "_CHO"])
            else:
                x[smCol + "_CHO"] = ""
                x[smCol + "_CHO_count"] = 0

            if len(dbe[smCol + "_CHOS"]) > 0:
                if len(dbe[smCol + "_CHOS"]) > 250:
                    x[smCol + "_CHOS"] = "Too many hits %d" % (len(dbe[smCol + "_CHOS"]))
                else:
                    x[smCol + "_CHOS"] = ", ".join([h[0] + str(h[1]) for h in dbe[smCol + "_CHOS"]])
                x[smCol + "_CHOS"+"_count"] = len(dbe[smCol + "_CHOS"])
            else:
                x[smCol + "_CHOS"] = ""
                x[smCol + "_CHOS"+"_count"] = 0

            if len(dbe[smCol + "_CHOP"]) > 0:
                if len(dbe[smCol + "_CHOP"]) > 250:
                    x[smCol + "_CHOP"] = "Too many hits %d" % (len(dbe[smCol + "_CHOP"]))
                else:
                    x[smCol + "_CHOP"] = ", ".join([h[0] + str(h[1]) for h in dbe[smCol + "_CHOP"]])
                x[smCol + "_CHOP"+"_count"] = len(dbe[smCol + "_CHOP"])
            else:
                x[smCol + "_CHOP"] = ""
                x[smCol + "_CHOP"+"_count"] = 0

            if len(dbe[smCol + "_CHON"]) > 0:
                if len(dbe[smCol + "_CHON"]) > 250:
                    x[smCol + "_CHON"] = "Too many hits %d" % (len(dbe[smCol + "_CHON"]))
                else:
                    x[smCol + "_CHON"] = ", ".join([h[0] + str(h[1]) for h in dbe[smCol + "_CHON"]])
                x[smCol + "_CHON"+"_count"] = len(dbe[smCol + "_CHON"])
            else:
                x[smCol + "_CHON"] = ""
                x[smCol + "_CHON"+"_count"] = 0

            if len(dbe[smCol + "_CHONP"]) > 0:
                if len(dbe[smCol + "_CHONP"]) > 250:
                    x[smCol + "_CHONP"] = "Too many hits %d" % (len(dbe[smCol + "_CHONP"]))
                else:
                    x[smCol + "_CHONP"] = ", ".join([h[0] + str(h[1]) for h in dbe[smCol + "_CHONP"]])
                x[smCol + "_CHONP"+"_count"] = len(dbe[smCol + "_CHONP"])
            else:
                x[smCol + "_CHONP"] = ""
                x[smCol + "_CHONP"+"_count"] = 0

            if len(dbe[smCol + "_CHONS"]) > 0:
                if len(dbe[smCol + "_CHONS"]) > 250:
                    x[smCol + "_CHONS"] = "Too many hits %d" % (len(dbe[smCol + "_CHONS"]))
                else:
                    x[smCol + "_CHONS"] = ", ".join([h[0] + str(h[1]) for h in dbe[smCol + "_CHONS"]])
                x[smCol + "_CHONS"+"_count"] = len(dbe[smCol + "_CHONS"])
            else:
                x[smCol + "_CHONS"] = ""
                x[smCol + "_CHONS"+"_count"] = 0

            if len(dbe[smCol + "_CHOPS"]) > 0:
                if len(dbe[smCol + "_CHOPS"]) > 250:
                    x[smCol + "_CHOPS"] = "Too many hits %d" % (len(dbe[smCol + "_CHOPS"]))
                else:
                    x[smCol + "_CHOPS"] = ", ".join([h[0] + str(h[1]) for h in dbe[smCol + "_CHOPS"]])
                x[smCol + "_CHOPS"+"_count"] = len(dbe[smCol + "_CHOPS"])
            else:
                x[smCol + "_CHOPS"] = ""
                x[smCol + "_CHOPS"+"_count"] = 0

            if len(dbe[smCol + "_CHONPS"]) > 0:
                if len(dbe[smCol + "_CHONPS"]) > 250:
                    x[smCol + "_CHONPS"] = "Too many hits %d" % (len(dbe[smCol + "_CHONPS"]))
                else:
                    x[smCol + "_CHONPS"] = ", ".join([h[0] + str(h[1]) for h in dbe[smCol + "_CHONPS"]])
                x[smCol + "_CHONPS"+"_count"] = len(dbe[smCol + "_CHONPS"])
            else:
                x[smCol + "_CHONPS"] = ""
                x[smCol + "_CHONPS"+"_count"] = 0

            x["all_"+smCol+"_count"] = len(dbe[smCol + "_CHO"])+len(dbe[smCol + "_CHOS"])+len(dbe[smCol + "_CHOP"])+len(dbe[smCol + "_CHON"])+len(dbe[smCol + "_CHONP"])+len(dbe[smCol + "_CHONS"])+ \
                                        len(dbe[smCol + "_CHOPS"])+ len(dbe[smCol + "_CHONPS"])

            return x

    x = formulaSearch(columns, adducts, ppm, useAtoms, atomsRange, useSevenGoldenRules=useSevenGoldenRules, useCn=useCn)
    table.applyFunction(x.updateSumFormulaCol, showProgress=True)

    table.addComment("## Sum formula generation. Adducts %s, ppm %.1f, atoms %s, atomsRange %s, seven golden rules %s, useXn %s"%(str(adducts), ppm, str(useAtoms), str(atomsRange), str(useSevenGoldenRules), str(useCn)))

    TableUtils.TableUtils.saveFile(table, toFile)


## adduct definition (name, m/z increment, polarity, charges, number of M
adductsP = [('+H', 1.007276, "+", 1, 1), ('+NH4', 18.033823, "+", 1, 1), ('+Na', 22.989218, "+", 1, 1), ('+K', 38.963158, "+", 1, 1)]
adductsN = [('-H', -1.007276, "-", 1, 1), ('+Na-2H', 20.974666, "-", 1, 1), ('+Cl', 34.969402, "-", 1, 1), ('+FA-H', 44.998201, "-", 1, 1), ('+Hac-H', 59.013851, '-', 1, 1), ('+Br', 78.918885, "-", 1, 1)]

import re

#taken from http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html
#use for a=["1", "2", "10", "11", "3"]
#natSort(a)
#for a=[("1", 1), ("2", 2), ("10", 3), ("11", 4), ("3", 5)]
#natSort(a, key=itemgetter(0))
def natSort(l, key=lambda ent: ent):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda ent, key=key: [convert(c) for c in re.split('([0-9]+)', str(key(ent)))]
    l.sort(key=alphanum_key)
    return l




def annotateResultsWithSumFormulas(resultsFile, useAtoms, atomsRange, Xn, useExactXn, ppm=5., adducts=adductsN+adductsP):

    processFile(resultsFile, toFile=resultsFile,
                columns={exID: exID, exMZ: exMZ, exRT: exRT, exAccMass: exAccMass, exXCount: exXCount,
                         exIonMode: exIonMode, exCharge: exCharge}, adducts=adducts, ppm=ppm,
                smCol="SFs",
                useAtoms=deepcopy(useAtoms), atomsRange=deepcopy(atomsRange), useCn=useExactXn,
                useSevenGoldenRules=True)
