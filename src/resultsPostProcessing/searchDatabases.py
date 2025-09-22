import sys

sys.path.append("C:/development/PyMetExtract")

from ..formulaTools import formulaTools
from ..utils import is_float

from .. import TableUtils

import csv

from copy import deepcopy
from math import floor, ceil


exID = "Num"
exMZ = "MZ"
exRT = "RT"
exAccMass = "M"
exXCount = "Xn"
exIonMode = "Ionisation_Mode"
exCharge = "Charge"


import logging
from .. import LoggingSetup

LoggingSetup.LoggingSetup.Instance().initLogging()


class DBEntry:
    def __init__(
        self,
        dbName,
        num,
        name,
        sumFormula,
        mass,
        rt_min,
        mz,
        polarity,
        additionalInfo,
        hitType="",
    ):
        self.dbName = dbName
        self.num = num
        self.name = name
        self.sumFormula = sumFormula
        self.mass = mass
        self.matchErrorPPM = -1
        self.matchErrorMass = -1
        self.rt_min = rt_min
        self.mz = mz
        self.polarity = polarity
        self.additionalInfo = additionalInfo

        self.hitType = hitType

    def __str__(self):
        return (
            "DB: "
            + str(self.dbName)
            + " Num: "
            + str(self.num)
            + " Name: "
            + str(self.name)
            + " SumFormula: "
            + str(self.sumFormula)
            + " Mass: "
            + str(self.mass)
            + " matchErrorPPM: "
            + str(self.matchErrorPPM)
            + " matchErrorMass: "
            + str(self.matchErrorMass)
            + " rtMin: "
            + str(self.rt_min)
            + " mz: "
            + str(self.mz)
            + " polarity: "
            + str(self.polarity)
            + " additionalInfo: "
            + str(self.additionalInfo)
        )


class DBSearch:
    def __init__(self):
        self.dbEntriesNeutral = []
        self.dbEntriesMZ = []

    def addEntriesFromFile(self, dbName, dbFile, callBackCheckFunction=None):
        imported = 0
        notImported = 0

        curEntriesCount = len(self.dbEntriesMZ) + len(self.dbEntriesNeutral)

        fT = formulaTools()

        with open(dbFile) as fIn:
            tsvin = csv.reader(fIn, delimiter="\t")

            headers = {}
            for rowi, row in enumerate(tsvin):
                if len(row) == 0:  # Skip empty lines
                    continue
                if row[0].startswith("#"):
                    continue
                if rowi == 0:
                    for j in range(len(row)):
                        headers[row[j]] = j

                else:
                    try:
                        num = row[headers["Num"]].strip().replace('"', "DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n", "RETURN").replace("\r", "CarrRETURN").replace("#", "HASH")
                        name = row[headers["Name"]].strip().replace('"', "DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n", "RETURN").replace("\r", "CarrRETURN").replace("#", "HASH")
                        sumFormula = row[headers["SumFormula"]].strip().replace('"', "DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n", "RETURN").replace("\r", "CarrRETURN").replace("#", "HASH")
                        rt_min = float(row[headers["Rt_min"]]) if row[headers["Rt_min"]] != "" else None
                        mz = float(row[headers["MZ"]]) if row[headers["MZ"]] != "" else None
                        polarity = row[headers["IonisationMode"]].strip().replace('"', "DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n", "RETURN").replace("\r", "CarrRETURN").replace("#", "HASH")
                        additionalInfo = {}
                        for header in headers.keys():
                            if header not in [
                                "Num",
                                "Name",
                                "SumFormula",
                                "Rt_min",
                                "MZ",
                                "IonisationMode",
                            ]:
                                additionalInfo[header] = row[headers[header]].replace('"', "DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n", "RETURN").replace("\r", "CarrRETURN").replace("#", "HASH")

                        mass = 0
                        if sumFormula != "":
                            try:
                                elems = fT.parseFormula(sumFormula)
                                mass = fT.calcMolWeight(elems)
                            except Exception as ex:
                                logging.error("DB import error (%s, row: %d): The sumformula (%s) of the entry %s '%s' could not be parsed" % (dbName, rowi, sumFormula, num, name))
                                notImported += 1

                        dbEntry = DBEntry(
                            dbName,
                            num,
                            name,
                            sumFormula,
                            mass,
                            rt_min,
                            mz,
                            polarity,
                            additionalInfo,
                        )

                        use = True
                        if callBackCheckFunction != None:
                            use = callBackCheckFunction(dbEntry)

                        if use:
                            if mass > 0:
                                self.dbEntriesNeutral.append(dbEntry)
                            else:
                                self.dbEntriesMZ.append(dbEntry)
                            imported += 1

                    except Exception as ex:
                        logging.error("DB import error: Could not import row %d (%s)" % (rowi, ex.message))
                        notImported += 1

        logging.info(
            "Imported DB %s with %d entries (Current number of entries: %d)"
            % (
                dbName,
                len(self.dbEntriesMZ) + len(self.dbEntriesNeutral) - curEntriesCount,
                len(self.dbEntriesMZ) + len(self.dbEntriesNeutral),
            )
        )
        if notImported > 0:
            logging.error("Not imported %d entries (see above errors)" % (notImported))
        return imported, notImported

    def optimizeDB(self):
        self.dbEntriesNeutral = sorted(self.dbEntriesNeutral, key=lambda x: x.mass)
        self.dbEntriesMZ = sorted(self.dbEntriesMZ, key=lambda x: x.mz)

    def _findGeneric(self, list, getValue, valueLeft, valueRight):
        if len(list) == 0:
            return (-1, -1)

        min = 0
        max = len(list)

        while min < max and (max - min) > 1:
            cur = int(ceil((max + min) / 2.0))

            if valueLeft <= getValue(list[cur]) <= valueRight:
                leftBound = cur
                while leftBound > 0 and getValue(list[leftBound - 1]) >= valueLeft:
                    leftBound -= 1

                rightBound = cur
                while (rightBound + 1) < len(list) and getValue(list[rightBound + 1]) <= valueRight:
                    rightBound += 1

                return leftBound, rightBound

            if getValue(list[cur]) > valueRight:
                max = cur
            else:
                min = cur

        cur = min
        if valueLeft <= getValue(list[cur]) <= valueRight:
            leftBound = cur
            while leftBound > 0 and getValue(list[leftBound - 1]) >= valueLeft:
                leftBound -= 1

            rightBound = cur
            while (rightBound + 1) < len(list) and getValue(list[rightBound + 1]) <= valueRight:
                rightBound += 1

            return leftBound, rightBound

        cur = len(list) - 1
        if valueLeft <= getValue(list[cur]) <= valueRight:
            leftBound = cur
            while leftBound > 0 and getValue(list[leftBound - 1]) >= valueLeft:
                leftBound -= 1

            rightBound = cur
            while (rightBound + 1) < len(list) and getValue(list[rightBound + 1]) <= valueRight:
                rightBound += 1

            return leftBound, rightBound

        return -1, -1

    def searchDBForMZ(
        self,
        mz,
        polarity,
        charges,
        ppm,
        rt_min=None,
        rt_error=0.1,
        checkXN="Exact",
        element="C",
        Xn=1,
        adducts=[
            ("+H", 1.007276, "+", 1, 1),
            ("+NH4", 18.033823, "+", 1, 1),
            ("+Na", 22.989218, "+", 1, 1),
            ("+K", 38.963158, "+", 1, 1),
            ("+CH3OH+H", 33.033489, "+", 1, 1),
            ("-H", -1.007276, "-", 1, 1),
            ("+Na-2H", 20.974666, "-", 1, 1),
            ("+Cl", 34.969402, "-", 1, 1),
            ("+Br", 78.918885, "-", 1, 1),
            ("-2H+", -2 * 1.007276, "-", 1, 1),
        ],
    ):
        possibleHits = []

        if checkXN not in ["Don't use", "Exact", "Minimum"] and not checkXN.startswith("PlusMinus_"):
            logging.error("Unknown option '%s' for parameter checkXN, must be either of 'Don't use', 'Exact', 'Minimum' or 'PlusMinux_X' where X is a positive integer" % (checkXN))
            raise Exception("Unknown option '%s' for parameter checkXN, must be either of 'Don't use', 'Exact', 'Minimum' or 'PlusMinux_X' where X is a positive integer" % (checkXN))

        ## search for non-charged DB entries by subtracting putative adducts from the provided mz value
        for adduct in adducts:
            if polarity == adduct[2] and charges == adduct[3]:
                mass = (mz - adduct[1]) * adduct[3] / adduct[4]

                ph = self._findGeneric(
                    self.dbEntriesNeutral,
                    lambda x: x.mass,
                    mass - mass * ppm / 1000000.0,
                    mass + mass * ppm / 1000000.0,
                )
                if ph[0] != -1:
                    for entryi in range(ph[0], ph[1] + 1):
                        entry = self.dbEntriesNeutral[entryi]
                        if rt_min == None or entry.rt_min == None or (abs(rt_min - entry.rt_min) <= rt_error):
                            elems = None
                            if entry.sumFormula != "":
                                fT = formulaTools()
                                elems = fT.parseFormula(entry.sumFormula)
                            if (
                                checkXN == "Don't use"
                                or elems is None
                                or (checkXN == "Exact" and element in elems.keys() and elems[element] == Xn)
                                or (checkXN == "Minimum" and element in elems.keys() and elems[element] >= Xn)
                                or (checkXN.startswith("PlusMinus_") and element in elems.keys() and abs(elems[element] - Xn) <= int(checkXN[10 : len(checkXN)]))
                            ):
                                entry = deepcopy(entry)
                                entry.hitType = "MZ with %s to DB-mass match" % (adduct[0])
                                entry.matchErrorPPM = (mass - entry.mass) * 1e6 / mass
                                entry.matchErrorMass = mass - entry.mass
                                possibleHits.append(entry)

        ## search for charged DB entries by the provided mz value

        ph = self._findGeneric(
            self.dbEntriesMZ,
            lambda x: x.mz,
            mz - mz * ppm / 1000000.0,
            mz + mz * ppm / 1000000.0,
        )
        if ph[0] != -1:
            for entryi in range(ph[0], ph[1] + 1):
                entry = self.dbEntriesMZ[entryi]
                print(entry)
                if entry.polarity == polarity and (rt_min == None or entry.rt_min == None or (abs(rt_min - entry.rt_min) <= rt_error)):
                    elems = None
                    if entry.sumFormula != "":
                        fT = formulaTools()
                        elems = fT.parseFormula(entry.sumFormula)
                    if (
                        checkXN == "Don't use"
                        or elems is None
                        or (checkXN == "Exact" and element in elems.keys() and elems[element] == Xn)
                        or (checkXN == "Minimum" and element in elems.keys() and elems[element] >= Xn)
                        or (checkXN.startswith("PlusMinus_") and element in elems.keys() and abs(elems[element] - Xn) <= int(checkXN[10 : len(checkXN)]))
                    ):
                        entry = deepcopy(entry)
                        entry.hitType = "MZ to DB-MZ match"
                        entry.matchErrorPPM = (mz - entry.mz) * 1e6 / mz
                        entry.matchErrorMass = mz - entry.mz
                        possibleHits.append(entry)

        return possibleHits

    def searchDBForMass(
        self,
        mass,
        polarity,
        charges,
        ppm,
        rt_min=None,
        rt_error=0.1,
        checkXN="Exact",
        element="C",
        Xn=1,
        adducts=[
            ("+H", 1.007276, "+", 1, 1),
            ("+NH4", 18.033823, "+", 1, 1),
            ("+Na", 22.989218, "+", 1, 1),
            ("+K", 38.963158, "+", 1, 1),
            ("+CH3OH+H", 33.033489, "+", 1, 1),
            ("-H", -1.007276, "-", 1, 1),
            ("+Na-2H", 20.974666, "-", 1, 1),
            ("+Cl", 34.969402, "-", 1, 1),
            ("+Br", 78.918885, "-", 1, 1),
            ("-2H+", -2 * 1.007276, "-", 1, 1),
        ],
    ):
        possibleHits = []

        ## search for charged DB entries by adding adducts to the provided mass
        for adduct in adducts:
            if charges == adduct[3]:
                mz = mass * adduct[4] / adduct[3] + adduct[1]

                ph = self._findGeneric(
                    self.dbEntriesMZ,
                    lambda x: x.mz,
                    mz - mz * ppm / 1000000.0,
                    mz + mz * ppm / 1000000.0,
                )
                if ph[0] != -1:
                    for entryi in range(ph[0], ph[1] + 1):
                        entry = self.dbEntriesMZ[entryi]
                        if entry.polarity == adduct[2] and (rt_min == None or entry.rt_min == None or (abs(rt_min - entry.rt_min) <= rt_error)):
                            elems = None
                            if entry.sumFormula != "":
                                fT = formulaTools()
                                elems = fT.parseFormula(entry.sumFormula)
                            if (
                                checkXN == "Don't use"
                                or elems is None
                                or (checkXN == "Exact" and element in elems.keys() and elems[element] == Xn)
                                or (checkXN == "Minimum" and element in elems.keys() and elems[element] >= Xn)
                                or (checkXN.startswith("PlusMinus_") and element in elems.keys() and abs(elems[element] - Xn) <= int(checkXN[10 : len(checkXN)]))
                            ):
                                entry = deepcopy(entry)
                                entry.hitType = "Mass to DB-MZ match with %s" % (adduct[0])
                                entry.matchErrorPPM = (mz - entry.mz) * 1e6 / mz
                                entry.matchErrorMass = mz - entry.mz
                                possibleHits.append(entry)

        ## search for non-charged DB entries by the provided mass
        ph = self._findGeneric(
            self.dbEntriesNeutral,
            lambda x: x.mass,
            mass - mass * ppm / 1000000.0,
            mass + mass * ppm / 1000000.0,
        )
        if ph[0] != -1:
            for entryi in range(ph[0], ph[1] + 1):
                entry = self.dbEntriesNeutral[entryi]
                if rt_min == None or entry.rt_min == None or (abs(rt_min - entry.rt_min) <= rt_error):
                    elems = None
                    if entry.sumFormula != "":
                        fT = formulaTools()
                        elems = fT.parseFormula(entry.sumFormula)
                    if (
                        checkXN == "Don't use"
                        or elems is None
                        or (checkXN == "Exact" and element in elems.keys() and elems[element] == Xn)
                        or (checkXN == "Minimum" and element in elems.keys() and elems[element] >= Xn)
                        or (checkXN.startswith("PlusMinus_") and element in elems.keys() and abs(elems[element] - Xn) <= int(checkXN[10 : len(checkXN)]))
                    ):
                        entry = deepcopy(entry)
                        entry.hitType = "Mass to DB-Mass match"
                        entry.matchErrorPPM = (mass - entry.mass) * 1e6 / mass
                        entry.matchErrorMass = mass - entry.mass
                        possibleHits.append(entry)

        return possibleHits

    def searchDB(
        self,
        mass,
        mz,
        polarity,
        charges,
        rt_min,
        ppm=5.0,
        rt_error=0.1,
        checkXN="Exact",
        element="C",
        Xn=1,
        adducts=[
            ("+H", 1.007276, "+", 1, 1),
            ("+NH4", 18.033823, "+", 1, 1),
            ("+Na", 22.989218, "+", 1, 1),
            ("+K", 38.963158, "+", 1, 1),
            ("+CH3OH+H", 33.033489, "+", 1, 1),
            ("-H", -1.007276, "-", 1, 1),
            ("+Na-2H", 20.974666, "-", 1, 1),
            ("+Cl", 34.969402, "-", 1, 1),
            ("+Br", 78.918885, "-", 1, 1),
            ("-2H+", -2 * 1.007276, "-", 1, 1),
        ],
    ):
        if mass != None:
            return self.searchDBForMass(
                mass,
                polarity,
                charges,
                ppm,
                rt_min,
                rt_error,
                checkXN,
                element,
                Xn,
                adducts,
            )
        else:
            return self.searchDBForMZ(
                mz,
                polarity,
                charges,
                ppm,
                rt_min,
                rt_error,
                checkXN,
                element,
                Xn,
                adducts,
            )


if False and __name__ == "__main__":
    db = DBSearch()
    db.addEntriesFromFile("SomeMets", "N:/iBAM/Christoph/Maria/DB_BiolPaper_Mets.txt")
    db.optimizeDB()

    print("hits from mass obtained from negative ion mode")
    for hit in db.searchDB(
        mass=290.1379484,
        mz=290.1379484 - 1.007276,
        polarity="-",
        charges=1,
        rt_min=None,
        ppm=5,
        rt_error=0.1,
        checkXN="Minimal",
        element="C",
        Xn=9,
    ):
        print(hit.sumFormula, hit.hitType, hit.name, str(hit.additionalInfo), str(hit))
    print("hits from mass obtained from positive ion mode")
    for hit in db.searchDB(
        mass=290.1379484,
        mz=290.1379484 - 1.007276,
        polarity="+",
        charges=1,
        rt_min=None,
        ppm=5,
        rt_error=0.1,
        checkXN="Minimal",
        element="C",
        Xn=9,
    ):
        print(hit.sumFormula, hit.hitType, hit.name, str(hit.additionalInfo), str(hit))

    print("mz hits obtained from neg mode ion")
    for hit in db.searchDB(
        mass=None,
        mz=290.1379484 - 1.007276,
        polarity="-",
        charges=1,
        rt_min=None,
        ppm=5,
        rt_error=0.1,
        checkXN="Minimal",
        element="C",
        Xn=9,
    ):
        print(hit.sumFormula, hit.hitType, hit.name, str(hit.additionalInfo), str(hit))

    print("mz hits obtained from pos mode ion")
    for hit in db.searchDB(
        mass=None,
        mz=290.1379484 + 1.007276,
        polarity="+",
        charges=1,
        rt_min=None,
        ppm=5,
        rt_error=0.1,
        checkXN="Minimal",
        element="C",
        Xn=9,
    ):
        print(hit.sumFormula, hit.hitType, hit.name, str(hit.additionalInfo), str(hit))


if True and __name__ == "__main__":
    db = DBSearch()
    # db.addEntriesFromFile("SOS", "C:/Users/cbueschl/Desktop/me2db_sos_db_250214.tsv")
    # db.addEntriesFromFile("primary_metabolites", "C:/Users/cbueschl/Desktop/me2db_primary_metabolites.tsv")
    # db.addEntriesFromFile("flaxCyc", "C:/Users/cbueschl/Desktop/me2db_flaxCyc.tsv")
    db.addEntriesFromFile("KEGG", "C:/development/MetaboliteDBs/DBsForMetExtractII/KEGG_compounds.tsv")
    db.optimizeDB()

    mass = None
    mz = 119.0503171

    polarity = "-"
    charges = 1

    hits = db.searchDB(
        mass=mass,
        mz=mz,
        polarity=polarity,
        charges=charges,
        rt_min=None,
        ppm=5,
        checkXN="Don't use",
    )

    import sys

    sys.stdout.flush()
    print("\n\n")
    print("DBName\tNum\tName\tChemicalFormula\tErrorPPM\tAdditionalInfo")
    for hit in hits:
        print(
            "%s\t%s\t%s\t%s\t%s\t%s"
            % (
                hit.dbName,
                hit.num,
                hit.name,
                hit.sumFormula,
                hit.matchErrorPPM,
                hit.additionalInfo,
            )
        )
