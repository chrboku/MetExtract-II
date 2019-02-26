import sys
sys.path.append("C:/PyMetExtract/PyMetExtract")

from formulaTools import formulaTools


import TableUtils

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
import LoggingSetup
LoggingSetup.LoggingSetup.Instance().initLogging()



class DBEntry:
    def __init__(self, dbName, num, name, sumFormula, mass, rt_min, mz, polarity, additionalInfo, hitType=""):
        self.dbName=dbName
        self.num=num
        self.name=name
        self.sumFormula=sumFormula
        self.mass=mass
        self.matchErrorPPM=-1
        self.rt_min=rt_min
        self.mz=mz
        self.polarity=polarity
        self.additionalInfo=additionalInfo

        self.hitType=hitType




class DBSearch:

    def __init__(self):
        self.dbEntriesNeutral=[]
        self.dbEntriesMZ=[]


    def addEntriesFromFile(self, dbName, dbFile, callBackCheckFunction=None):
        fT=formulaTools()

        with open(dbFile) as fIn:
            tsvin=csv.reader(fIn, delimiter="\t")

            headers={}
            for rowi, row in enumerate(tsvin):
                if row[0].startswith("#"):
                    continue
                if rowi==0:
                    for j in range(len(row)):
                        headers[row[j]]=j

                else:
                    try:
                        num=row[headers["Num"]].strip()
                        name=row[headers["Name"]].strip().replace("\"","DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n","RETTURN").replace("\r","CarrRETURN").replace("#","HASH")
                        sumFormula=row[headers["SumFormula"]].strip()
                        rt_min=float(row[headers["Rt_min"]]) if row[headers["Rt_min"]]!="" else None
                        mz=float(row[headers["MZ"]]) if row[headers["MZ"]]!="" else None
                        polarity=row[headers["IonisationMode"]].strip()
                        additionalInfo={}
                        for header in headers.keys():
                            if header not in ["Num", "Name", "SumFormula", "Rt_min", "MZ", "IonisationMode"]:
                                additionalInfo[header]=row[headers[header]]

                        mass=0
                        if sumFormula!="":
                            try:
                                elems=fT.parseFormula(sumFormula)
                                mass=fT.calcMolWeight(elems)
                            except Exception as ex:
                                logging.error("DB import error (%s, row: %d): The sumformula (%s) of the entry %s '%s' could not be parsed"%(dbName, rowi, sumFormula, num, name))
                                continue

                        dbEntry=DBEntry(dbName, num, name, sumFormula, mass, rt_min, mz, polarity, additionalInfo)

                        use=True
                        if callBackCheckFunction!=None:
                            use=callBackCheckFunction(dbEntry)

                        if use:
                            if mass>0:
                                self.dbEntriesNeutral.append(dbEntry)
                            else:
                                self.dbEntriesMZ.append(dbEntry)
                    except Exception as ex:
                        logging.error("DB import error: Could not import row %d (%s)"%(rowi, ex.message))
                        continue
        print "Imported DB %s (Current number of entries: %d)"%(dbName, len(self.dbEntriesMZ)+len(self.dbEntriesNeutral))


    def optimizeDB(self):
        self.dbEntriesNeutral=sorted(self.dbEntriesNeutral, key=lambda x: x.mass)
        self.dbEntriesMZ=sorted(self.dbEntriesMZ, key=lambda x: x.mz)

    def _findGeneric(self, list, getValue, valueLeft, valueRight):
        if len(list) == 0:
            return (-1, -1)

        min = 0
        max = len(list)


        while min < max and (max-min)>1:
            cur = int(ceil((max + min) / 2.))

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

        cur=min
        if valueLeft <= getValue(list[cur]) <= valueRight:
            leftBound = cur
            while leftBound > 0 and getValue(list[leftBound - 1]) >= valueLeft:
                leftBound -= 1

            rightBound = cur
            while (rightBound + 1) < len(list) and getValue(list[rightBound + 1]) <= valueRight:
                rightBound += 1

            return leftBound, rightBound

        cur = len(list)-1
        if valueLeft <= getValue(list[cur]) <= valueRight:
            leftBound = cur
            while leftBound > 0 and getValue(list[leftBound - 1]) >= valueLeft:
                leftBound -= 1

            rightBound = cur
            while (rightBound + 1) < len(list) and getValue(list[rightBound + 1]) <= valueRight:
                rightBound += 1

            return leftBound, rightBound


        return -1, -1

    def searchDBForMZ(self, mz, polarity, charges, ppm, rt_min=None, rt_error=0.1, checkXN="Exact", element="C", Xn=1, adducts=
    [('+H', 1.007276, '+', 1, 1), ('+NH4', 18.033823, '+', 1, 1),
     ('+Na', 22.989218, '+', 1, 1), ('+K', 38.963158, '+', 1, 1), ('+CH3OH+H', 33.033489, '+', 1, 1),
     ('-H', -1.007276, '-', 1, 1), ('+Na-2H', 20.974666, '-', 1, 1),
     ('+Cl', 34.969402, '-', 1, 1), ('+Br', 78.918885, '-', 1, 1),
     ('-2H+', -2 * 1.007276, '-', 1, 1)]):
        possibleHits=[]

        ## search for non-charged DB entries by subtracting putative adducts from the provided mz value
        for adduct in adducts:
            if polarity==adduct[2] and charges==adduct[3]:
                mass=(mz - adduct[1])*adduct[3]/adduct[4]

                ph=self._findGeneric(self.dbEntriesNeutral, lambda x: x.mass, mass-mass*ppm/1000000., mass+mass*ppm/1000000.)
                if ph[0]!=-1:
                    for entryi in range(ph[0], ph[1]+1):
                        entry = self.dbEntriesNeutral[entryi]
                        if rt_min==None or entry.rt_min==None or (abs(rt_min-entry.rt_min)<=rt_error):
                            elems=None
                            if entry.sumFormula!="":
                                fT=formulaTools()
                                elems=fT.parseFormula(entry.sumFormula)
                            if checkXN=="Don't use" or elems is None or (checkXN=="Exact" and elems[element]==Xn) or (checkXN=="Minimum" and elems[element]>=Xn) or (checkXN.startswith("PlusMinus_") and abs(elems[element]-Xn)<=int(checkXN[10:len(checkXN)])):
                                entry=deepcopy(entry)
                                entry.hitType="Calc. Adduct: %s"%adduct[0]
                                entry.matchErrorPPM=(mass-entry.mass)*1E6/mz
                                possibleHits.append(entry)


        ## search for charged DB entries by the provided mz value
        ph=self._findGeneric(self.dbEntriesMZ, lambda x: x.mz, mz-mz*ppm/1000000., mz+mz*ppm/1000000.)
        if ph[0]!=-1:
            for entryi in range(ph[0], ph[1]+1):
                entry=self.dbEntriesMZ[entryi]
                if rt_min==None or entry.rt_min==None or (abs(rt_min-entry.rt_min)<=rt_error):
                    elems=None
                    if entry.sumFormula!="":
                        fT=formulaTools()
                        elems=fT.parseFormula(entry.sumFormula)
                    if checkXN=="Don't use" or elems is None or (checkXN=="Exact" and elems[element]==Xn) or (checkXN=="Minimum" and elems[element]>=Xn) or (checkXN.startswith("PlusMinus_") and abs(elems[element]-Xn)<=int(checkXN[10:len(checkXN)])):
                        entry=deepcopy(entry)
                        entry.hitType="MZ match"
                        entry.matchErrorPPM=(mz-entry.mz)*1E6/mz
                        possibleHits.append(entry)

        return possibleHits


    def searchDBForMass(self, mass, polarity, charges, ppm, rt_min=None, rt_error=0.1, checkXN="Exact", element="C", Xn=1, adducts=
    [('+H', 1.007276, '+', 1, 1), ('+NH4', 18.033823, '+', 1, 1),
     ('+Na', 22.989218, '+', 1, 1), ('+K', 38.963158, '+', 1, 1), ('+CH3OH+H', 33.033489, '+', 1, 1),
     ('-H', -1.007276, '-', 1, 1), ('+Na-2H', 20.974666, '-', 1, 1),
     ('+Cl', 34.969402, '-', 1, 1), ('+Br', 78.918885, '-', 1, 1),
     ('-2H+', -2 * 1.007276, '-', 1, 1)]):
        possibleHits = []

        ## search for charged DB entries by adding adducts to the provided mass
        for adduct in adducts:
            if polarity == adduct[2] and charges == adduct[3]:
                mz = mass * adduct[4]/adduct[3]-adduct[1]

                ph = self._findGeneric(self.dbEntriesMZ, lambda x: x.mz, mz - mz * ppm / 1000000.,
                                       mz + mz * ppm / 1000000.)
                if ph[0]!=-1:
                    for entryi in range(ph[0], ph[1]+1):
                        entry=self.dbEntriesMZ[entryi]
                        if rt_min==None or entry.rt_min==None or (abs(rt_min - entry.rt_min) <= rt_error):
                            elems=None
                            if entry.sumFormula!="":
                                fT=formulaTools()
                                elems=fT.parseFormula(entry.sumFormula)
                            if checkXN=="Don't use" or elems is None or (checkXN=="Exact" and elems[element]==Xn) or (checkXN=="Minimum" and elems[element]>=Xn) or (checkXN.startswith("PlusMinus_") and abs(elems[element]-Xn)<=int(checkXN[10:len(checkXN)])):
                                entry = deepcopy(entry)
                                entry.hitType = "Calc. Adduct: %s" % adduct[0]
                                entry.matchErrorPPM=(mz-entry.mz)*1E6/mz
                                possibleHits.append(entry)

        ## search for non-charged DB entries by the provided mass
        ph = self._findGeneric(self.dbEntriesNeutral, lambda x: x.mass, mass - mass * ppm / 1000000., mass + mass * ppm / 1000000.)
        if ph[0]!=-1:
            for entryi in range(ph[0], ph[1]+1):
                entry=self.dbEntriesNeutral[entryi]
                if rt_min==None or entry.rt_min==None or (abs(rt_min - entry.rt_min) <= rt_error):
                    elems=None
                    if entry.sumFormula!="":
                        fT=formulaTools()
                        elems=fT.parseFormula(entry.sumFormula)
                    if checkXN=="Don't use" or elems is None or (checkXN=="Exact" and elems[element]==Xn) or (checkXN=="Minimum" and elems[element]>=Xn) or (checkXN.startswith("PlusMinus_") and abs(elems[element]-Xn)<=int(checkXN[10:len(checkXN)])):
                        entry = deepcopy(entry)
                        entry.hitType = "M match"
                        entry.matchErrorPPM=(mass-entry.mass)*1E6/mass
                        possibleHits.append(entry)

        return possibleHits



    def searchDB(self, mass, mz, polarity, charges, rt_min, ppm=5., rt_error=0.1, checkXN="Exact", element="C", Xn=1, adducts=
    [('+H', 1.007276, '+', 1, 1), ('+NH4', 18.033823, '+', 1, 1),
     ('+Na', 22.989218, '+', 1, 1), ('+K', 38.963158, '+', 1, 1), ('+CH3OH+H', 33.033489, '+', 1, 1),
     ('-H', -1.007276, '-', 1, 1), ('+Na-2H', 20.974666, '-', 1, 1),
     ('+Cl', 34.969402, '-', 1, 1), ('+Br', 78.918885, '-', 1, 1),
     ('-2H+', -2 * 1.007276, '-', 1, 1)]):
        if mass!=None:
            return self.searchDBForMass(mass, polarity, charges, ppm, rt_min, rt_error, checkXN, element, Xn, adducts)
        else:
            return self.searchDBForMZ(mz, polarity, charges, ppm, rt_min, rt_error, checkXN, element, Xn, adducts)


if False and __name__=="__main__":
    db=DBSearch()

    db.addEntriesFromFile("SomeMets", "C:/PyMetExtract/implTest/TestMetabolites.tsv")
    db.optimizeDB()




    print "hits"
    for hit in db.searchDB(mass=None, mz=297.132612435, polarity="+", charges=1, rt_min=None, ppm=150, rt_error=0.1,
                           checkXN="Exact", element="C", Xn=15):
        print hit.name, hit.hitType




if False and __name__=="__main__":
    ppm=5.
    rtError=0.5

    db=DBSearch()

    db.addEntriesFromFile("PlantCyc13", "C:\PyMetExtract\MetaboliteDBs\DBsForMetExtractII\plantCyc13_export.tsv")
    db.optimizeDB()

    db2=DBSearch()

    db2.addEntriesFromFile("CCD", "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/CCD_Dec_2013.tsv")
    db2.optimizeDB()

    with open("E:/FH/141107_158_maize_FV/results.OGrp.tsv", "rb") as fin:
        with open("E:/FH/141107_158_maize_FV/results.OGrp.newDBSearch.tsv", "wb") as fout:

            csvReader=csv.reader(fin, delimiter="\t")
            csvWriter=csv.writer(fout, delimiter="\t")

            headers={}
            for rowi, row in enumerate(csvReader):

                if rowi==0:
                    for celli, cell in enumerate(row):
                        headers[cell]=celli
                    csvWriter.writerow(row+["PlantCyc13", "CCD_13"])
                elif not row[0].startswith("#"):
                    try:
                        mz=float(row[headers["MZ"]])
                        rt=float(row[headers["RT"]])
                        xn=int(row[headers["Xn"]])
                        ionMode=row[headers["Ionisation_Mode"]]
                        charge=int(row[headers["Charge"]])

                        cnMatch = "Exact"

                        hits=[]
                        for hit in db.searchDB(mass=None, mz=mz, polarity=ionMode, charges=charge, rt_min=None, ppm=ppm,
                                               rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                            hits.append(hit)

                        hits2=[]
                        for hit in db2.searchDB(mass=None, mz=mz, polarity=ionMode, charges=charge, rt_min=None, ppm=ppm,
                                               rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                            hits2.append(hit)

                        csvWriter.writerow(row+\
                                           [";".join(["%s {Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits])]+\
                                           [";".join(["%s {Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits2])])

                    except Exception as ex:
                        print ex.message
                        csvWriter.writerow(row + ["", ""])
                else:
                    csvWriter.writerow(row)

            csvWriter.writerow(["###### Database search"])
            csvWriter.writerow(["## ppm=%.1f"%(ppm)])

    with open("E:/FH/141107_158_maize_FV/results.OGrp.mInt.tsv", "rb") as fin:
        with open("E:/FH/141107_158_maize_FV/results.OGrp.mInt.newDBSearch.tsv", "wb") as fout:

            csvReader=csv.reader(fin, delimiter="\t")
            csvWriter=csv.writer(fout, delimiter="\t")

            headers={}
            for rowi, row in enumerate(csvReader):

                if rowi==0:
                    for celli, cell in enumerate(row):
                        headers[cell]=celli
                    csvWriter.writerow(row+["PlantCyc13", "CCD_13"])
                elif not row[0].startswith("#"):
                    try:
                        mz=float(row[headers["MZ"]])
                        rt=float(row[headers["RT"]])
                        xn=int(row[headers["Xn"]])
                        ionMode=row[headers["Ionisation_Mode"]]
                        charge=int(row[headers["Charge"]])

                        cnMatch = "Exact"

                        hits=[]
                        for hit in db.searchDB(mass=None, mz=mz, polarity=ionMode, charges=charge, rt_min=None, ppm=ppm,
                                               rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                            hits.append(hit)

                        hits2=[]
                        for hit in db2.searchDB(mass=None, mz=mz, polarity=ionMode, charges=charge, rt_min=None, ppm=ppm,
                                               rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                            hits2.append(hit)

                        csvWriter.writerow(row+\
                                           [";".join(["%s {Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits])]+\
                                           [";".join(["%s {Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits2])])

                    except Exception as ex:
                        print ex.message
                        csvWriter.writerow(row + ["", ""])
                else:
                    csvWriter.writerow(row)

            csvWriter.writerow(["###### Database search"])
            csvWriter.writerow(["## ppm=%.1f"%(ppm)])




if True and __name__=="__main__":
    ppm=5.
    rtError=0.5

    db=DBSearch()

    #db.addEntriesFromFile("DB_Auxins", "H:/190117_531_TAM_feeded Fg_Thomas/DBs/DB_Auxins_SF_RF_MZ.tsv")
    db.addEntriesFromFile("DBChemSpider_Trp_metabolites", "H:/190117_531_TAM_feeded Fg_Thomas/DBs/DBChemspider_Trp_metabolites.tsv")
    db.optimizeDB()
