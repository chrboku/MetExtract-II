import sys
sys.path.append("C:/PyMetExtract/PyMetExtract")

from formulaTools import formulaTools


import TableUtils

import csv

from copy import deepcopy
from math import floor


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
        self.rt_min=rt_min
        self.mz=mz
        self.polarity=polarity
        self.additionalInfo=additionalInfo

        self.hitType=hitType




class DBSearch:

    def __init__(self):
        self.dbEntriesNeutral=[]
        self.dbEntriesMZ=[]


    def addEntriesFromFile(self, dbName, dbFile):
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
                    num=row[headers["Num"]]
                    name=row[headers["Name"]]
                    sumFormula=row[headers["SumFormula"]]
                    rt_min=float(row[headers["Rt_min"]]) if row[headers["Rt_min"]]!="" else None
                    mz=float(row[headers["MZ"]]) if row[headers["MZ"]]!="" else None
                    polarity=row[headers["IonisationMode"]]
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
                            logging.error("DB import error: The sumformula (%s) of the entry %s could not be parsed"%(sumFormula, num))
                            continue

                    dbEntry=DBEntry(dbName, num, name, sumFormula, mass, rt_min, mz, polarity, additionalInfo)

                    if mass>0:
                        self.dbEntriesNeutral.append(dbEntry)
                    else:
                        self.dbEntriesMZ.append(dbEntry)

    def optimizeDB(self):
        self.dbEntriesNeutral=sorted(self.dbEntriesNeutral, key=lambda x: x.mass)
        self.dbEntriesMZ=sorted(self.dbEntriesMZ, key=lambda x: x.mz)

    def _findGeneric(self, list, getValue, valueLeft, valueRight):
        if len(list) == 0:
            return (-1, -1)

        min = 0
        max = len(list)

        while min < max:
            cur = int(floor((max + min) / 2.))

            if valueLeft <= getValue(list[cur]) <= valueRight:
                leftBound = cur
                while leftBound > 0 and getValue(list[leftBound - 1]) >= valueLeft:
                    leftBound -= 1

                rightBound = cur
                while (rightBound + 1) < len(list) and getValue(list[rightBound + 1]) <= valueRight:
                    rightBound += 1

                return leftBound, rightBound

            if getValue(list[cur]) > valueRight:
                max = cur - 1
            else:
                min = cur + 1

        return -1, -1

    def searchDBForMZ(self, mz, polarity, charges, ppm, rt_min=None, rt_error=0.1, checkXN="Exact", element="C", Xn=1, adducts=
    [('+H', 1.007276, '+', 1), ('+NH4', 18.033823, '+', 1),
     ('+Na', 22.989218, '+', 1), ('+K', 38.963158, '+', 1),
     ('-H', -1.007276, '-', 1), ('+Na-2H', 20.974666, '-', 1),
     ('+Cl', 34.969402, '-', 1), ('+Br', 78.918885, '-', 1),
     ('-2H+', -2*1.007276, '-', 1)]
                      ):
        possibleHits=[]

        ## search for non-charged DB entries by subtracting putative adducts from the provided mz value
        for adduct in adducts:
            if polarity==adduct[2] and charges==adduct[3]:
                mass=mz-adduct[1]
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
                                possibleHits.append(entry)


        ## search for charged DB entries by the provided mz value
        ph=self._findGeneric(self.dbEntriesMZ, lambda x: x.mz, mz-mz*ppm/1000000., mz+mz*ppm/1000000.)
        if ph[0]!=-1:
            for entryi in ph:
                entry=self.dbEntriesMZ[entryi]
                if rt_min==None or entry.rt_min==None or (abs(rt_min-entry.rt_min)<=rt_error):
                    elems=None
                    if entry.sumFormula!="":
                        fT=formulaTools()
                        elems=fT.parseFormula(entry.sumFormula)
                    if checkXN=="Don't use" or elems is None or (checkXN=="Exact" and elems[element]==Xn) or (checkXN=="Minimum" and elems[element]>=Xn) or (checkXN.startswith("PlusMinus_") and abs(elems[element]-Xn)<=int(checkXN[10:len(checkXN)])):
                        entry=deepcopy(entry)
                        entry.hitType="MZ match"
                        possibleHits.append(entry)

        return possibleHits


    def searchDBForMass(self, mass, polarity, charges, ppm, rt_min=None, rt_error=0.1, checkXN="Exact", element="C", Xn=1, adducts=
    [('+H', 1.007276, '+', 1), ('+NH4', 18.033823, '+', 1),
     ('+Na', 22.989218, '+', 1), ('+K', 38.963158, '+', 1), ('+CHOH+H', 33.033489, '+', 1),
     ('-H', -1.007276, '-', 1), ('+Na-2H', 20.974666, '-', 1),
     ('+Cl', 34.969402, '-', 1), ('+Br', 78.918885, '-', 1),
     ('-2H+', -2 * 1.007276, '-', 1)]
                      ):
        possibleHits = []

        ## search for charged DB entries by adding adducts to the provided mass
        for adduct in adducts:
            if polarity == adduct[2] and charges == adduct[3]:
                mz = mass + adduct[1]
                ph = self._findGeneric(self.dbEntriesMZ, lambda x: x.mz, mz - mz * ppm / 1000000.,
                                       mz + mz * ppm / 1000000.)
                if ph[0]!=-1:
                    for entryi in ph:
                        entry=self.dbEntriesMZ[entryi]
                        if rt_min==None or entry.rt_min==None or (abs(rt_min - entry.rt_min) <= rt_error):
                            elems=None
                            if entry.sumFormula!="":
                                fT=formulaTools()
                                elems=fT.parseFormula(entry.sumFormula)
                            if checkXN=="Don't use" or elems is None or (checkXN=="Exact" and elems[element]==Xn) or (checkXN=="Minimum" and elems[element]>=Xn) or (checkXN.startswith("PlusMinus_") and abs(elems[element]-Xn)<=int(checkXN[10:len(checkXN)])):
                                entry = deepcopy(entry)
                                entry.hitType = "Calc. Adduct: %s" % adduct[0]
                                possibleHits.append(entry)

        ## search for non-charged DB entries by the provided mass
        ph = self._findGeneric(self.dbEntriesNeutral, lambda x: x.mass, mass - mass * ppm / 1000000., mass + mass * ppm / 1000000.)
        if ph[0]!=-1:
            for entryi in ph:
                entry=self.dbEntriesNeutral[entryi]
                if rt_min==None or entry.rt_min==None or (abs(rt_min - entry.rt_min) <= rt_error):
                    elems=None
                    if entry.sumFormula!="":
                        fT=formulaTools()
                        elems=fT.parseFormula(entry.sumFormula)
                    if checkXN=="Don't use" or elems is None or (checkXN=="Exact" and elems[element]==Xn) or (checkXN=="Minimum" and elems[element]>=Xn) or (checkXN.startswith("PlusMinus_") and abs(elems[element]-Xn)<=int(checkXN[10:len(checkXN)])):
                        entry = deepcopy(entry)
                        entry.hitType = "M match"
                        possibleHits.append(entry)

        return possibleHits



    def searchDB(self, mass, mz, polarity, charges, rt_min, ppm=5., rt_error=0.1, checkXN="Exact", element="C", Xn=1, adducts=
    [('+H', 1.007276, '+', 1), ('+NH4', 18.033823, '+', 1),
     ('+Na', 22.989218, '+', 1), ('+K', 38.963158, '+', 1), ('+CHOH+H', 33.033489, '+', 1),
     ('-H', -1.007276, '-', 1), ('+Na-2H', 20.974666, '-', 1),
     ('+Cl', 34.969402, '-', 1), ('+Br', 78.918885, '-', 1),
     ('-2H+', -2 * 1.007276, '-', 1)]):
        if mass!=None:
            return self.searchDBForMass(mass, polarity, charges, ppm, rt_min, rt_error, checkXN, element, Xn, adducts)
        else:
            return self.searchDBForMZ(mz, polarity, charges, ppm, rt_min, rt_error, checkXN, element, Xn, adducts)


if __name__=="__main__":
    db=DBSearch()

    db.addEntriesFromFile("SomeMets", "C:/PyMetExtract/implTest/TestMetabolites.tsv")
    db.optimizeDB()




    print "hits"
    for hit in db.searchDB(mass=None, mz=297.132612435, polarity="+", charges=1, rt_min=None, ppm=150, rt_error=0.1,
                           checkXN="Exact", element="C", Xn=15):
        print hit.name, hit.hitType
