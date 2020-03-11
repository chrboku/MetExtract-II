import sys
sys.path.append("C:/PyMetExtract/PyMetExtract")

from formulaTools import formulaTools
from utils import is_float

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
        self.matchErrorMass=-1
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
        imported=0
        notImported=0

        curEntriesCount=len(self.dbEntriesMZ)+len(self.dbEntriesNeutral)

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
                        num=row[headers["Num"]].strip().replace("\"","DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n","RETURN").replace("\r","CarrRETURN").replace("#","HASH")
                        name=row[headers["Name"]].strip().replace("\"","DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n","RETURN").replace("\r","CarrRETURN").replace("#","HASH")
                        sumFormula=row[headers["SumFormula"]].strip().replace("\"","DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n","RETURN").replace("\r","CarrRETURN").replace("#","HASH")
                        rt_min=float(row[headers["Rt_min"]]) if row[headers["Rt_min"]]!="" else None
                        mz=float(row[headers["MZ"]]) if row[headers["MZ"]]!="" else None
                        polarity=row[headers["IonisationMode"]].strip().replace("\"","DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n","RETURN").replace("\r","CarrRETURN").replace("#","HASH")
                        additionalInfo={}
                        for header in headers.keys():
                            if header not in ["Num", "Name", "SumFormula", "Rt_min", "MZ", "IonisationMode"]:
                                additionalInfo[header]=row[headers[header]].replace("\"","DOURBLEPRIME").replace("'", "PRIME").replace("\t", "TAB").replace("\n","RETURN").replace("\r","CarrRETURN").replace("#","HASH")

                        mass=0
                        if sumFormula!="":
                            try:
                                elems=fT.parseFormula(sumFormula)
                                mass=fT.calcMolWeight(elems)
                            except Exception as ex:
                                logging.error("DB import error (%s, row: %d): The sumformula (%s) of the entry %s '%s' could not be parsed"%(dbName, rowi, sumFormula, num, name))
                                notImported+=1

                        dbEntry=DBEntry(dbName, num, name, sumFormula, mass, rt_min, mz, polarity, additionalInfo)

                        use=True
                        if callBackCheckFunction!=None:
                            use=callBackCheckFunction(dbEntry)

                        if use:
                            if mass>0:
                                self.dbEntriesNeutral.append(dbEntry)
                            else:
                                self.dbEntriesMZ.append(dbEntry)
                            imported+=1

                    except Exception as ex:
                        logging.error("DB import error: Could not import row %d (%s)"%(rowi, ex.message))
                        notImported+=1

        print "Imported DB %s with %d entries (Current number of entries: %d)"%(dbName, len(self.dbEntriesMZ)+len(self.dbEntriesNeutral)-curEntriesCount, len(self.dbEntriesMZ)+len(self.dbEntriesNeutral))
        return imported, notImported

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
                                entry.hitType=adduct[0]
                                entry.matchErrorPPM=(mass-entry.mass)*1E6/mass
                                entry.matchErrorMass=mass-entry.mass
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
                        entry.matchErrorMass=mz-entry.mz
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
                                entry.hitType = adduct[0]
                                entry.matchErrorPPM=(mz-entry.mz)*1E6/mz
                                entry.matchErrorMass=mz-entry.mz
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
                        entry.matchErrorMass=mass-entry.mass
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

    db.addEntriesFromFile("Antibase2017",                              "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017.tsv")
    db.addEntriesFromFile("Antibase2017_Gliocladium",                  "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Gliocladium.tsv")
    db.addEntriesFromFile("Antibase2017_Trichoderma",                  "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Trichoderma.tsv")
    db.addEntriesFromFile("Antibase2017_Trichorzianins",               "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/Trichorzianins.tsv")

    db.addEntriesFromFile("Antibase2017_Trichoderma_viride",           "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Trichoderma_viride.tsv")
    db.addEntriesFromFile("Antibase2017_Trichoderma_asperellum",       "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Trichoderma_asperellum.tsv")
    db.addEntriesFromFile("Antibase2017_Trichoderma_atroviride",       "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Trichoderma_atroviride.tsv")

    db.addEntriesFromFile("Antibase2017_Trichoderma_hamatum",          "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Trichoderma_hamatum.tsv")
    db.addEntriesFromFile("Antibase2017_Trichoderma_harzianum",        "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Trichoderma_harzianum.tsv")
    db.addEntriesFromFile("Antibase2017_Trichoderma_koningii",         "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Trichoderma_koningii.tsv")

    db.addEntriesFromFile("Antibase2017_Trichoderma_longibrachiatum",  "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Trichoderma_longibrachiatum.tsv")
    db.addEntriesFromFile("Antibase2017_Trichoderma_reesei",           "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Trichoderma_reesei.tsv")
    db.addEntriesFromFile("Antibase2017_Trichoderma_virens",           "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/AntiBase_2017_Trichoderma_virens.tsv")
    db.optimizeDB()

    dbs=["DBs_Antibase2017", "DBs_Antibase2017_Gliocladium", "DBs_Antibase2017_Trichoderma", "DBs_Antibase2017_Trichorzianins",
         "DBs_Antibase2017_Trichoderma_viride","DBs_Antibase2017_Trichoderma_asperellum","DBs_Antibase2017_Trichoderma_atroviride",
         "DBs_Antibase2017_Trichoderma_hamatum","DBs_Antibase2017_Trichoderma_harzianum","DBs_Antibase2017_Trichoderma_koningii",
         "DBs_Antibase2017_Trichoderma_longibrachiatum","DBs_Antibase2017_Trichoderma_reesei","DBs_Antibase2017_Trichoderma_virens"]


    annoFiles=[("H:/180829_435_Trichoderma_TOR_SusanneZeilinger_MainExperiment/EVAL_InitialDataEval13C/results_noDB.tsv",
                "C:/temp/results_mainExp_SFs_DBs_temp.tsv"),
               ("H:/180829_435_Trichoderma_TOR_SusanneZeilinger_MainExperiment/EVAL_ComparisonMedia/results.tsv",
                "C:/temp/results_mediaComparsison_SFs_DBs_temp.tsv"),
               ("H:/180907_445_Trichoderma_TOR_Glutamine_SideExperiment/results.tsv",
                "C:/temp/results_GlutamineSideExperiment_SFs_DBs_temp.tsv")]

    for annoFile in annoFiles:
        with open(annoFile[0], "rb") as fin:
            with open(annoFile[1], "wb") as fout:

                csvReader=csv.reader(fin, delimiter="\t")
                csvWriter=csv.writer(fout, delimiter="\t")

                headers={}
                for rowi, row in enumerate(csvReader):

                    if rowi==0:
                        for celli, cell in enumerate(row):
                            headers[cell]=celli
                        csvWriter.writerow(row+[d+"_count" for d in dbs]+[d for d in dbs])
                    elif not row[0].startswith("#"):
                        try:
                            m=row[headers["M"]]
                            mz=float(row[headers["MZ"]])
                            rt=float(row[headers["RT"]])
                            xn=int(row[headers["Xn"]])
                            ionMode=row[headers["Ionisation_Mode"]]
                            charge=int(row[headers["Charge"]])

                            cnMatch = "Exact"



                            hits=[]

                            if m!="" and "," not in m:
                                for hit in db.searchDB(mass=float(m), mz=None, polarity=ionMode, charges=charge, rt_min=None,
                                                       ppm=ppm,
                                                       rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                                    hits.append(hit)
                            else:
                                for hit in db.searchDB(mass=None, mz=mz, polarity=ionMode, charges=charge, rt_min=None,
                                                       ppm=ppm,
                                                       rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                                    hits.append(hit)

                            addCount=[]
                            add=[]

                            for dbName in dbs:
                                dbName=dbName.replace("DBs_", "")
                                a=["%s {Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == dbName]

                                addCount.append("%d"%len(a) if len(a)>0 else "")
                                add.append(";".join(a))

                            csvWriter.writerow(row+addCount+add)


                        except Exception as ex:
                            print ex.message
                            csvWriter.writerow(row + ["" for d in dbs]+["" for d in dbs])
                    else:
                        csvWriter.writerow(row)

                csvWriter.writerow(["###### Database search"])
                csvWriter.writerow(["## ppm=%.1f"%(ppm)])


    ppm=150.
    rtError=0.5

    db=DBSearch()

    db.addEntriesFromFile("trichoatrokontins", "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/Trichoatrokontins.tsv")
    db.optimizeDB()

    dbs=["DBs_trichoatrokontins"]


    annoFiles=[("C:/temp/results_mainExp_SFs_DBs_temp.tsv",
                "H:/180829_435_Trichoderma_TOR_SusanneZeilinger_MainExperiment/results_mainExp_SFs_DBs.tsv"),
               ("C:/temp/results_mediaComparsison_SFs_DBs_temp.tsv",
                "H:/180829_435_Trichoderma_TOR_SusanneZeilinger_MainExperiment/results_mediaComparsison_SFs_DBs.tsv"),
               ("C:/temp/results_GlutamineSideExperiment_SFs_DBs_temp.tsv",
                "H:/180829_435_Trichoderma_TOR_SusanneZeilinger_MainExperiment/results_GlutamineSideExperiment_SFs_DBs.tsv")]

    for annoFile in annoFiles:
        with open(annoFile[0], "rb") as fin:
            with open(annoFile[1], "wb") as fout:

                csvReader=csv.reader(fin, delimiter="\t")
                csvWriter=csv.writer(fout, delimiter="\t")

                headers={}
                for rowi, row in enumerate(csvReader):

                    if rowi==0:
                        for celli, cell in enumerate(row):
                            headers[cell]=celli
                        csvWriter.writerow(row+[d+"_count" for d in dbs]+[d for d in dbs])
                    elif not row[0].startswith("#"):
                        try:
                            m=row[headers["M"]]
                            mz=float(row[headers["MZ"]])
                            rt=float(row[headers["RT"]])
                            xn=int(row[headers["Xn"]])
                            ionMode=row[headers["Ionisation_Mode"]]
                            charge=int(row[headers["Charge"]])

                            cnMatch = "Exact"



                            hits=[]

                            if m!="" and "," not in m:
                                for hit in db.searchDB(mass=float(m), mz=None, polarity=ionMode, charges=charge, rt_min=None,
                                                       ppm=ppm,
                                                       rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                                    hits.append(hit)
                            else:
                                for hit in db.searchDB(mass=None, mz=mz, polarity=ionMode, charges=charge, rt_min=None,
                                                       ppm=ppm,
                                                       rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                                    hits.append(hit)

                            addCount=[]
                            add=[]

                            for dbName in dbs:
                                dbName=dbName.replace("DBs_", "")
                                a=["%s {Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == dbName]

                                addCount.append("%d"%len(a) if len(a)>0 else "")
                                add.append(";".join(a))

                            csvWriter.writerow(row+addCount+add)


                        except Exception as ex:
                            print ex.message
                            csvWriter.writerow(row + ["" for d in dbs]+["" for d in dbs])
                    else:
                        csvWriter.writerow(row)

                csvWriter.writerow(["###### Database search"])
                csvWriter.writerow(["## ppm=%.1f"%(ppm)])





if False and __name__=="__main__":
    ppm=5.
    rtError=0.5

    db=DBSearch()

    db.addEntriesFromFile("Auxin_metabolites", "D:/KangKang/180116/DBs/Auxin_metabolites.tsv")
    db.addEntriesFromFile("Auxin_TAM_conjugates_calculated", "D:/KangKang/180116/DBs/Auxin_TAM_conjugates_calculated.tsv")
    db.addEntriesFromFile("DB_Auxins_SF_RF_MZ", "D:/KangKang/180116/DBs/DB_Auxins_SF_RF_MZ.tsv")
    db.addEntriesFromFile("DBChemspider_Trp_metabolites", "D:/KangKang/180116/DBs/DBChemspider_Trp_metabolites.tsv")

    db.addEntriesFromFile("TAM_Metabolites_Apogee", "D:/KangKang/180116/DBs/TAM_Metabolites_Apogee.tsv")
    db.addEntriesFromFile("Trp_Kegg_MetExtractt", "D:/KangKang/180116/DBs/Trp_Kegg_MetExtractt.tsv")
    db.addEntriesFromFile("Trp_metabolites_399_Durum", "D:/KangKang/180116/DBs/Trp_metabolites_399_Durum.tsv")
    db.addEntriesFromFile("WheatDB", "D:/KangKang/180116/DBs/WheatDB.tsv")
    db.optimizeDB()


    with open("D:/KangKang/180116\EVAL_Trp/results.tsv", "rb") as fin:
        with open("D:/KangKang/180116\EVAL_Trp/results_DBs.tsv", "wb") as fout:

            csvReader=csv.reader(fin, delimiter="\t")
            csvWriter=csv.writer(fout, delimiter="\t")

            headers={}
            for rowi, row in enumerate(csvReader):

                if rowi==0:
                    for celli, cell in enumerate(row):
                        headers[cell]=celli
                    csvWriter.writerow(row+["Auxin_metabolites",
                                            "Auxin_TAM_conjugates_calculated",
                                            "DB_Auxins_SF_RF_MZ",
                                            "DBChemspider_Trp_metabolites",

                                            "TAM_Metabolites_Apogee",
                                            "Trp_Kegg_MetExtractt",
                                            "Trp_metabolites_399_Durum",
                                            "WheatDB",

                                            "Auxin_metabolites_count",
                                            "Auxin_TAM_conjugates_calculated_count",
                                            "DB_Auxins_SF_RF_MZ_count",
                                            "DBChemspider_Trp_metabolites_count",

                                            "TAM_Metabolites_Apogee_count",
                                            "Trp_Kegg_MetExtractt_count",
                                            "Trp_metabolites_399_Durum_count",
                                            "WheatDB_count"])
                elif not row[0].startswith("#"):
                    try:
                        m=row[headers["M"]]
                        mz=float(row[headers["MZ"]])
                        rt=float(row[headers["RT"]])
                        xn=int(row[headers["Xn"]])
                        ionMode=row[headers["Ionisation_Mode"]]
                        charge=int(row[headers["Charge"]])

                        cnMatch = "Minimum"

                        hits=[]

                        if m!="" and "," not in m:
                            for hit in db.searchDB(mass=float(m), mz=None, polarity=ionMode, charges=charge, rt_min=rt,
                                                   ppm=ppm, rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                                hits.append(hit)
                        else:
                            for hit in db.searchDB(mass=None, mz=mz, polarity=ionMode, charges=charge, rt_min=rt,
                                                   ppm=ppm, rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                                hits.append(hit)

                        a=["%s {RT: %.2f min, Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.rt_min if is_float(hit.rt_min) else -1, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == "Auxin_metabolites"]
                        b=["%s {RT: %.2f min, Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.rt_min if is_float(hit.rt_min) else -1, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == "Auxin_TAM_conjugates_calculated"]
                        c=["%s {RT: %.2f min, Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.rt_min if is_float(hit.rt_min) else -1, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == "DB_Auxins_SF_RF_MZ"]
                        d=["%s {RT: %.2f min, Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.rt_min if is_float(hit.rt_min) else -1, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == "DBChemspider_Trp_metabolites"]

                        e=["%s {RT: %.2f min, Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.rt_min if is_float(hit.rt_min) else -1, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == "TAM_Metabolites_Apogee"]
                        f=["%s {RT: %.2f min, Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.rt_min if is_float(hit.rt_min) else -1, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == "Trp_Kegg_MetExtractt"]
                        g=["%s {RT: %.2f min, Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.rt_min if is_float(hit.rt_min) else -1, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == "Trp_metabolites_399_Durum"]
                        h=["%s {RT: %.2f min, Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.rt_min if is_float(hit.rt_min) else -1, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == "WheatDB"]

                        csvWriter.writerow(row+\
                                            [str(len(a)) if len(a)>0 else ""]+ \
                                            [str(len(b)) if len(b)>0 else ""] + \
                                            [str(len(c)) if len(c)>0 else ""] + \
                                            [str(len(d)) if len(d)>0 else ""] + \

                                           [str(len(e)) if len(e) > 0 else ""] + \
                                           [str(len(f)) if len(f) > 0 else ""] + \
                                           [str(len(g)) if len(g) > 0 else ""] + \
                                           [str(len(h)) if len(h) > 0 else ""] + \

                                            [";".join(a)] + \
                                            [";".join(b)] + \
                                            [";".join(c)] + \
                                            [";".join(d)] + \

                                            [";".join(e)] + \
                                            [";".join(f)] + \
                                            [";".join(g)] + \
                                            [";".join(h)])


                    except Exception as ex:
                        import traceback
                        traceback.print_exc()
                        logging.error(str(traceback))

                        print ex.message
                        csvWriter.writerow(row + ["", "", "", "",   "", "", "", "",       "", "", "", "",   "", "", "", ""])
                else:
                    csvWriter.writerow(row)

            csvWriter.writerow(["###### Database search"])
            csvWriter.writerow(["## ppm=%.1f"%(ppm)])






if False and __name__=="__main__":
    ppm=5.
    rtError=0.5

    db=DBSearch()

    db.addEntriesFromFile("Triticumaestivum", "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/Triticum_aestivum_literature.tsv")
    db.addEntriesFromFile("Flavonoids",       "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/Flavonoids.tsv")
    db.addEntriesFromFile("PPAs",             "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/PPAs.tsv")
    db.addEntriesFromFile("Auxins",           "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/DB_Auxins_SF_RF_MZ.tsv")
    db.addEntriesFromFile("HCAs",             "C:/PyMetExtract/MetaboliteDBs/DBsForMetExtractII/HCAs_HCAAs_metabolites.tsv")




    db.optimizeDB()

    dbs=["DBs_Triticumaestivum", "DBs_Flavonoids", "DBs_PPAs", "DBs_Auxins", "DBs_HCAs"]


    annoFiles=[("H:/180713_382_Labelboxpaper/EVAL_PrimaryDataEvaluation/results_KD_KD.tsv",
                "H:/180713_382_Labelboxpaper/EVAL_PrimaryDataEvaluation/results_KD_KD_DBs.tsv"),
               ("H:/180713_382_Labelboxpaper/EVAL_PrimaryDataEvaluation/results_KD_xxx.tsv",
                "H:/180713_382_Labelboxpaper/EVAL_PrimaryDataEvaluation/results_KD_xxx_DBs.tsv"),
               ("H:/180713_382_Labelboxpaper/EVAL_PrimaryDataEvaluation/results_KM_KM.tsv",
                "H:/180713_382_Labelboxpaper/EVAL_PrimaryDataEvaluation/results_KM_KM_DBs.tsv"),
               ("H:/180713_382_Labelboxpaper/EVAL_PrimaryDataEvaluation/results_xx_KDMm.tsv",
                "H:/180713_382_Labelboxpaper/EVAL_PrimaryDataEvaluation/results_xx_KDMm_DBs.tsv")]

    for annoFile in annoFiles:
        with open(annoFile[0], "rb") as fin:
            with open(annoFile[1], "wb") as fout:

                csvReader=csv.reader(fin, delimiter="\t")
                csvWriter=csv.writer(fout, delimiter="\t")
                csvWriter=csv.writer(fout, delimiter="\t")

                headers={}
                for rowi, row in enumerate(csvReader):

                    if rowi==0:
                        for celli, cell in enumerate(row):
                            headers[cell]=celli
                        csvWriter.writerow(row+[d for d in dbs]+[d+"_extended" for d in dbs]+[d+"_count" for d in dbs]+[d for d in dbs])
                    elif not row[0].startswith("#"):
                        try:
                            m=row[headers["M"]]
                            mz=float(row[headers["MZ"]])
                            rt=float(row[headers["RT"]])
                            xn=int(row[headers["Xn"]])
                            ionMode=row[headers["Ionisation_Mode"]]
                            charge=int(row[headers["Charge"]])

                            cnMatch = "Exact"



                            hits=[]

                            if m!="" and "," not in m:
                                for hit in db.searchDB(mass=float(m), mz=None, polarity=ionMode, charges=charge, rt_min=None,
                                                       ppm=ppm,
                                                       rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                                    hits.append(hit)
                            else:
                                for hit in db.searchDB(mass=None, mz=mz, polarity=ionMode, charges=charge, rt_min=None,
                                                       ppm=ppm,
                                                       rt_error=rtError, checkXN=cnMatch, element="C", Xn=xn):
                                    hits.append(hit)

                            addCount=[]
                            add=[]
                            addextended=[]

                            for dbName in dbs:
                                dbName=dbName.replace("DBs_", "")
                                a=["%s {Type: %s, Num: %s, Sumformula: %s, AdditionalInfo: %s}" % (hit.name, hit.hitType, hit.num, hit.sumFormula, hit.additionalInfo) for hit in hits if hit.dbName == dbName]
                                b=["%s" % hit.name for hit in hits if hit.dbName == dbName]
                                addCount.append("%d"%len(a) if len(a)>0 else "")
                                addextended.append(";".join(a))
                                add.append(";".join(b))

                            csvWriter.writerow(row+add+addextended+addCount)


                        except Exception as ex:
                            print ex.message
                            csvWriter.writerow(row + ["" for d in dbs]+["" for d in dbs]+["" for d in dbs])
                    else:
                        csvWriter.writerow(row)

                csvWriter.writerow(["###### Database search"])
                csvWriter.writerow(["## ppm=%.1f"%(ppm)])

