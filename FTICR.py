
from utils import Bunch, getRatio
from formulaTools import formulaTools

from Chromatogram import Chromatogram

from HCA import *

from MSScan import MSScan
from SGR import SGRGenerator

from copy import deepcopy





class FTICRProcessing:

    def __init__(self, log=False):
        self._log=log
        self.history=[]

        self.fT=formulaTools()

    def addToHistory(self, mes):
        self.history.append(mes)

    def log(self, mes):
        if self._log:
            print mes

    def importMSScan(self, file, colNameMZ="MZ", colNameIntensity="Int"):

        if file.lower().endswith(".mzxml"):
            chrom=Chromatogram()
            chrom.parse_file(file)

            if len(chrom.MS1_list)!=1:
                raise Exception("There are more than 1 MSscans in the file %s"%file)

            return chrom.MS1_list[0]

        elif file.lower().endswith(".tsv"):
            msscan = MSScan()
            with open(file, "rb") as fin:

                headers={}
                for rowi, row in enumerate(fin):
                    if rowi == 0:
                        cells=row.strip().split("\t")
                        for celli, cell in enumerate(cells):
                            headers[cell]=celli
                    else:
                        cells = row.split("\t")
                        msscan.mz_list.append(float(cells[headers[colNameMZ]]))
                        msscan.intensity_list.append(float(cells[headers[colNameIntensity]]))

                self.log("File %s read. %dsignals in the file.."%(file, len(msscan.mz_list)))

            return msscan


        raise Exception("Unknown file type")

    def matchPairs(self, msscan, file, ppm=0.5, mzDelta=1.00335484, nativePurity=0.9893, labeledPurity=0.99, maxEnrDeviation=.015, cns=[3,60], isotopologsCanBeOmitted=False, intensityThreshold=0):

        res=[]

        for i in range(len(msscan.mz_list)):
            mzi=msscan.mz_list[i]
            inti=msscan.intensity_list[i]

            if inti>=intensityThreshold:

                for cn in range(cns[0], cns[1]):
                    mzj=mzi+cn*mzDelta
                    bounds=msscan.findMZ(mzj, ppm)

                    if bounds[0]==-1:
                        continue

                    self.log("found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (%.5f) in the range with the highest abundance is used"%(mzj))
                    j=msscan.getMostIntensePeak(bounds[0], bounds[1])

                    intj=msscan.intensity_list[j]

                    if intj>=intensityThreshold:
                        self.log("Peak pair found: Cn: %2d,    MZs %9.5f / %9.5f,    Intensities: %8.1f / %8.1f  (%7.3f),     "%(cn, mzi, mzj, inti, intj, inti/intj))

                        ipo=msscan.findMZ(mzi+mzDelta, ppm)
                        jmo=msscan.findMZ(mzj-mzDelta, ppm)

                        self.log("    --> found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (M+1) in the range with the highest abundance is used")
                        ipo=msscan.getMostIntensePeak(ipo[0], ipo[1])
                        self.log("    --> found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (M'-1) in the range with the highest abundance is used")
                        jmo=msscan.getMostIntensePeak(jmo[0], jmo[1])

                        ipoRat=0
                        if ipo!=-1:
                            ipoRat=msscan.intensity_list[ipo]/inti
                            self.log("    -->  M +1 found. Ratio is %8.3f (should be %8.3f)"%(ipoRat, getRatio(nativePurity, cn, 1)))
                        jmoRat=0
                        if jmo!=-1:
                            jmoRat=msscan.intensity_list[jmo]/intj
                            self.log("    -->  M'-1 found. Ratio is %8.3f (should be %8.3f)"%(jmo, getRatio(labeledPurity, cn, 1)))

                        imo=msscan.findMZ(mzi-mzDelta, ppm)
                        jpo=msscan.findMZ(mzj+mzDelta, ppm)

                        self.log("    --> found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (M+1) in the range with the highest abundance is used")
                        imo=msscan.getMostIntensePeak(imo[0], imo[1])
                        self.log("    --> found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (M'-1) in the range with the highest abundance is used")
                        jpo=msscan.getMostIntensePeak(jpo[0], jpo[1])

                        imoRat=0
                        if imo!=-1:
                            imoRat=msscan.intensity_list[imo]/inti
                            self.log("    -->  M -1 found. Ratio is %8.3f"%(imoRat))
                        jpoRat=0
                        if jpo!=-1:
                            jpoRat=msscan.intensity_list[jpo]/intj
                            self.log("    -->  M'+1 found. Ratio is %8.3f"%(jpoRat))

                        if ( (    isotopologsCanBeOmitted and ((ipoRat > 0 and abs(ipoRat - getRatio(nativePurity, cn, 1)) < maxEnrDeviation) or ipoRat == 0) and ((jmoRat > 0 and abs(jmoRat - getRatio(labeledPurity, cn, 1)) < maxEnrDeviation) or jmoRat == 0)) or \
                             (not isotopologsCanBeOmitted and ipoRat > 0 and abs(ipoRat - getRatio(nativePurity, cn, 1)) < maxEnrDeviation and jmoRat > 0 and abs(jmoRat - getRatio(labeledPurity, cn, 1)) < maxEnrDeviation)  ) and \
                                jpoRat < 0.05 and imoRat < 0.05:
                            res.append(Bunch(sampe=file, mz=mzi, mzl=mzj, intensity=inti, intensityl=intj, cn=cn, ipoRatio=ipoRat, jmoRatio=jmoRat))


        return res



    def annotateWithSFs(self, allRes, msScans, matchPPM=0.5, clusteringPPM=2, annotationPPM=0.5,
                        nativePurity=0.9893, labeledPurity=0.99,
                        atoms=["C", "N", "H", "O", "P", "S"], atomsRange=[(0, 0), (0, 0), (0, 0), (0, 0), [0, 0], [0, 0]],
                        adducts={"[M-H]-": -1.007276, "[M+Cl]-": 36.969402},
                        useSevenGoldenRules=True):
        foundSFs = []
        sgr = SGRGenerator()
        for cn in set([r.cn for r in allRes]):

            theoRatNative=getRatio(nativePurity, cn, 1)
            theoRatLabeled=getRatio(labeledPurity, cn, 1)

            hc = HierarchicalClustering(sorted([r for r in allRes if r.cn == cn], key=lambda x: x.mz),
                                        dist=lambda x, y: x.getValue() - y.getValue(),
                                        val=lambda x: x.mz, mean=lambda x, y: x / y, add=lambda x, y: x + y)

            t = cutTreeSized(hc.getTree(), ppm=clusteringPPM)
            for c in t:
                mzs = [r.getObject().mz for r in c.getKids()]

                sfs=[]
                meanmz = sum(mzs) / len(mzs)

                atoms=deepcopy(atoms)
                atomsRange=deepcopy(atomsRange)

                for elemi, elem in enumerate(atoms):
                    if elem=="C":
                        atomsRange[elemi]=(cn, cn)

                for adductName, offset in adducts.items():

                    a = sgr.findFormulas(meanmz-offset,
                                         useAtoms=deepcopy(atoms), atomsRange=deepcopy(atomsRange),
                                         useSevenGoldenRules=useSevenGoldenRules, useSecondRule=True, ppm=annotationPPM)
                    for t in a:
                        mzTheo=self.fT.calcMolWeight(self.fT.parseFormula(t))+offset
                        ppmError=(meanmz-mzTheo)*1E6/meanmz
                        sfs.append(Bunch(sf=t, ion=adductName, theoMZ=mzTheo, ppmError=ppmError))

                ipos={"native":{},"labeled":{}, "mzError":{}}

                for file in msScans.keys():
                    msScan=msScans[file]
                    b  = msScan.findMZ(meanmz, ppm=matchPPM)
                    b1 = msScan.findMZ(meanmz+1.00335484, ppm=matchPPM)
                    if b[0]!=-1 and b1[0]!=-1:
                        b  = msScan.getMostIntensePeak(b[0], b[1])
                        b1 = msScan.getMostIntensePeak(b1[0], b1[1])
                        ipos["native"][file] = msScan.intensity_list[b1] / msScan.intensity_list[b] - theoRatNative

                        if len(sfs)==1:
                            ipos["mzError"][file]=(msScan.mz_list[b]-sfs[0].theoMZ)*1E6/sfs[0].theoMZ

                    b = msScan.findMZ(meanmz+cn*1.00335484, ppm=matchPPM)
                    b1 = msScan.findMZ(meanmz+cn*1.00335484 - 1.00335484, ppm=matchPPM)
                    if b[0] != -1 and b1[0] != -1:
                        b = msScan.getMostIntensePeak(b[0], b[1])
                        b1 = msScan.getMostIntensePeak(b1[0], b1[1])
                        ipos["labeled"][file] = msScan.intensity_list[b1] / msScan.intensity_list[b] - theoRatLabeled


                foundSFs.append(Bunch(sfs=sfs, dbs=[], meanMZ=meanmz, cn=cn, files={}, isotopologRatios=ipos))

        return foundSFs


    def annotateWithDBs(self, foundSFs, annotationPPM=0.5, dbs=[], adducts={"[M-H]-": -1.007276, "[M+Cl]-": 36.969402}):
        for db in dbs:
            ## Bunch(DBName=, sumFormula=, name=)
            db.m=self.fT.calcMolWeight(self.fT.parseFormula(db.sumFormula))
            db.cn=self.fT.parseFormula(db.sumFormula)["C"]

        for foundSF in foundSFs:
            for adductName, offset in adducts.items():
                m=foundSF.meanMZ-offset

                for db in dbs:
                    ppmError=abs(m-db.m)*1E6/m
                    if ppmError<annotationPPM and foundSF.cn==db.cn:
                        foundSF.dbs.append(Bunch(ion=adductName, dbName=db.DBName, name=db.name, ppmError=ppmError))

        return foundSFs



    def applyAverageMZErrorPerFile(self, foundSFs, msScans):
        avMZErrorPerFile={}

        for file in msScans:
            allMZErrors=[]
            for foundSF in foundSFs:
                if file in foundSF.isotopologRatios["mzError"]:
                    allMZErrors.append(foundSF.isotopologRatios["mzError"][file])

            if len(allMZErrors)>0:
                avMZError=sum(allMZErrors)/len(allMZErrors)
                avMZErrorPerFile[file]=avMZError
                self.log("  -- Average mz error for file %s: %.3f"%(file, avMZError))
                self.addToHistory("  -- Average mz error for file %s: %.3f"%(file, avMZError))

                for mzi, mz in enumerate(msScans[file].mz_list):
                    msScans[file].mz_list[mzi]=mz-mz*avMZError/1E6

        return avMZErrorPerFile




    def processMSScans(self, msScans,
                       enr12C=0.9893, enr13C=0.99, maxEnrDeviation=0.15,
                       ppm=0.5, clusteringPPM=2, annotationPPM=0.5,
                       intensityThreshold=0,
                       atoms=["C", "N", "H", "O", "P", "S"], atomsRange=[(0, 0), (0, 0), (0, 0), (0, 0), [0, 0], [0, 0]],
                       adducts={"[M-H]-": -1.007276, "[M+Cl]-": 36.969402},
                       recalibrateWithUniqueSumFormulas=True,
                       useSevenGoldenRules=True,
                       dbs=[],
                       setFunctionMax=None, setFunctionValue=None, setFunctionText=None):

        self.addToHistory("Data processing with the new FT-ICR modul of MetExtract II")
        self.addToHistory("----------------------------------------------------------")
        self.addToHistory("Parameters")
        self.addToHistory("Enrichment native metabolite forms: %.3f%%"%(enr12C*100))
        self.addToHistory("Enrichment 13C-labeled metabolite forms: %.3f%%"%(enr13C*100))
        self.addToHistory("Minimum peak intensity (for M and M'): %d"%intensityThreshold)
        self.addToHistory("Maximum enrichment deviation: +/- %.3f"%(maxEnrDeviation))
        self.addToHistory("Search for signal pairs maximum allowed error: %.1f ppm"%ppm)
        self.addToHistory("Clustering of found signals pairs among samples maximum allowed cluster size: %.1f ppm"%clusteringPPM)
        self.addToHistory("Annotation with sum formulas and databases maximum mz deviation: %.1f ppm"%annotationPPM)
        self.addToHistory("Sum formula generation: %s, use Seven Golden Rules (Kind et al. 2007) %s"%("".join(["%s(%d-%d)" % (atoms[i], atomsRange[i][0], atomsRange[i][1]) for i in range(len(atoms))]), useSevenGoldenRules))
        self.addToHistory("Note: The number of carbon atoms detected for the signal pairs will be used")
        self.addToHistory("Used adducts for annotating signal pair ions: %s"%str(adducts))
        self.addToHistory("Calibrate MS spectra after sum formula generation: %s"%recalibrateWithUniqueSumFormulas)
        self.addToHistory("Using %d database hits for signal pair annotation"%(len(dbs)))

        self.addToHistory("----------------------------------------------------------")
        allRes=[]
        if setFunctionMax!=None: setFunctionMax(len(msScans)*2 + 4)
        if setFunctionValue!=None: setFunctionValue(0)
        for i, fileName in enumerate(msScans.keys()):
            if setFunctionText!=None: setFunctionText("Processing file %s"%fileName)
            msScan=msScans[fileName]

            res=self.matchPairs(msScan, fileName, ppm=ppm, nativePurity=enr12C, labeledPurity=enr13C, maxEnrDeviation=maxEnrDeviation, intensityThreshold=intensityThreshold)
            self.log("%d signals pairs found in %s"%(len(res), fileName))
            self.addToHistory("%d signals pairs found in %s"%(len(res), fileName))
            if setFunctionValue!=None: setFunctionValue(i)
            allRes.extend(res)

        if setFunctionText != None: setFunctionText("Combining results and generating sum formulas")
        s = "".join(["%s(%d-%d)" % (atoms[i], atomsRange[i][0], atomsRange[i][1]) for i in range(len(atoms))])
        self.log("Annotating elements with %s" % s)
        foundSFs=self.annotateWithSFs(allRes, msScans, matchPPM=ppm, clusteringPPM=clusteringPPM, annotationPPM=annotationPPM,
                                 nativePurity=enr12C, labeledPurity=enr13C,
                                 atoms=deepcopy(atoms), atomsRange=deepcopy(atomsRange),
                                 adducts=adducts,
                                 useSevenGoldenRules=useSevenGoldenRules)

        averageMZErrors={}
        if recalibrateWithUniqueSumFormulas:
            if setFunctionValue != None: setFunctionValue(len(msScans) + 1)
            if setFunctionText != None: setFunctionText("Re-calibrating files")
            self.log("Applying calibration based on generated unique sum formulas..")
            self.addToHistory("Spectra will be calibrated based")
            averageMZErrors=self.applyAverageMZErrorPerFile(foundSFs, msScans)

            self.log("re-doing analysis with calibrated MS scans")

            allRes=[]
            if setFunctionValue!=None: setFunctionValue(len(msScans)+2)
            for i, fileName in enumerate(msScans.keys()):
                if setFunctionText!=None: setFunctionText("Processing file %s"%fileName)
                msScan=msScans[fileName]

                res=self.matchPairs(msScan, fileName, ppm=ppm, nativePurity=enr12C, labeledPurity=enr13C, maxEnrDeviation=maxEnrDeviation, intensityThreshold=intensityThreshold)
                self.log("%d signals pairs found in %s"%(len(res), fileName))
                self.addToHistory("%d signals pairs found in %s after calibration" % (len(res), fileName))
                if setFunctionValue!=None: setFunctionValue(len(msScans)+2+i)
                allRes.extend(res)

            if setFunctionText != None: setFunctionText("Combining results and generating sum formulas")
            s = "".join(["%s(%d-%d)" % (atoms[i], atomsRange[i][0], atomsRange[i][1]) for i in range(len(atoms))])
            self.log("Annotating elements with %s" % s)
            foundSFs=self.annotateWithSFs(allRes, msScans, matchPPM=ppm, clusteringPPM=clusteringPPM, annotationPPM=annotationPPM,
                                     nativePurity=enr12C, labeledPurity=enr13C,
                                     atoms=deepcopy(atoms), atomsRange=deepcopy(atomsRange),
                                     adducts=adducts,
                                     useSevenGoldenRules=useSevenGoldenRules)


        foundSFs=self.annotateWithDBs(foundSFs, annotationPPM=annotationPPM, dbs=dbs)
        if setFunctionValue!=None: setFunctionValue(len(msScans) + 1)

        return foundSFs, averageMZErrors



    def generateDataMatrix(self, msScans, foundSFs, ppm=0.5):
        for file in msScans.keys():
            msscan=msScans[file]
            for f in foundSFs:
                b=msscan.findMZ(f.meanMZ, ppm=ppm)
                if b[0]==-1:
                    f.files[file]=""
                else:
                    self.log("found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (%.5f) in the range with the highest abundance is used" % (f.meanMZ))
                    b = msscan.getMostIntensePeak(b[0], b[1])

                    f.files[file]=msscan.intensity_list[b]

                b = msscan.findMZ(f.meanMZ+f.cn*1.00335484, ppm=ppm)
                if b[0] == -1:
                    f.files[file+"_L"] = ""
                else:
                    self.log("found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (%.5f) in the range with the highest abundance is used" % (f.meanMZ+f.cn*1.00335484))
                    b = msscan.getMostIntensePeak(b[0], b[1])

                    f.files[file+"_L"] = msscan.intensity_list[b]

        return foundSFs


    def writeMatrixToFile(self, outFile, foundSFs, msScans, ppm):
        with open(outFile, "wb") as fout:

            filesOrder=sorted(msScans.keys())
            fout.write("\t".join(["Num", "MeanMZ", "Cn", "SF_all", "SF", "SFErrorPPM", "DB_all", "DB"]+[a[a.rfind("/")+1:].replace(".tsv", "") for a in filesOrder]+[a[a.rfind("/")+1:].replace(".tsv", "")+"_L" for a in filesOrder]))
            fout.write("\n")

            ind=1
            for f in foundSFs:
                fout.write("\t".join(["FTICR_%d"%ind, str(f.meanMZ), "%d"%f.cn, "; ".join([str(sf) for sf in f.sfs]), "" if len(f.sfs)!=1 else f.sfs[0].sf, "" if len(f.sfs)!=1 else "%.3f"%((f.meanMZ-f.sfs[0].theoMZ)*1E6/f.sfs[0].theoMZ), ";".join([str(db) for db in f.dbs]), "" if len(f.dbs)!=1 else f.dbs[0].name]))
                fout.write("\t")
                fout.write("\t".join([str(f.files[i]) for i in filesOrder]))
                fout.write("\t")
                fout.write("\t".join([str(f.files[i+"_L"]) for i in filesOrder]))
                fout.write("\n")
                ind+=1


            for hist in self.history:
                fout.write("## %s"%hist)
                fout.write("\n")



if __name__=="__main__":


    files=[
        "H:/190318_538_AlternariaII/FT Daten/STAA1_37.tsv",
        "H:/190318_538_AlternariaII/FT Daten/STAA2_35.tsv",
        "H:/190318_538_AlternariaII/FT Daten/STAA3_39.tsv",
        "H:/190318_538_AlternariaII/FT Daten/STAA4_40.tsv",
        "H:/190318_538_AlternariaII/FT Daten/STAA5_43.tsv",

        "H:/190318_538_AlternariaII/FT Daten/STAS1_38.tsv",
        "H:/190318_538_AlternariaII/FT Daten/STAS2_41.tsv",
        "H:/190318_538_AlternariaII/FT Daten/STAS3_42.tsv",
        "H:/190318_538_AlternariaII/FT Daten/STAS4_36.tsv",
        "H:/190318_538_AlternariaII/FT Daten/STAS5_44.tsv",

        "H:/190318_538_AlternariaII/FT Daten/MDAA_46.tsv",
        "H:/190318_538_AlternariaII/FT Daten/MDAS_47.tsv",

        "H:/190318_538_AlternariaII/FT Daten/MDblanks_45.tsv"
    ]

    fticr=FTICRProcessing()

    msScans={}
    for file in files:
        msscan=fticr.importMSScan(file)
        msScans[file]=msscan
        print "Read file '%s' with %d signals"%(file, len(msScans[file].mz_list))

    allRes=[]
    ppm=0.5
    clusteringPPM=2
    foundSFs=fticr.processMSScans(msScans, ppm=ppm, clusteringPPM=clusteringPPM)
    print "%d signal pairs were annotated with unique sum formulas. %d non-unique sum formulas"%(len(set(foundSFs)), len([f for f in foundSFs if len(f.sfs)>1]))

    fticr.generateDataMatrix(msScans, foundSFs, ppm=ppm)
    fticr.writeMatrixToFile("H:/190318_538_AlternariaII/FT Daten/results_MEII.txt", foundSFs, msScans, ppm=ppm)





