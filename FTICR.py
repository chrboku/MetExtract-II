
from utils import Bunch, getRatio
from formulaTools import formulaTools

from HCA import *

from MSScan import MSScan
from SGR import SGRGenerator

from copy import deepcopy







def importMSScan(file, colNameMZ="MZ", colNameIntensity="Int", log=False):

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

        if log: print "File", file, "read.", len(msscan.mz_list), "signals in the file.."
        if log: print " ..",

    return msscan















def matchPairs(msscan, file, ppm=0.5, mzDelta=1.00335484, nativePurity=0.9893, labeledPurity=0.99, cns=[3,60], isotopologsCanBeOmitted=False, log=False):

    res=[]

    for i in range(len(msscan.mz_list)):
        mzi=msscan.mz_list[i]
        inti=msscan.intensity_list[i]

        for cn in range(cns[0], cns[1]):
            mzj=mzi+cn*mzDelta
            bounds=msscan.findMZ(mzj, ppm)

            if bounds[0]==-1:
                continue

            if bounds[0]!=bounds[1]:
                print "found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (%.5f) in the range with the highest abundance is used"%(mzj)
                bounds=msscan.getMostIntensePeak(bounds[0], bounds[1])

            j=bounds[0]
            intj=msscan.intensity_list[j]
            if log: print "Peak pair found: Cn: %2d,    MZs %9.5f / %9.5f,    Intensities: %8.1f / %8.1f  (%7.3f),     "%(cn, mzi, mzj, inti, intj, inti/intj)




            ipo=msscan.findMZ(mzi+mzDelta, ppm)
            jmo=msscan.findMZ(mzj-mzDelta, ppm)

            if ipo[0]!=ipo[1]:
                print "    --> found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (M+1) in the range with the highest abundance is used"
                ipo=msscan.getMostIntensePeak(ipo[0], ipo[1])
            if jmo[0]!=jmo[1]:
                print "    --> found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (M'-1) in the range with the highest abundance is used"
                jmo=msscan.getMostIntensePeak(jmo[0], jmo[1])

            ipoRat=0
            if ipo[0]!=-1:
                ipoRat=msscan.intensity_list[ipo[0]]/inti
                if log: print "    -->  M +1 found. Ratio is %8.3f (should be %8.3f)"%(ipoRat, getRatio(nativePurity, cn, 1))
            jmoRat=0
            if jmo[0]!=-1:
                jmoRat=msscan.intensity_list[jmo[0]]/intj
                if log: print "    -->  M'-1 found. Ratio is %8.3f (should be %8.3f)"%(jmo, getRatio(labeledPurity, cn, 1))




            imo=msscan.findMZ(mzi-mzDelta, ppm)
            jpo=msscan.findMZ(mzj+mzDelta, ppm)

            if imo[0]!=imo[1]:
                print "    --> found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (M+1) in the range with the highest abundance is used"
                imo=msscan.getMostIntensePeak(imo[0], imo[1])
            if jpo[0]!=jpo[1]:
                print "    --> found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (M'-1) in the range with the highest abundance is used"
                jpo=msscan.getMostIntensePeak(jpo[0], jpo[1])

            imoRat=0
            if imo[0]!=-1:
                imoRat=msscan.intensity_list[imo[0]]/inti
                if log: print "    -->  M -1 found. Ratio is %8.3f"%(imoRat)
            jpoRat=0
            if jpo[0]!=-1:
                jpoRat=msscan.intensity_list[jpo[0]]/intj
                if log: print "    -->  M'+1 found. Ratio is %8.3f"%(jpoRat)




            if ( (    isotopologsCanBeOmitted and ((ipoRat > 0 and abs(ipoRat - getRatio(nativePurity, cn, 1)) < 0.15) or ipoRat == 0) and ((jmoRat > 0 and abs(jmoRat - getRatio(labeledPurity, cn, 1)) < 0.15) or jmoRat == 0)) or \
                 (not isotopologsCanBeOmitted and ipoRat > 0 and abs(ipoRat - getRatio(nativePurity, cn, 1)) < 0.15 and jmoRat > 0 and abs(jmoRat - getRatio(labeledPurity, cn, 1)) < 0.15)  ) and \
                    jpoRat < 0.05 and imoRat < 0.05:
                res.append(Bunch(sampe=file, mz=mzi, mzl=mzj, intensity=inti, intensityl=intj, cn=cn, ipoRatio=ipoRat, jmoRatio=jmoRat))


    return res



def annotateWithSFs(allRes, clusteringPPM=2, annotationPPM=0.5,
                    atoms=["C", "N", "H", "O", "P", "S"], atomsRange=[(0, 0), (0, 0), (0, 0), (0, 0), [0, 0], [0, 0]],
                    useSevenGoldenRules=True):
    foundSFs = []
    sgr = SGRGenerator()
    fT=formulaTools()
    for cn in set([r.cn for r in allRes]):

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

            a = sgr.findFormulas(meanmz + 1.007276, useAtoms=deepcopy(atoms),
                                 atomsRange=deepcopy(atomsRange),
                                 useSevenGoldenRules=useSevenGoldenRules, useSecondRule=True, ppm=annotationPPM)
            for t in a:
                mz=fT.calcMolWeight(fT.parseFormula(t))-1.007276
                ppmError=(meanmz-mz)*1E6/meanmz
                sfs.append(Bunch(sf=t, ion="[M-H]-", theoMZ=mz, ppmError=ppmError))

            a = sgr.findFormulas(meanmz - 36.969402, useAtoms=deepcopy(atoms),
                                 atomsRange=deepcopy(atomsRange),
                                 useSevenGoldenRules=useSevenGoldenRules, useSecondRule=True, ppm=annotationPPM)
            for t in a:
                mz=fT.calcMolWeight(fT.parseFormula(t))+36.969402
                ppmError=(meanmz-mz)*1E6/meanmz
                sfs.append(Bunch(sf=t, ion="[M+Cl]-", theoMZ=mz, ppmError=ppmError))

            foundSFs.append(Bunch(sfs=sfs, meanMZ=meanmz, cn=cn, files={}))

    return foundSFs




def processMSScans(files, ppm=0.5, clusteringPPM=2, annotationPPM=0.5,
                   atoms=["C", "N", "H", "O", "P", "S"], atomsRange=[(0, 0), (0, 0), (0, 0), (0, 0), [0, 0], [0, 0]],
                   useSevenGoldenRules=True,
                   setFunctionMax=None, setFunctionValue=None, setFunctionText=None):
    allRes=[]
    if setFunctionMax!=None: setFunctionMax(len(files)+1)
    if setFunctionValue!=None: setFunctionValue(0)
    for i, fileName in enumerate(files.keys()):
        if setFunctionText!=None: setFunctionText("Processing file %s"%fileName)
        msScan=files[fileName]

        res=matchPairs(msScan, fileName, ppm=ppm)
        print "%d signals pairs found in %s"%(len(res), fileName)
        if setFunctionValue!=None: setFunctionValue(i)
        allRes.extend(res)

    if setFunctionText != None: setFunctionText("Combining results and generating sum formulas")
    s = "".join(["%s(%d-%d)" % (atoms[i], atomsRange[i][0], atomsRange[i][1]) for i in range(len(atoms))])
    print "Annotating elements with %s" % s
    foundSFs=annotateWithSFs(allRes, clusteringPPM=clusteringPPM, annotationPPM=annotationPPM,
                             atoms=deepcopy(atoms), atomsRange=deepcopy(atomsRange),
                             useSevenGoldenRules=useSevenGoldenRules)
    if setFunctionValue!=None: setFunctionValue(len(files)+1)


    return foundSFs



def generateDataMatrix(msScans, foundSFs, ppm=0.5):
    for file in msScans.keys():
        msscan=msScans[file]
        for f in foundSFs:
            b=msscan.findMZ(f.meanMZ, ppm=ppm)
            if b[0]==-1:
                f.files[file]=""
            else:
                if b[0] != b[1]:
                    print "found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (%.5f) in the range with the highest abundance is used" % (f.meanMZ)
                    b = msscan.getMostIntensePeak(b[0], b[1])

                f.files[file]=msscan.intensity_list[b[0]]

            b = msscan.findMZ(f.meanMZ+f.cn*1.00335484, ppm=ppm)
            if b[0] == -1:
                f.files[file+"_L"] = ""
            else:
                if b[0] != b[1]:
                    print "found several peaks. As this usually never happens with FT-ICR data, it is not considered here yet. The signal (%.5f) in the range with the highest abundance is used" % (f.meanMZ+f.cn*1.00335484)
                    b = msscan.getMostIntensePeak(b[0], b[1])

                f.files[file+"_L"] = msscan.intensity_list[b[0]]

    return foundSFs


def writeMatrixToFile(outFile, foundSFs, msScans, ppm):
    with open(outFile, "wb") as fout:

        fout.write("\t".join(["Num", "MeanMZ", "Cn", "SF_all", "SF"]+[a[a.rfind("/")+1:].replace(".tsv", "") for a in msScans.keys()]+[a[a.rfind("/")+1:].replace(".tsv", "")+"_L" for a in msScans.keys()]))
        fout.write("\n")

        ind=1
        for f in foundSFs:
            fout.write("\t".join(["FTICR_%d"%ind, str(f.meanMZ), "%d"%f.cn, ";".join([str(sf) for sf in f.sfs]), "" if len(f.sfs)!=1 else f.sfs[0].sf]))
            fout.write("\t")
            fout.write("\t".join([str(f.files[i]) for i in msScans.keys()]))
            fout.write("\t")
            fout.write("\t".join([str(f.files[i+"_L"]) for i in msScans.keys()]))
            fout.write("\n")
            ind+=1

        fout.write("## Data processing with new MEII-Spectra modul")
        fout.write("\n")
        fout.write("## ppm: %.1f (for matching of mz signals and sumformula generation)"%ppm)
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

    msScans={}
    for file in files:
        msscan=importMSScan(file)
        msScans[file]=msscan
        print "Read file '%s' with %d signals"%(file, len(msScans[file].mz_list))

    allRes=[]
    ppm=0.5
    clusteringPPM=2
    foundSFs=processMSScans(msScans, ppm=ppm, clusteringPPM=clusteringPPM)
    print "%d signal pairs were annotated with unique sum formulas. %d non-unique sum formulas"%(len(set(foundSFs)), len([f for f in foundSFs if len(f.sfs)>1]))

    generateDataMatrix(msScans, foundSFs, ppm=ppm)
    writeMatrixToFile("H:/190318_538_AlternariaII/FT Daten/results_MEII.txt", foundSFs, msScans, ppm=ppm)





