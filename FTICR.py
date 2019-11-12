
from utils import Bunch, getRatio

from HCA import *

from MSScan import MSScan
from SGR import SGRGenerator



def matchPairs(scan, file, ppm=0.5, mzDelta=1.00335484, nativePurity=0.9893, labeledPurity=0.99, isotopologs=[], log=False):

    res=[]

    for i in range(len(msscan.mz_list)):
        mzi=msscan.mz_list[i]
        inti=msscan.intensity_list[i]

        for cn in range(3, 60):
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
                if log: print "    -->  M +1 found. Ratio is %8.3f (should be %8.3f)"%(msscan.intensity_list[ipo[0]]/inti, getRatio(nativePurity, cn, 1))
            jmoRat=0
            if jmo[0]!=-1:
                jmoRat=msscan.intensity_list[jmo[0]]/intj
                if log: print "    -->  M'-1 found. Ratio is %8.3f (should be %8.3f)"%(msscan.intensity_list[jmo[0]]/intj, getRatio(labeledPurity, cn, 1))

            #if ((ipoRat > 0 and abs(ipoRat - getRatio(nativePurity, cn, 1)) < 0.15) or ipoRat == 0) and ((jmoRat > 0 and abs(jmoRat - getRatio(labeledPurity, cn, 1)) < 0.15) or jmoRat == 0):
            if ipoRat > 0 and abs(ipoRat - getRatio(nativePurity, cn, 1)) < 0.15 and jmoRat > 0 and abs(jmoRat - getRatio(labeledPurity, cn, 1)) < 0.15:
                res.append(Bunch(sampe=file, mz=mzi, mzl=mzj, intensity=inti, intensityl=intj, cn=cn, ipoRatio=ipoRat, jmoRatio=jmoRat))


    return res


def importMSScan(file):

    msscan = MSScan()
    with open(file, "rb") as fin:

        for rowi, row in enumerate(fin):
            if rowi == 0:
                pass
            else:
                a = row.split("\t")

                msscan.mz_list.append(float(a[0]))
                msscan.intensity_list.append(float(a[1]))

        print "File", file, "read.", len(msscan.mz_list), "signals in the file.."
        print " ..",

    return msscan


def annotateWithSFs(allRes, ppm=2):
    foundSFs = []
    sgr = SGRGenerator()
    for cn in set([r.cn for r in allRes]):
        # print "\nCn:", cn

        hc = HierarchicalClustering(sorted([r for r in allRes if r.cn == cn], key=lambda x: x.mz),
                                    dist=lambda x, y: x.getValue() - y.getValue(),
                                    val=lambda x: x.mz, mean=lambda x, y: x / y, add=lambda x, y: x + y)

        t = cutTreeSized(hc.getTree(), ppm=2)
        # printTree(hc.getTree(), inlet="    ")
        # print " -----"
        for c in t:
            mzs = [r.getObject().mz for r in c.getKids()]
            # print "  %5.3f ppm error (%8.5f-%8.5f): "%(abs(min(mzs)-max(mzs))*1000000./max(mzs), min(mzs), max(mzs)), mzs

            meanmz = sum(mzs) / len(mzs)
            sfs = sgr.findFormulas(meanmz + 1.007276, useAtoms=["C", "N", "H", "O", "P", "S"],
                                   atomsRange=[(cn, cn), (0, 500), (0, 10000), (0, 400), [0, 10], [0, 20]],
                                   useSevenGoldenRules=True, useSecondRule=True, ppm=ppm)

            sfs.extend(sgr.findFormulas(meanmz - 36.969402, useAtoms=["C", "N", "H", "O", "P", "S"],
                                        atomsRange=[(cn, cn), (0, 500), (0, 10000), (0, 400), [0, 10], [0, 20]],
                                        useSevenGoldenRules=True, useSecondRule=True, ppm=ppm))

            # for sf in sfs:
            #    print "      --> ", sf

            if len(sfs) > 0:
                foundSFs.append(Bunch(sfs=sfs, meanMZ=meanmz, cn=cn, files={}))

    return foundSFs


if __name__=="__main__":

    allRes=[]
    msScans={}

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


    ppm=0.5
    for file in files:

        msscan=importMSScan(file)
        msScans[file]=msscan

        res=matchPairs(msscan, file, ppm=ppm)
        allRes.extend(res)

        print "\r   --> %d pairs found"%(len(res))

    print "\n\nHCA\n\n"
    foundSFs=annotateWithSFs(allRes, ppm=ppm)


    print "%d signal pairs were annotated with unique sum formulas. %d non-unique sum formulas"%(len(set(foundSFs)), len([f for f in foundSFs if len(f.sfs)>1]))


    for file in files:
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


    with open("H:/190318_538_AlternariaII/FT Daten/results_MEII.txt", "wb") as fout:

        fout.write("\t".join(["Num", "SF", "Cn"]+[a[a.rfind("/")+1:].replace(".tsv", "") for a in files]+[a[a.rfind("/")+1:].replace(".tsv", "")+"_L" for a in files]))
        fout.write("\n")

        ind=1
        for f in foundSFs:
            fout.write("\t".join(["FTICR_%d"%ind, ";".join(f.sfs), "%d"%f.cn]))
            fout.write("\t")
            fout.write("\t".join([str(f.files[i]) for i in files]))
            fout.write("\t")
            fout.write("\t".join([str(f.files[i+"_L"]) for i in files]))
            fout.write("\n")
            ind+=1

        fout.write("## Data processing with new MEII-Spectra modul")
        fout.write("\n")
        fout.write("## ppm: %.1f (for matching of mz signals and sumformula generation)"%ppm)
        fout.write("\n")





