import logging
import LoggingSetup
LoggingSetup.LoggingSetup.Instance().initLogging()

import matplotlib
import matplotlib.pyplot as plt

from Chromatogram import Chromatogram
from utils import Bunch, merge_dicts


from utils import createTableFromBunch, SQLInsert, writeObjectAsSQLInsert, SQLSelectAsObject
import sqlite3
import os
import os.path
import pickle
import base64
from sumFormula import sumFormulaGenerator
from formulaTools import formulaTools

from MSScan import MS2Scan

from copy import deepcopy


from math import floor

from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.charts.lineplots import ScatterPlot
from reportlab.lib import pagesizes
from reportlab.lib import colors
from reportlab.lib.colors import Color
from reportlab.lib.units import mm
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.platypus.flowables import Flowable
from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF
from reportlab.lib.styles import getSampleStyleSheet







def cropEICs(eic, times, scanIDs, startRT, stopRT):
    eicNew=[]
    timesNew=[]
    scanIDsNew=[]
    for i in range(len(times)):
        if startRT<=times[i]<=stopRT:
            eicNew.append(eic[i])
            timesNew.append(times[i])
            scanIDsNew.append(scanIDs[i])
    return eicNew, timesNew, scanIDsNew






# helper function for ReportLab
def coord(x, y):
    return x, y

# helper function for ReportLab
def noLabel(a):
    return ""

# helper class for ReportLab
class TTR(Flowable):  #TableTextRotate
    '''Rotates a tex in a table cell.'''

    def __init__(self, text):
        Flowable.__init__(self)
        self.text = text

    def draw(self):
        canvas = self.canv
        canvas.rotate(90)
        canvas.drawString(0, -1, self.text)





from collections import defaultdict as ddict
from math import sqrt
def similarity(spec1, spec2):

    vector1 = ddict(int)
    vector2 = ddict(int)
    mzs = set()
    for mz, i in spec1:
        vector1[round(mz,1)] += i
        mzs.add(round(mz,1))
    for mz, i in spec2:
        vector2[round(mz,1)] += i
        mzs.add(round(mz,1))

    z = 0
    n_v1 = 0
    n_v2 = 0

    for mz in mzs:
        int1 = vector1[mz]
        int2 = vector2[mz]
        z += int1*int2
        n_v1 += int1*int1
        n_v2 += int2*int2
    try:
        cosine = z / (sqrt(n_v1) * sqrt(n_v2))
    except:
        cosine = 0.0
    return cosine


class AnnotatedMSMSSpectra:
    def __init__(self, scanMS2Native, scanMS2Labelled, charge=1, ionMode='+'):
        self.nativeMz_list=scanMS2Native.mz_list
        self.nativeIntensity_list=scanMS2Native.intensity_list

        self.labelledMz_list=scanMS2Labelled.mz_list
        self.labelledIntensity_list=scanMS2Labelled.intensity_list

        self.peakAnnotations=[]

        self.charge=charge
        self.ionMode=ionMode

    def getAnnotationsForPeakInScan(self, native=False, labelled=False, index=-1):
        getInd=None

        if native:
            getInd=lambda x:x.NIndex
        if labelled:
            getInd=lambda x:x.LIndex

        return [anno for anno in self.peakAnnotations if getInd(anno)==index]
    def delAnnotationsForPeakInScan(self, native=False, labelled=False, index=-1):
        if self.hasAnnotationForPeakInScan(native=native, labelled=labelled, index=index):
            todel=[]

            getInd=None

            if native:
                getInd=lambda x:x.NIndex
            if labelled:
                getInd=lambda x:x.LIndex

            for i in range(len(self.peakAnnotations)):
                if getInd(self.peakAnnotations[i])==index:
                    todel.append(i)

            for i in reversed(sorted(todel)):
                self.peakAnnotations.pop(i)
    def hasAnnotationForPeakInScan(self, native=False, labelled=False, index=-1):
        if not(native or labelled) or (native and labelled):
            return

        c=self.getAnnotationsForPeakInScan(native=native, labelled=labelled, index=index)
        return len(self.getAnnotationsForPeakInScan(native=native, labelled=labelled, index=index))>0
    def newAnnotationForPeakInScan(self, nativeIndex=None, labelledIndex=None, anno=None):
        if anno==None:
            anno=Bunch()

        self.peakAnnotations.append(anno)

        return anno

    def getMS2NativeScan(self):
        return (self.nativeMz_list, self.nativeIntensity_list)
    def getMS2LabelledScan(self):
        return (self.labelledMz_list, self.labelledIntensity_list)
    def getAnnotatedMSScan(self, native=False, labelled=False):
        if not(native or labelled) or (native and labelled):
            return

        usePeaks=[]
        for i in range(len(self.nativeMz_list)):
            if self.hasAnnotationForPeakInScan(native=native, labelled=labelled, index=i):
                usePeaks.append(i)
        usePeaks=sorted(list(set(usePeaks)))

        if native:
            return ([self.nativeMz_list[i] for i in usePeaks], [self.nativeIntensity_list[i] for i in usePeaks], [self.getAnnotationsForPeakInScan(native=native, labelled=labelled, index=i) for i in usePeaks])
        if labelled:
            return ([self.labelledMz_list[i] for i in usePeaks], [self.labelledIntensity_list[i] for i in usePeaks], [self.getAnnotationsForPeakInScan(native=native, labelled=labelled, index=i) for i in usePeaks])
    def getMzCorrectedComplementaryMSScanForAnnotatedMSScan(self, native=False, labelled=False, charge=1, isotopeOffset=1.00335):
        if not(native or labelled) or (native and labelled):
            return

        mzs, ints, annos=self.getAnnotatedMSScan(native=not(native), labelled=not(labelled))

        for i, mz, annotations in zip(range(len(mzs)), mzs, annos):
            mzs[i]=mz +(1. if labelled else -1.)*annotations[0].Cn*isotopeOffset/charge

        return (mzs, ints)







class ProcessTarget:
    def __init__(self, targets, chromatogramFile,
                 fullScanEICppm, fullScanThreshold, minMSMSPeakIntensityScaled, scalePrecursorMZ,
                 minXn, useZeroLabelingAtoms, useTracExtractAnnotation,
                 labellingOffset, matchingPPM, maxRelError, annotationElements, annotationPPM, useParentFragmentConsistencyRule,
                 saveAsTSV, saveAsPDF,
                 lock, queue, pID, feVersion):

        self.targets=targets
        self.chromatogramFile=chromatogramFile


        self.labellingOffset=labellingOffset
        self.fullScanEICppm=fullScanEICppm
        self.fullScanThreshold=fullScanThreshold
        self.minMSMSPeakIntensityScaled=minMSMSPeakIntensityScaled
        self.scalePrecursorMZ=scalePrecursorMZ

        self.matchingPPM=matchingPPM
        self.maxRelError=maxRelError

        self.minXn=minXn
        self.useZeroLabelingAtoms=useZeroLabelingAtoms
        self.useTracExtractAnnotation=useTracExtractAnnotation

        self.annotationElements=annotationElements
        self.annotationPPM=annotationPPM
        self.useParentFragmentConsistencyRule=useParentFragmentConsistencyRule

        self.saveAsTSV=saveAsTSV
        self.saveAsPDF=saveAsPDF

        self.lock = lock
        self.queue = queue
        self.pID = pID

        self.feVersion=feVersion

    # Thread safe printing function
    def printMessage(self, message, type="info"):
        if self.lock is not None:
            self.lock.acquire()
            if type.lower()=="info":
                logging.info("   %d: %s" % (self.pID, message))
            elif type.lower()=="error":
                logging.error("   %d: %s" % (self.pID, message))
            elif type.lower()=="warning":
                logging.warning("   %d: %s" % (self.pID, message))
            else:
                logging.info("   %d: %s" % (self.pID, message))
            self.lock.release()

    # helper function used to update the status of the current processing in the Process Dialog
    def postMessageToProgressWrapper(self, mes, val="", forPID=None):
        if forPID is None:
            forPID=self.pID
        if forPID != -1 and self.queue is not None:
            if mes.lower() == "text":
                self.queue.put(Bunch(pid=forPID, mes="text", val="%d: %s" % (forPID, val)))
            elif mes == "value" or mes == "max":
                self.queue.put(Bunch(pid=forPID, mes=mes, val=val))
            elif mes == "start" or mes == "end" or mes == "failed":
                self.queue.put(Bunch(pid=forPID, mes=mes))

    def findClosestScanTo(self, times, startTime):
        bestMatch=-1
        bestMatchRTDiff=100000
        for i in range(len(times)):
            if abs(times[i]-startTime)<bestMatchRTDiff:
                bestMatch=i
                bestMatchRTDiff=abs(times[i]-startTime)

        return bestMatch

    def scaleMSScan(self, msScan, precursorMZ=-1, minVal=0., maxVal=100., minMZ=0, maxMZ=10000000, ppm=5.):
        if len(msScan.intensity_list)==0:
            return msScan
        if len(msScan.intensity_list)==1:
            msScan.intensity_list[0]=100.
            return msScan

        if precursorMZ==-1 or not self.scalePrecursorMZ:
            minInt=minVal
            maxInt=max([msScan.intensity_list[index] for index in range(len(msScan.intensity_list)) if minMZ<=msScan.mz_list[index]<=maxMZ])
        else:
            minInt=minVal
            maxInt=max([msScan.intensity_list[index] for index in range(len(msScan.intensity_list)) if precursorMZ*(1.-ppm/1000000.)<=msScan.mz_list[index]<=precursorMZ*(1.+ppm/1000000.)])

        for i in range(len(msScan.intensity_list)):
            msScan.intensity_list[i]=minVal+(maxVal-minVal)*msScan.intensity_list[i]*1./(maxInt-minInt)

        return msScan

    def removePeaksAboveMZ(self, msScan, maxMZ):
        msScan.intensity_list=[msScan.intensity_list[i] for i in range(len(msScan.mz_list)) if msScan.mz_list[i]<=maxMZ]
        msScan.mz_list=[msScan.mz_list[i] for i in range(len(msScan.mz_list)) if msScan.mz_list[i]<=maxMZ]
        return msScan

    def selectMSMSScansForProcessing(self, mzxml, Cn, minRT, maxRT, nativeMZ, metaboliteCharge, fullScanEvent, nativeMS2ScanEvent, labelledMS2ScanEvent):
        eic, timesFS, scanIDsFS=mzxml.getEIC(nativeMZ, ppm=self.fullScanEICppm, filterLine=fullScanEvent)
        eicL, timesFS, scanIDsFS=mzxml.getEIC(nativeMZ+self.labellingOffset*Cn/metaboliteCharge, ppm=self.fullScanEICppm, filterLine=fullScanEvent)

        TIC, timesMS2Native, scanIdsMS2Native=mzxml.getTIC(filterLine=nativeMS2ScanEvent, useMS2=True)
        TICL, timesMS2Labelled, scanIdsMS2Labelled=mzxml.getTIC(filterLine=labelledMS2ScanEvent, useMS2=True)

        # limit eic to predefined time window
        #print eic
        for i in range(len(eic)):
            if minRT <= timesFS[i]/60. <= maxRT:
                pass
            else:
                eic[i]=0
        maxIntFullScanIndex=max(zip(range(len(eic)), eic), key=lambda x:x[1])[0]

        maxIntMS2NativeIndex=self.findClosestScanTo(timesMS2Native, timesFS[maxIntFullScanIndex])
        maxIntMS2LabelledIndex=self.findClosestScanTo(timesMS2Labelled, timesFS[maxIntFullScanIndex])


        ######################################
        ######################################
        ######################################
        if False:
            plt.plot([t/60. for t in timesFS], eic, color='black')
            plt.plot([t/60. for t in timesFS], [-i for i in eicL], color='black')

            plt.plot([t/60. for t in timesMS2Native], TIC, color='red')
            plt.plot([t/60. for t in timesMS2Labelled], [-i for i in TICL], color='red')

            plt.axvline(x=timesFS[maxIntFullScanIndex]/60., ymin=-10000, ymax = 10000, linewidth=2, color='black')
            try:
                plt.axvline(x=timesMS2Native[maxIntMS2NativeIndex]/60., ymin=-10000, ymax = 10000, linewidth=2, color='yellow')
                plt.axvline(x=timesMS2Labelled[maxIntMS2LabelledIndex]/60., ymin=-10000, ymax = 10000, linewidth=2, color='green')
            except:
                pass

            plt.show()

        return (maxIntFullScanIndex, maxIntMS2NativeIndex, maxIntMS2LabelledIndex,
                timesFS[maxIntFullScanIndex], timesMS2Native[maxIntMS2NativeIndex], timesMS2Labelled[maxIntMS2LabelledIndex],
                scanIDsFS[maxIntFullScanIndex], scanIdsMS2Native[maxIntMS2NativeIndex], scanIdsMS2Labelled[maxIntMS2LabelledIndex])

    def removePeaksBelowThreshold(self, msScan, intThreshold):
        usePeaks=[]
        for i, intensity in zip(range(len(msScan.intensity_list)), msScan.intensity_list):
            if intensity>=intThreshold:
                usePeaks.append(i)

        msScan.mz_list=[msScan.mz_list[i] for i in usePeaks]
        msScan.intensity_list=[msScan.intensity_list[i] for i in usePeaks]

        return msScan

    def cleanScanFromIsotopoes(self, msScan, maxPPM, metaboliteCharge, isotopeOffset=1.00335, use=[+1,-1]):
        if len(msScan.mz_list) in [0,1]:
            return msScan

        toDel=[]
        for i in range(len(msScan.mz_list)-1):
            for j in range(i+1, len(msScan.mz_list)):
                amz=msScan.mz_list[i]
                aint=msScan.intensity_list[i]

                bmz=msScan.mz_list[j]
                bint=msScan.intensity_list[j]


                if (abs(bmz-amz-isotopeOffset/metaboliteCharge)*1000000./amz)<=maxPPM:
                    if +1 in use and bint<aint:
                        toDel.append(j)
                    if -1 in use and aint<bint:
                        toDel.append(i)

        msScan.mz_list=[msScan.mz_list[i] for i in range(len(msScan.mz_list)) if i not in toDel]
        msScan.intensity_list=[msScan.intensity_list[i] for i in range(len(msScan.intensity_list)) if i not in toDel]

        return msScan

    def calculateCn(self, metaboliteCharge, scanMS2Native, scanMS2Labelled, matchingPPM, maxRelError, maxCn, isotopeOffset=1.00335):

        retMS=AnnotatedMSMSSpectra(scanMS2Native=scanMS2Native, scanMS2Labelled=scanMS2Labelled,
                                   charge=metaboliteCharge, ionMode=scanMS2Native.polarity)


        for i in range(len(scanMS2Native.mz_list)):
            bestMatch=Bunch(nInd=-1, lInd=-1, Cn=-1)
            bestRat=10000.
            for j in range(len(scanMS2Labelled.mz_list)):


                nMz=scanMS2Native.mz_list[i]
                nInt=scanMS2Native.intensity_list[i]

                lMz=scanMS2Labelled.mz_list[j]
                lInt=scanMS2Labelled.intensity_list[j]

                atoms=range(self.minXn, maxCn+1)
                if self.useZeroLabelingAtoms:
                    atoms.insert(0, 0)

                _debug=True

                if nMz<(lMz+1):

                    for n in atoms:

                        if n < nMz/12:
                            ppmDiff=abs(lMz-nMz-n*isotopeOffset/metaboliteCharge)*1e6/nMz
                            if ppmDiff < 200 and _debug: print "mz native: {0: >10.3f}, mz labeled: {1: >10.3f}, Cn: {2: >2d}, ppmdiff: {3: >10.3f}, nInt: {4: >10.3f}, lInt: {5: >10.3f}".format(nMz, lMz, n, ppmDiff, nInt, lInt),
                            # check if delta mz is explained by n carbon atoms
                            if ppmDiff<=matchingPPM:

                                # check if intensity ratios are within expected error window
                                if abs(lInt-nInt)<=maxRelError:

                                    if abs(lInt-nInt)<bestRat:
                                        bestMatch=Bunch(nInd=i, lInd=j, Cn=n)
                                        bestRat=abs(lInt-nInt)
                                        if ppmDiff < 200 and _debug: print " ***",
                            if ppmDiff < 200 and _debug: print ""

            if bestMatch.nInd!=-1:
                retMS.newAnnotationForPeakInScan(nativeIndex=bestMatch.nInd,labelledIndex=bestMatch.lInd,
                                                 anno=Bunch(Cn=bestMatch.Cn, NIndex=bestMatch.nInd, LIndex=bestMatch.lInd, ratioError=bestRat))

        return retMS
    def calculateOptimalMatch(self, scanMS2Native, scanMS2Labelled, matchingPPM, maxRelError, maxCn, charge, isotopeOffset=1.00335):
        scanMS2Annotated=self.calculateCn(charge, scanMS2Native, scanMS2Labelled, matchingPPM=matchingPPM, maxRelError=maxRelError, maxCn=maxCn, isotopeOffset=isotopeOffset)

        for lUsedPeakIndex in [anno.LIndex for anno in scanMS2Annotated.peakAnnotations]:
            annosFor=[anno for anno in scanMS2Annotated.peakAnnotations if anno.LIndex==lUsedPeakIndex]
            bAnno=min(annosFor, key=lambda x:x.ratioError)
            t=list(set(scanMS2Annotated.peakAnnotations).difference(set(annosFor)))
            t.append(bAnno)
            scanMS2Annotated.peakAnnotations=t

        return scanMS2Annotated

    def removeNonAnnotatedPeaks(self, scanAnnotated):
        usedNativeIndices=sorted(list(set([anno.NIndex for anno in scanAnnotated.peakAnnotations])))
        usedLabelledIndices=sorted(list(set([anno.LIndex for anno in scanAnnotated.peakAnnotations])))

        scanAnnotated.nativeMz_list=[scanAnnotated.nativeMz_list[i] for i in usedNativeIndices]
        scanAnnotated.nativeIntensity_list=[scanAnnotated.nativeIntensity_list[i] for i in usedNativeIndices]
        newMappingsNative={}
        for i in range(len(scanAnnotated.nativeMz_list)):
            newMappingsNative[usedNativeIndices[i]]=i

        scanAnnotated.labelledMz_list=[scanAnnotated.labelledMz_list[i] for i in usedLabelledIndices]
        scanAnnotated.labelledIntensity_list=[scanAnnotated.labelledIntensity_list[i] for i in usedLabelledIndices]
        newMappingsLabelled={}
        for i in range(len(scanAnnotated.labelledMz_list)):
            newMappingsLabelled[usedLabelledIndices[i]]=i

        for anno in scanAnnotated.peakAnnotations:
            anno.NIndex=newMappingsNative[anno.NIndex]
            anno.LIndex=newMappingsLabelled[anno.LIndex]

        return scanAnnotated

    def annotatePeaksWithSumFormulas(self, adductObj, scanAnnotated, useAtoms, atomsRange, fixed, ppm, useSevenGoldenRules=True):
        sfg=sumFormulaGenerator()
        ft=formulaTools()

        for annotation in scanAnnotated.peakAnnotations:
            annotation.generatedSumFormulas={}

            uA=deepcopy(useAtoms)
            ar=deepcopy(atomsRange)
            f=deepcopy(fixed)

            if uA[0]=="C":
                if self.useTracExtractAnnotation:
                    ar[0]=(annotation.Cn, int(floor((scanAnnotated.nativeMz_list[annotation.NIndex]-adductObj.mzoffset/adductObj.charge)/12.)))
                    del f[f.index('C')]
                else:
                    ar[0]=annotation.Cn

            sfs=sfg.findFormulas(scanAnnotated.nativeMz_list[annotation.NIndex]-adductObj.mzoffset/adductObj.charge, ppm=ppm,
                                 useAtoms=uA, atomsRange=ar, fixed=f, useSevenGoldenRules=useSevenGoldenRules, useSecondRule=False)


            if len(sfs)>0:
                annotation.generatedSumFormulas[adductObj.name]=[]
                for sf in sfs:
                    elems=ft.parseFormula(sf)
                    massSF=ft.calcMolWeight(elems)-adductObj.mzoffset/adductObj.charge
                    dppm=(scanAnnotated.nativeMz_list[annotation.NIndex]-massSF)*1000000./massSF

                    annotation.generatedSumFormulas[adductObj.name].append(Bunch(sumFormula=sf, deltaPPM=dppm, neutralLossToParent=""))


    def calcNeutralLosses(self, scanAnnotated, parentSumFormula, ppm):
        assert isinstance(parentSumFormula, str) and parentSumFormula!=""

        fT=formulaTools()

        elemsParent=fT.parseFormula(parentSumFormula)

        for annotation in scanAnnotated.peakAnnotations:
            for add, sfs in annotation.generatedSumFormulas.items():
                for sf in sfs:
                    elemsFragment=fT.parseFormula(sf.sumFormula)
                    elemsNeutralLoss=fT.calcDifferenceBetweenElemDicts(elemsFragment, elemsParent)
                    if len(elemsNeutralLoss)==0 and fT.flatToString(elemsParent)==fT.flatToString(elemsFragment):
                        sf.neutralLossToParent="Parent"
                    else:
                        mParent=fT.calcMolWeight(elemsParent)
                        mFragment=fT.calcMolWeight(elemsFragment)
                        mNeutralLoss=fT.calcMolWeight(elemsNeutralLoss)

                        diffUnexplained=mParent-(mFragment+mNeutralLoss)

                        if abs(diffUnexplained*1000000./mParent) <= 2*ppm:
                            sf.neutralLossToParent=fT.flatToString(elemsNeutralLoss)
                        else:
                            sf.neutralLossToParent="Unknown"


    def calcParentSumFormulas(self, nativeMZ, Cn, adductObj, useAtoms, atomsRange, fixed, ppm, useSevenGoldenRules=True):
        sfg=sumFormulaGenerator()

        mParent=nativeMZ-adductObj.mzoffset/adductObj.charge

        uA=deepcopy(useAtoms)
        ar=deepcopy(atomsRange)
        f=deepcopy(fixed)
        if uA[0]=="C":
            ar[0]=Cn

        sfs=sfg.findFormulas(mParent, ppm=ppm, useAtoms=uA, atomsRange=ar, fixed=f, useSevenGoldenRules=useSevenGoldenRules, useSecondRule=False)

        return sfs

    def checkIfSubset(self, elemsFragment, elemsParent):
        isPossible=True
        for element, count in elemsFragment.items():
            if element in elemsParent.keys() and count <=elemsParent[element]:
                pass
            else:
                isPossible=False
        return isPossible

    def checkForFragmentParentSumFormulaConsistency(self, scanMS2Annotated, parentSumFormula, removeImpossibleAnnotations=True):
        fT=formulaTools()
        ret=[]

        annos=[scanMS2Annotated.getAnnotationsForPeakInScan(native=True, index=i)[0] for i in range(len(scanMS2Annotated.nativeMz_list))]

        if not isinstance(parentSumFormula, list):
            parentSumFormula=[parentSumFormula]

        if removeImpossibleAnnotations:
            for anno in annos:
                for adduct in anno.generatedSumFormulas.keys():

                    todel=[]
                    for index, annoSF in enumerate(anno.generatedSumFormulas[adduct]):
                        annoSFElems=fT.parseFormula(annoSF.sumFormula)

                        useFragmentAnnotation=False
                        for ps in parentSumFormula:
                            psElems=fT.parseFormula(ps)
                            isPossible=self.checkIfSubset(annoSFElems, psElems)
                            useFragmentAnnotation=useFragmentAnnotation or isPossible

                        if not useFragmentAnnotation:
                            todel.append(index)

                    if len(todel)>0:
                        todel=list(set(todel))
                        todel=sorted(todel)
                        todel.reverse()

                        for delIndex in todel:
                            anno.generatedSumFormulas[adduct].pop(delIndex)
                        if len(anno.generatedSumFormulas[adduct])==0:
                            del anno.generatedSumFormulas[adduct]



        return ret


    def saveResultsToTSV(self, target, curs):
        annotatedSpectrumNative=[p for p in SQLSelectAsObject(curs,
                            selectStatement="SELECt mzs, ints, annos AS annos FROM MSSpectra WHERE forTarget=%d AND type='native_cleaned'"%target.id)][0]
        annotatedSpectrumLabeled=[p for p in SQLSelectAsObject(curs,
                            selectStatement="SELECt mzs, ints, annos AS annos FROM MSSpectra WHERE forTarget=%d AND type='labelled_cleaned'"%target.id)][0]

        if len(annotatedSpectrumNative.mzs)>0:
            annotatedSpectrumNative.mzs=[float(mz) for mz in annotatedSpectrumNative.mzs.split(",")]
        else:
            annotatedSpectrumNative.mzs=[]

        if len(annotatedSpectrumNative.ints)>0:
            annotatedSpectrumNative.ints=[float(i) for i in annotatedSpectrumNative.ints.split(",")]
        else:
            annotatedSpectrumNative.ints=[]

        if len(annotatedSpectrumNative.annos)>0:
            annotatedSpectrumNative.annos=pickle.loads(base64.b64decode(annotatedSpectrumNative.annos))
        else:
            annotatedSpectrumNative.annos=[]

        if len(annotatedSpectrumLabeled.mzs)>0:
            annotatedSpectrumLabeled.mzs=[float(mz) for mz in annotatedSpectrumLabeled.mzs.split(",")]
        else:
            annotatedSpectrumLabeled.mzs=[]

        if len(annotatedSpectrumLabeled.ints)>0:
            annotatedSpectrumLabeled.ints=[float(i) for i in annotatedSpectrumLabeled.ints.split(",")]
        else:
            annotatedSpectrumLabeled.ints=[]

        if len(annotatedSpectrumLabeled.annos)>0:
            annotatedSpectrumLabeled.annos=pickle.loads(base64.b64decode(annotatedSpectrumLabeled.annos))
        else:
            annotatedSpectrumLabeled.annos=[]

        with open(target.lcmsmsFile + "." + target.targetName + ".tsv", "wb") as fOut:

            fOut.write("\t".join(["Num", "MZ", "L_MZ", "D_MZ_ppm", "relInt", "L_relInt", "Cn", "sumFormula", "Adduct", "NeutralLoss"]))
            fOut.write("\n")

            lines=[]
            mzs=[]
            for index in range(len(annotatedSpectrumNative.mzs)):
                line=[]
                indexLab=annotatedSpectrumNative.annos[index].LIndex
                b = Bunch(index=index, mz=annotatedSpectrumNative.mzs[index], relInt=annotatedSpectrumNative.ints[index], l_mz=annotatedSpectrumLabeled.mzs[indexLab], l_relInt=annotatedSpectrumLabeled.ints[indexLab],
                          Cn=annotatedSpectrumNative.annos[index].Cn, sumFormulas=annotatedSpectrumNative.annos[index].generatedSumFormulas)
                b.dmzppm=(b.l_mz-b.mz-b.Cn*1.00335)*1000000./b.mz
                ##(generatedSumFormulas:{'[M]+': [(neutralLossToParent:,deltaPPM:-12.0126539494,sumFormula:C4H5O2N3)]},ratioError:6.78271579018,NIndex:0,Cn:0,LIndex:0)

                totSumForms=0
                for key, val in b.sumFormulas.items():
                    totSumForms+=len(val)
                    if len(val)==0:
                        del b.sumFormulas[key]

                line.append("\t".join([str(b.index), "%.4f"%b.mz, "%.4f"%b.l_mz, "%.4f"%b.dmzppm,  "%.1f"%b.relInt, "%.1f"%b.l_relInt, "%d"%b.Cn]))
                line.append("\t")
                mzs.append(b.mz)

                if totSumForms==0:
                    pass
                elif totSumForms==1:
                    line.append("\t".join([str(list(b.sumFormulas.values())[0][0].sumFormula),
                                          # "%.2f"%list(b.sumFormulas.values())[0][0].deltaPPM,
                                          str(list(b.sumFormulas.keys())[0]),
                                          str(list(b.sumFormulas.values())[0][0].neutralLossToParent)]))
                else:
                    line.append("*")
                    for adduct in b.sumFormulas.keys():
                        if len(b.sumFormulas[adduct])>0:
                            for sumForm in b.sumFormulas[adduct]:
                                line.append("\n")
                                line.append("\t".join(["", "", "", "", "", "", "", sumForm.sumFormula, adduct, sumForm.neutralLossToParent]))
                lines.append(line)

            for i, v in sorted(enumerate(mzs), key=lambda x:x[1]):
                fOut.write("".join(lines[i]))
                fOut.write("\n")



    def saveResultsToPDF(self, target, curs):

        pdf = canvas.Canvas(target.lcmsmsFile + "." + target.targetName + ".pdf", pagesize=pagesizes.A4)

        fT=formulaTools()

        ### write processing settings
        def writeKeyValuePair(key, value, pdf, currentHeight):
            pdf.drawString(50, currentHeight, key)
            pdf.drawString(240, currentHeight, value)
            return currentHeight - 18

        currentHeight = 800
        currentHeight=writeKeyValuePair("ID", str(target.id), pdf, currentHeight)
        currentHeight=writeKeyValuePair("Target name", target.targetName, pdf, currentHeight)
        currentHeight=writeKeyValuePair("LC-HRMS/MS file", target.lcmsmsFile[(target.lcmsmsFile.rfind("/")+1):], pdf, currentHeight)
        currentHeight=writeKeyValuePair("Sum formula", target.parentSumFormula, pdf, currentHeight)
        currentHeight=writeKeyValuePair("Adduct", "%s (%.4f)"%(target.adduct, target.adductMZOffset), pdf, currentHeight)

        pdf.line(50, currentHeight + 20 - 4, 540, currentHeight + 20 - 4); currentHeight-=10

        currentHeight=writeKeyValuePair("Native precursor m/z", "%.4f"%target.precursorMZ, pdf, currentHeight)
        currentHeight=writeKeyValuePair("Cn", str(target.Cn), pdf, currentHeight)
        currentHeight=writeKeyValuePair("Charge count", str(target.chargeCount), pdf, currentHeight)
        currentHeight=writeKeyValuePair("Retention time [min]", "%.2f - %.2f"%(target.startRT, target.stopRT), pdf, currentHeight)
        currentHeight=writeKeyValuePair("MS scan event", target.scanEventMS1, pdf, currentHeight)
        currentHeight=writeKeyValuePair("Native MS/MS scan event", target.scanEventMS2Native, pdf, currentHeight)
        currentHeight=writeKeyValuePair("Labelled MS/MS scan event", target.scanEventMS2Labelled, pdf, currentHeight)
        currentHeight=writeKeyValuePair("Selected native MS/MS scan id", str(target.scanIDNativeRaw), pdf, currentHeight)
        currentHeight=writeKeyValuePair("Selected labelled MS/MS scan id", str(target.scanIDLabelledRaw), pdf, currentHeight)

        pdf.line(50, currentHeight + 20 - 4, 540, currentHeight + 20 - 4); currentHeight-=10

        currentHeight=writeKeyValuePair("FragExtract version",  str(self.feVersion), pdf, currentHeight)




        nativeRawSpectrum=[p for p in SQLSelectAsObject(curs, selectStatement="SELECT mzs, ints FROM MSSpectra WHERE forTarget=%d AND type='native_raw'"%target.id)][0]
        if len(nativeRawSpectrum.mzs)>0:
            nativeRawSpectrum.mzs=[float(mz) for mz in nativeRawSpectrum.mzs.split(",")]
        else:
            nativeRawSpectrum.mzs=[]
        if len(nativeRawSpectrum.ints)>0:
            nativeRawSpectrum.ints=[float(i) for i in nativeRawSpectrum.ints.split(",")]
        else:
            nativeRawSpectrum.ints=[]

        labelledRawSpectrum=[p for p in SQLSelectAsObject(curs, selectStatement="SELECT mzs, ints FROM MSSpectra WHERE forTarget=%d AND type='labelled_raw'"%target.id)][0]
        if len(labelledRawSpectrum.mzs)>0:
            labelledRawSpectrum.mzs=[float(mz) for mz in labelledRawSpectrum.mzs.split(",")]
        else:
            labelledRawSpectrum.mzs=[]
        if len(labelledRawSpectrum.ints)>0:
            labelledRawSpectrum.ints=[float(i) for i in labelledRawSpectrum.ints.split(",")]
        else:
            labelledRawSpectrum.ints=[]

        annotatedSpectrum=[p for p in SQLSelectAsObject(curs,
                        selectStatement="SELECT mzs, ints, annos AS annos FROM MSSpectra WHERE forTarget=%d AND type='native_cleaned'"%target.id)][0]
        if len(annotatedSpectrum.mzs)>0:
            annotatedSpectrum.mzs=[float(mz) for mz in annotatedSpectrum.mzs.split(",")]
        else:
            annotatedSpectrum.mzs=[]
        if len(annotatedSpectrum.ints)>0:
            annotatedSpectrum.ints=[float(i) for i in annotatedSpectrum.ints.split(",")]
        else:
            annotatedSpectrum.ints=[]
        if len(annotatedSpectrum.annos)>0:
            annotatedSpectrum.annos=pickle.loads(base64.b64decode(annotatedSpectrum.annos))
        else:
            annotatedSpectrum.annos=[]

        eics=[eic for eic in
              SQLSelectAsObject(curs,
                                selectStatement="SELECT intensityList, timesList, forMZ, type FROM EICs WHERE forTarget=%d"%target.id)]

        data=[]

        data.append(["Num", "mz", "relInt", "Cn", "sumFormula", "deviation [ppm]", "Adduct", "NeutralLoss"])
        widths=[40,60,30,30,80,50,80]
        for index in range(len(annotatedSpectrum.mzs)):
            b = Bunch(index=index, mz=annotatedSpectrum.mzs[index], relInt=annotatedSpectrum.ints[index],
                      Cn=annotatedSpectrum.annos[index].Cn, sumFormulas=annotatedSpectrum.annos[index].generatedSumFormulas)

            totSumForms=0
            for key, val in b.sumFormulas.items():
                totSumForms+=len(val)
                if len(val)==0:
                    del b.sumFormulas[key]

            dl=[str(b.index), "%.4f"%b.mz, "%.1f"%b.relInt, "%d"%b.Cn]

            if totSumForms==0:
                dl.extend(["","","",""])
                data.append(dl)
            elif totSumForms==1:
                m=fT.calcMolWeight(fT.parseFormula(list(b.sumFormulas.values())[0][0].sumFormula))+target.adductObj.mzoffset/target.adductObj.charge
                dl.extend([str(list(b.sumFormulas.values())[0][0].sumFormula),
                           "%.2f"%((m-b.mz)*1000000./b.mz),
                           str(list(b.sumFormulas.keys())[0]),
                           str(list(b.sumFormulas.values())[0][0].neutralLossToParent)])
                data.append(dl)
            else:
                dl.extend(["*","","",""])
                data.append(dl)
                for adduct in b.sumFormulas.keys():
                    if len(b.sumFormulas[adduct])>0:
                        for sumForm in b.sumFormulas[adduct]:
                            m=fT.calcMolWeight(fT.parseFormula(sumForm.sumFormula))+target.adductObj.mzoffset/target.adductObj.charge
                            dl=["", "", "", "", sumForm.sumFormula, "%.2f"%((m-b.mz)*1000000./b.mz), adduct, sumForm.neutralLossToParent]
                            data.append(dl)

        style = [('FONTSIZE', (0, 0), (-1, -1), 8),
                 ('LEFTPADDING', (0, 0), (-1, -1), 0),
                 ('RIGHTPADDING', (0, 0), (-1, -1), 0),
                 ('TOPPADDING', (0, 0), (-1, -1), 0),
                 ('BOTTOMPADDING', (0, 0), (-1, -1), 0),
                 ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                 ('ALIGN', (0, 0), (-1, -1), 'LEFT')]
        table = Table(data, colWidths=widths, style=style)

        table.wrapOn(pdf, 2 * len(data), 2 * len(data))
        table.drawOn(pdf, 50, currentHeight-30-len(data)*9)

        ## Plot MSMS spectra - native MSMS
        drawing = Drawing(500, 260)
        lp = LinePlot()
        lp.x = 20
        lp.y = 0
        lp.height = 260
        lp.width = 400

        dd=[]
        colors=[]
        strokeWidth=[]

        ## intensity threshold
        dd.append([(min(min(nativeRawSpectrum.mzs), min(labelledRawSpectrum.mzs)), self.minMSMSPeakIntensityScaled), (max(max(nativeRawSpectrum.mzs), max(labelledRawSpectrum.mzs)), self.minMSMSPeakIntensityScaled)])
        colors.append(Color(100 / 255., 100 / 255., 100 / 255.))
        strokeWidth.append(0.3)
        dd.append([(min(min(nativeRawSpectrum.mzs), min(labelledRawSpectrum.mzs)), -self.minMSMSPeakIntensityScaled), (max(max(nativeRawSpectrum.mzs), max(labelledRawSpectrum.mzs)), -self.minMSMSPeakIntensityScaled)])
        colors.append(Color(50 / 255., 50 / 255., 50 / 255.))
        strokeWidth.append(0.3)

        ## parent ions
        dd.append([(target.precursorMZ, 110), (target.precursorMZ, 120)])
        colors.append(Color(107 / 255., 142 / 255., 35 / 255.))
        strokeWidth.append(2.1)
        dd.append([(target.precursorMZ+1.00335*target.Cn/target.chargeCount, -110), (target.precursorMZ+1.00335*target.Cn/target.chargeCount, -120)])
        colors.append(Color(107 / 255., 142 / 255., 35 / 255.))
        strokeWidth.append(2.1)
        dd.append([(target.precursorMZ, 0), (target.precursorMZ, 120)])
        colors.append(Color(107 / 255., 142 / 255., 35 / 255.))
        strokeWidth.append(0.1)
        dd.append([(target.precursorMZ+1.00335*target.Cn/target.chargeCount, 0), (target.precursorMZ+1.00335*target.Cn/target.chargeCount, -120)])
        colors.append(Color(107 / 255., 142 / 255., 35 / 255.))
        strokeWidth.append(0.1)

        ## matched fragment peaks
        for mz, it, anno in zip(annotatedSpectrum.mzs, annotatedSpectrum.ints, annotatedSpectrum.annos):
            dd.append([(mz, 0), (mz, it)])
            colors.append(Color(178 / 255., 34 / 255., 34 / 255.))
            strokeWidth.append(1.2)
            dd.append([(mz+1.00335*anno.Cn/target.chargeCount, 0), (mz+1.00335*anno.Cn/target.chargeCount, -it)])
            colors.append(Color(178 / 255., 34 / 255., 34 / 255.))
            strokeWidth.append(1.2)

        ## recorded ions
        for mz, it in zip(nativeRawSpectrum.mzs, nativeRawSpectrum.ints):
            dd.append([(mz, 0), (mz, it)])
            colors.append(Color(47 / 255., 79 / 255., 79 / 255.))
            strokeWidth.append(.35)
        for mz, it in zip(labelledRawSpectrum.mzs, labelledRawSpectrum.ints):
            dd.append([(mz, 0), (mz, -it)])
            colors.append(Color(47 / 255., 79 / 255., 79 / 255.))
            strokeWidth.append(.35)

        lp.data = dd
        lp.joinedLines = 1
        for i in range(len(dd)):
            lp.lines[i].strokeColor = colors[i]
            lp.lines[i].strokeWidth = strokeWidth[i]

        drawing.add(lp)
        renderPDF.draw(drawing, pdf, 15 , 25)
        pdf.drawString(30, 295, "Mass spectra")




        ## Plot EICs
        drawing = Drawing(500, 260)
        lp = LinePlot()
        lp.x = 460
        lp.y = 0
        lp.height = 260
        lp.width = 100

        dd=[]
        colors=[]
        strokeWidth=[]

        for eic in eics:
            times=[float(f) for f in eic.timesList.split(",")]
            intensities=[float(f) for f in eic.intensityList.split(",")]

            times=[t/60. for t in times]

            m=max(intensities)
            if m==0:
                m=1
            intensities=[i/m for i in intensities]

            if eic.type=="Labeled":
                intensities=[-i for i in intensities]

            dd.append(zip(times, intensities))
            colors.append(Color(47 / 255., 79 / 255., 79 / 255.))
            strokeWidth.append(.35)

        lp.data = dd
        lp.joinedLines = 1
        for i in range(len(dd)):
            lp.lines[i].strokeColor = colors[i]
            lp.lines[i].strokeWidth = strokeWidth[i]

        drawing.add(lp)
        renderPDF.draw(drawing, pdf, 15 , 25)
        pdf.drawString(470, 295, "EICs")

        pdf.save()





    def process(self):
        try:

            for index, target in enumerate(self.targets):
                self.pID=target.pID

                if index==0:
                    self.postMessageToProgressWrapper("start", forPID=target.pID)
                    self.postMessageToProgressWrapper("max", 100., forPID=target.pID)
                    self.postMessageToProgressWrapper("text", "Parsing file %s"%(self.chromatogramFile[(self.chromatogramFile.rfind("/")+1):]), forPID=target.pID)


            # region parse mzxml file
            chromatogram=Chromatogram()
            chromatogram.parse_file(self.chromatogramFile)
            # endregion
            # region create results db
            conn=sqlite3.connect(self.chromatogramFile+".identified.sqlite")
            curs=conn.cursor()
            #endregion

            target=self.targets[0]
            createTableFromBunch("Targets", target, cursor=curs, primaryKeys=["id"], autoIncrements=["id"], ifNotExists=True)
            createTableFromBunch("MSSpectra", Bunch(id=int, forTarget=int, mzs=str, ints=str, annos=str, scanTime=float, scanID=str, type=str), cursor=curs,
                                 primaryKeys=["id"], autoIncrements=["id"], ifNotExists=True)
            createTableFromBunch("EICs", Bunch(id=int, forTarget=int, forMZ=float, type=str, timesList=str, intensityList=str), cursor=curs, primaryKeys=["id"], autoIncrements=["id"], ifNotExists=True)
            conn.commit()

            for index, target in enumerate(self.targets):
                self.pID=target.pID

                if index!=0:
                    self.postMessageToProgressWrapper("start", forPID=target.pID)
                    self.postMessageToProgressWrapper("max", 100., forPID=target.pID)
                self.postMessageToProgressWrapper("value", 40.+60.*index/len(self.targets), forPID=target.pID)
                self.postMessageToProgressWrapper("text", "Processing file %s, target %.4f"%(self.chromatogramFile[(self.chromatogramFile.rfind("/")+1):], target.precursorMZ), forPID=target.pID)

                self._process(chromatogram, target, conn, curs)

                self.postMessageToProgressWrapper("end", forPID=target.pID)

        except Exception as ex:
            import traceback
            traceback.print_exc()
            logging.error(str(traceback))

            self.printMessage("Error in %s: %s" % (self.chromatogramFile, str(ex)))
            self.postMessageToProgressWrapper("failed")

    def _process(self, chromatogram, target, conn, curs):
        self.printMessage("processing file: %s; target: %.4f"%(self.chromatogramFile, target.precursorMZ))

        useAtoms=[]
        atomsRange=[]
        fixed=[]
        if "C" in [ce.name for ce in self.annotationElements]:
            useAtoms.append("C")
            atomsRange.append(-1)
            fixed=["C"]
        for ce in self.annotationElements:
            if ce.name != "C":
                useAtoms.append(ce.name)
                atomsRange.append((ce.minCount, ce.maxCount))


        # region 1. Processing step: get scan events and select the user specified ones
        scanEventsMS2=sorted(list(chromatogram.getFilterLines(includeMS1=False, includeMS2=True, includePosPolarity=True, includeNegPolarity=True)))
        scanEventsMS1=sorted(list(chromatogram.getFilterLines(includeMS1=True, includeMS2=False, includePosPolarity=True, includeNegPolarity=True)))
        nativeScanEventString=scanEventsMS2[target.nativeScanEventNum]
        labelledScanEventString=scanEventsMS2[target.labelledScanEventNum]
        fullScanEventString=scanEventsMS1[target.fullScanEventNum]
        target.scanEventMS1=fullScanEventString
        target.scanEventMS2Native=nativeScanEventString
        target.scanEventMS2Labelled=labelledScanEventString
        # endregion

        # region 2. create results tables
        # write target to results db
        id=writeObjectAsSQLInsert(curs, "Targets", target, writeFields=["targetName", "lcmsmsFile", "precursorMZ", "Cn", "chargeCount", "startRT", "parentSumFormula",
                                                                       "stopRT", "scanEventMS1", "scanEventMS2Native", "scanEventMS2Labelled",
                                                                       "adduct", "adductMZOffset"])
        target.id=id
        # endregion

        # region 3. Processing step: select the most abundant FullScan and select the two following MS2 scans
        t=self.selectMSMSScansForProcessing(chromatogram, target.Cn, target.startRT, target.stopRT, target.precursorMZ, target.chargeCount, fullScanEventString, nativeScanEventString, labelledScanEventString)
        maxIntFullScanIndex, maxIntMS2NativeIndex, maxIntMS2LabelledIndex, \
            maxIntFullScanTime, maxIntMS2NativeTime, maxIntMS2LabelledTime,\
            maxIntFullScanID, maxIntMS2NativeID, maxIntMS2LabelledID=t

        scanFS=deepcopy(chromatogram.getScanByID(maxIntFullScanID))
        scanMS2Native=deepcopy(chromatogram.getScanByID(maxIntMS2NativeID))
        scanMS2Labelled=deepcopy(chromatogram.getScanByID(maxIntMS2LabelledID))
        # endregion

        if scanMS2Native.polarity != scanMS2Labelled.polarity:
            logging.error(" -> Error: Polarity of native and labelled scan are not the same")
            return
        if scanMS2Native.precursorCharge != scanMS2Labelled.precursorCharge:
            logging.error(" -> Error: Charge of native and labelled scan are not the same")
            return

        # region 4. Processing step: remove everything with a higher mass than the precursor
        #scanFS=self.removePeaksAboveMZ(scanFS, target.precursorMZ+self.labellingOffset*target.Cn/target.chargeCount)
        #scanMS2Native=self.removePeaksAboveMZ(scanMS2Native, target.precursorMZ)
        #scanMS2Labelled=self.removePeaksAboveMZ(scanMS2Labelled, target.precursorMZ+self.labellingOffset*target.Cn/target.chargeCount)
        # endregion

        # region 5. Processing step: scale MS scans
        scanMS2Native=self.scaleMSScan(scanMS2Native, maxMZ=target.precursorMZ+0.5, precursorMZ=target.precursorMZ, ppm=self.matchingPPM)
        scanMS2Labelled=self.scaleMSScan(scanMS2Labelled, maxMZ=target.precursorMZ+self.labellingOffset*target.Cn/target.chargeCount+.5, precursorMZ=target.precursorMZ+self.labellingOffset*target.Cn/target.chargeCount, ppm=self.matchingPPM)

        assert maxIntFullScanTime==scanFS.retention_time and maxIntMS2NativeTime==scanMS2Native.retention_time and maxIntMS2LabelledTime==scanMS2Labelled.retention_time
        assert isinstance(scanMS2Native, MS2Scan) and isinstance(scanMS2Labelled, MS2Scan)

        # endregion
        # region write native spectra to results db
        writeObjectAsSQLInsert(curs, "MSSpectra", Bunch(forTarget=target.id, mzs=",".join([str(mz) for mz in scanMS2Native.mz_list]),
                                                        ints=",".join([str(i) for i in scanMS2Native.intensity_list]),
                                                        scanTime=scanMS2Native.retention_time, scanID=str(maxIntMS2NativeID), type="native_raw"),
                               writeFields=["forTarget", "mzs", "ints", "scanTime", "type"])
        # write labelled spectra to results db
        writeObjectAsSQLInsert(curs, "MSSpectra", Bunch(forTarget=target.id, mzs=",".join([str(mz) for mz in scanMS2Labelled.mz_list]),
                                                        ints=",".join([str(i) for i in scanMS2Labelled.intensity_list]),
                                                        scanTime=scanMS2Labelled.retention_time, scanID=str(maxIntMS2LabelledID), type="labelled_raw"),
                               writeFields=["forTarget", "mzs", "ints", "scanTime", "type"])

        curs.execute("UPDATE Targets SET scanIDNativeRaw='%s', scanIDLabelledRaw='%s' WHERE id=%d"%(str(maxIntMS2NativeID), str(maxIntMS2LabelledID), target.id))
        target.scanIDNativeRaw=str(maxIntMS2NativeID)
        target.scanIDLabelledRaw=str(maxIntMS2LabelledID)

        conn.commit() 
        # endregion

        # region 6. Processing step: apply threshold and remove isotopolog peaks
        scanMS2Native=self.removePeaksBelowThreshold(scanMS2Native, intThreshold=self.minMSMSPeakIntensityScaled)
        scanMS2Native=self.cleanScanFromIsotopoes(scanMS2Native, maxPPM=self.matchingPPM, metaboliteCharge=target.chargeCount, isotopeOffset=1.00335, use=[+1])

        scanMS2Labelled=self.removePeaksBelowThreshold(scanMS2Labelled, intThreshold=self.minMSMSPeakIntensityScaled)
        scanMS2Labelled=self.cleanScanFromIsotopoes(scanMS2Labelled, maxPPM=self.matchingPPM, metaboliteCharge=target.chargeCount, isotopeOffset=1.00335, use=[-1])
        # endregion

        # region 7. Processing step: calculate Cn for each fragment peak
        scanMS2Annotated=self.calculateOptimalMatch(scanMS2Native, scanMS2Labelled, self.matchingPPM, self.maxRelError, target.Cn, target.chargeCount, isotopeOffset=1.00335)
        # endregion

        # region 8. Processing step: select only matched peaks
        self.removeNonAnnotatedPeaks(scanMS2Annotated)
        # endregion

        # region 9. Processing step: annotate peaks with sum formulas
        self.annotatePeaksWithSumFormulas(target.adductObj, scanMS2Annotated, useAtoms, atomsRange, fixed, ppm=self.annotationPPM)
        # endregion

        # region 10. Processing step: calculate putative sum formulas for parent
        genSumForms=self.calcParentSumFormulas(target.precursorMZ, target.Cn, target.adductObj, useAtoms, atomsRange, fixed, ppm=self.annotationPPM)

        # endregion

        # region
        # use either generated sum formulas or check, if the provided is among the generated ones
        if target.parentSumFormula=="" or target.parentSumFormula==None:
            if len(genSumForms)==1:
                target.parentSumFormula=genSumForms[0]
                hasUniqueParentSumFormula=True
            else:
                target.parentSumFormula=genSumForms
                hasUniqueParentSumFormula=False

            if self.useParentFragmentConsistencyRule: # check for fragment-parent-sumFormula consistency
                self.checkForFragmentParentSumFormulaConsistency(scanMS2Annotated, target.parentSumFormula, removeImpossibleAnnotations=hasUniqueParentSumFormula)
        else:
            fT=formulaTools()
            elemsProvided=fT.parseFormula(target.parentSumFormula)
            target.parentSumFormula=fT.flatToString(elemsProvided)
            alsoGenerated=False

            for sf in genSumForms:
                elemsSF=fT.parseFormula(sf)

                if fT.flatToString(elemsSF)==fT.flatToString(elemsProvided):
                    alsoGenerated=True

            if not alsoGenerated:
                logging.warning("WARNING: The provided sum formula (%s) for this target was not generated by the seven golden rules"%(fT.flatToString(elemsProvided)), ",".join(genSumForms))

            hasUniqueParentSumFormula=True

            if self.useParentFragmentConsistencyRule: # check for fragment-parent-sumFormula consistency
                self.checkForFragmentParentSumFormulaConsistency(scanMS2Annotated, target.parentSumFormula, removeImpossibleAnnotations=True)

        target.parentSumFormula=str(target.parentSumFormula)


        sumFormsDisplay=target.parentSumFormula.replace("'", "").replace("\"", "").replace("[", "").replace("]", "")
        curs.execute("UPDATE Targets SET parentSumFormula='%s' WHERE id=%d"%(sumFormsDisplay, target.id))
        # endregion

        # region 11. Processing step: Calculate neutral losses for MSMS peaks
        if hasUniqueParentSumFormula and isinstance(target.parentSumFormula, str) and target.parentSumFormula!="":  ## neutral losses will only be calculated for provided or unique sum formulas
            self.calcNeutralLosses(scanMS2Annotated, target.parentSumFormula, ppm=self.annotationPPM)
        # endregion

        # region write cleaned spectra to results db
        writeObjectAsSQLInsert(curs, "MSSpectra", Bunch(forTarget=target.id, mzs=",".join([str(mz) for mz in scanMS2Annotated.nativeMz_list]),
                                                        ints=",".join([str(i) for i in scanMS2Annotated.nativeIntensity_list]),
                                                        annos=base64.b64encode(pickle.dumps([scanMS2Annotated.getAnnotationsForPeakInScan(native=True, index=i)[0] for i in range(len(scanMS2Annotated.nativeMz_list))])),
                                                        scanTime=scanMS2Native.retention_time, type="native_cleaned"),
                               writeFields=["forTarget", "mzs", "ints", "annos", "scanTime", "type"])
        # write labelled spectra to results db
        writeObjectAsSQLInsert(curs, "MSSpectra", Bunch(forTarget=target.id, mzs=",".join([str(mz) for mz in scanMS2Annotated.labelledMz_list]),
                                                        ints=",".join([str(i) for i in scanMS2Annotated.labelledIntensity_list]),
                                                        annos=base64.b64encode(pickle.dumps([scanMS2Annotated.getAnnotationsForPeakInScan(labelled=True, index=i)[0] for i in range(len(scanMS2Annotated.nativeMz_list))])),
                                                        scanTime=scanMS2Labelled.retention_time, type="labelled_cleaned"),
                               writeFields=["forTarget", "mzs", "ints", "annos", "scanTime", "type"])

        for i in range(len(scanMS2Annotated.nativeMz_list)):
            mz=scanMS2Annotated.nativeMz_list[i]
            mzL=scanMS2Annotated.labelledMz_list[i]

            eicN, timesN, scanIdsN=chromatogram.getEIC(mz=mz,  ppm=self.matchingPPM, filterLine=nativeScanEventString, removeSingles=False, useMS1=False, useMS2=True)
            eicL, timesL, scanIdsL=chromatogram.getEIC(mz=mzL, ppm=self.matchingPPM, filterLine=labelledScanEventString, removeSingles=False, useMS1=False, useMS2=True)



            eicN, timesN, scanIdsN=cropEICs(eicN, timesN, scanIdsN, target.startRT*60, target.stopRT*60)
            eicL, timesL, scanIdsL=cropEICs(eicL, timesL, scanIdsL, target.startRT*60, target.stopRT*60)


            writeObjectAsSQLInsert(curs, "EICs", Bunch(forTarget=target.id,
                                                       forMZ=mz,
                                                       type="Native",
                                                       timesList=",".join([str(t) for t in timesN]),
                                                       intensityList=",".join([str(i) for i in eicN])),
                                   writeFields=["forTarget", "forMZ", "type", "timesList", "intensityList"])

            writeObjectAsSQLInsert(curs, "EICs", Bunch(forTarget=target.id,
                                                       forMZ=mz,
                                                       type="Labeled",
                                                       timesList=",".join([str(t) for t in timesL]),
                                                       intensityList=",".join([str(i) for i in eicL])),
                                   writeFields=["forTarget", "forMZ", "type", "timesList", "intensityList"])



        conn.commit()
        # endregion

        if self.saveAsPDF:
            self.saveResultsToPDF(target=target, curs=curs)
        if self.saveAsTSV:
            self.saveResultsToTSV(target=target, curs=curs)

    # endregion


