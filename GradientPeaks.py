import numpy
from math import floor
from utils import Bunch




class Peak:
    def __init__(self, peakIndex, peakLeftFlank=-1, peakRightFlank=-1, leftInflection=-1, rightInflection=-1, peakArea=-1):
        self.peakIndex = peakIndex
        self.peakLeftFlank = peakLeftFlank
        self.peakRightFlank = peakRightFlank
        self.leftInflection = leftInflection
        self.rightInflection = rightInflection
        self.peakArea = peakArea

    def __str__(self):
        return "Center: %d leftFlank: %d rightFlank: %d Area:%.0f"%(self.peakIndex, self.peakLeftFlank, self.peakRightFlank, self.peakArea)


class GradientPeaks:
    def __init__(self, minInt=1000, minIntFlanks=10, minIncreaseRatio=.05, expTime=[5,45], minFlankToCenterRatio=1):
        self.minInt = minInt
        self.minIntFlanks=minIntFlanks
        self.minFlankToCenterRatio=minFlankToCenterRatio
        self.minIncreaseRatio = minIncreaseRatio
        self.expTime = expTime

    ### Find any local maximum in a series. Checking and verification of that peak will be performed later
    def getLocalMaxima(self, y, mInt=0):
        maxima=[]

        ## include left start of series if possible
        if y[0]>=mInt and y[0]>y[1]:
            maxima.append(0)

        ## check any signal except the start and ending
        for i in range(1, len(y)-1):
            if y[i]>=mInt and y[i-1]<y[i]>y[i+1]:
                maxima.append(i)

        ## include right ending of series if possible
        if y[len(y)-1]>=mInt and y[len(y)-2]<=y[len(y)-1]:
            maxima.append(len(y)-1)

        return maxima


    ### Add up all intensities between the peak's left and right flank
    def calculatePeakArea(self, x, y, peak):
        peak.peakArea=sum([y[i] for i in range(peak.peakIndex-peak.peakLeftFlank, peak.peakIndex+peak.peakRightFlank+1)])

        return True

    ### expand a local maxima by following its gradient to the bottom
    def expandPeak(self, peak, x, y):

        # find left border
        peak.peakLeftFlank=0
        goOn=True
        goOnInflection=True
        lastDelta=0
        while goOn and (peak.peakIndex-peak.peakLeftFlank-1)>0:
            rat=y[peak.peakIndex-peak.peakLeftFlank-1]/y[peak.peakIndex-peak.peakLeftFlank]
            delta=(y[peak.peakIndex-peak.peakLeftFlank]-y[peak.peakIndex-peak.peakLeftFlank-1])/(x[peak.peakIndex-peak.peakLeftFlank]-x[peak.peakIndex-peak.peakLeftFlank-1])

            if goOnInflection and delta<=lastDelta:
                peak.leftInflection=peak.peakLeftFlank
                goOnInflection=False
            else:
                lastDelta=delta


            if rat>self.minIncreaseRatio and rat<1. and y[peak.peakIndex-peak.peakLeftFlank-1]>=self.minIntFlanks:
                peak.peakLeftFlank+=1
            else:
                goOn=False

        #find right border
        peak.peakRightFlank=0
        goOn=True
        goOnInflection=True
        lastDelta=0
        while goOn and (peak.peakIndex+peak.peakRightFlank+1)<len(x):
            rat=y[peak.peakIndex+peak.peakRightFlank+1]/y[peak.peakIndex+peak.peakRightFlank]
            delta=(y[peak.peakIndex+peak.peakRightFlank]-y[peak.peakIndex+peak.peakRightFlank+1])/(x[peak.peakIndex+peak.peakRightFlank+1]-x[peak.peakIndex+peak.peakRightFlank])

            if goOnInflection and delta<=lastDelta:
                peak.rightInflection=peak.peakRightFlank
                goOnInflection=False
            else:
                lastDelta=delta


            if rat>self.minIncreaseRatio and rat<1. and y[peak.peakIndex+peak.peakRightFlank+1]>=self.minIntFlanks:
                peak.peakRightFlank+=1
            else:
                goOn=False

        peakWidthValid=self.expTime[0]<=peak.peakRightFlank<=self.expTime[1] or self.expTime[0]<=peak.peakLeftFlank<=self.expTime[1]

        return peakWidthValid

    def calcPeakProperties(self, x, y, peak):
        peak.leftDelta=(y[peak.peakIndex]-y[peak.peakIndex-peak.peakLeftFlank])/(x[peak.peakIndex]-x[peak.peakIndex-peak.peakLeftFlank])
        peak.rightDelta=(y[peak.peakIndex]-y[peak.peakIndex+peak.peakRightFlank])/(x[peak.peakIndex+peak.peakRightFlank]-x[peak.peakIndex])
        peak.leftInflectionDelta=(y[peak.peakIndex]-y[peak.peakIndex-peak.leftInflection])/(x[peak.peakIndex]-x[peak.peakIndex-peak.leftInflection])
        peak.rightInflectionDelta=(y[peak.peakIndex]-y[peak.peakIndex+peak.rightInflection])/(x[peak.peakIndex+peak.rightInflection]-x[peak.peakIndex])

        return True




    def findPeaks(self, x, y):
        localMaxs=self.getLocalMaxima(y, self.minInt)

        # check local maxima for peaks
        peaks = []
        for localMax in localMaxs:
            peak = Peak(peakIndex=localMax)
            # calculate peak borders and expand peak, check if within expected chromatographic width
            usePeak = self.expandPeak(peak, x, y)
            usePeak = usePeak and self.calculatePeakArea(x, y, peak)
            usePeak = usePeak and self.calcPeakProperties(x, y, peak)

            if usePeak:
                peaks.append(peak)

        return peaks

    ######################################################################################
    ### NOTE: Parameters only kept for compatibility. Will be removed in later version ###
    ######################################################################################
    def getPeaksFor(self, times, eic, scales=None, snrTh=0.1, startIndex=None, endIndex=None):

        peaks=self.findPeaks(times, eic)

        ret=[]
        for peak in peaks:
            ret.append(Bunch(peakIndex=peak.peakIndex, peakScale=(peak.peakRightFlank+peak.peakLeftFlank)/2, peakSNR=100, peakArea=peak.peakArea, peakLeftFlank=peak.peakLeftFlank, peakRightFlank=peak.peakRightFlank))


        return ret

if __name__=="__main__":
    from Chromatogram import Chromatogram
    from utils import smoothDataSeries

    t = Chromatogram()

    if True:
        mzXMLFile="F:/160112_238_posneg_Labelled_wheat_experiment/160112_posneg_236_12C13C_fullyLab_Remus_DON_1.mzXML"
        mz = 375.14439
        ppm = 8.
        cn=20

    if False:
        mzXMLFile="E:/Cambridge/Jan_2015/wetransfer-cd97f3/QE_PR493_PV_051014-AIF-02.mzXML"
        mzXMLFile="E:/Cambridge/Apr_2015/wetransfer-2d73b0/RNA-elegans-large_150328005031.mzXML"
        mzXMLFile="E:/Cambridge/Apr_2015/wetransfer-2d73b0/RNA-elegans-small.mzXML"
        mzXMLFile="E:/Cambridge/Apr_2015/wetransfer-2d73b0/RNA-elegans-large13C.mzXML"
        mz=268.103685
        mz=365.166
        cn=5
        ppm=5.


    t.parse_file(mzXMLFile)

    filterLine=[i for i in t.getFilterLines(includeMS2=False)][0]
    print filterLine
    eicN, times, scanIds = t.getEIC(mz, ppm, filterLine=filterLine)
    eicL, times, scanIds = t.getEIC(mz+cn*1.00335, ppm, filterLine)
    timesMin=[t/60. for t in times]

    import matplotlib.pyplot as plt
    plt.plot(timesMin, eicN)
    plt.plot(timesMin, [-e for e in eicL])

    eicN = smoothDataSeries(times, eicN, windowLen=2, window="Gaussian")
    eicL = smoothDataSeries(times, eicL, windowLen=2, window="Gaussian")

    eic=eicN



    plt.plot(timesMin, eicN)
    plt.plot(timesMin, [-e for e in eicL])

    gp=GradientPeaks(minInt=100, minIntFlanks=10, minIncreaseRatio=.05, expTime=[3, 45])
    #gp=GradientPeaks(minInt=50000, minIntFlanks=5000, minIncreaseRatio=.05, expTime=[5,250], minFlankToCenterRatio=1)
    peaks=gp.findPeaks(times, eicN)
    peaksL=gp.findPeaks(times, eicL)

    from utils import printObjectsAsTable

    "GradientPeaks"
    for peak in peaks: peak.minRT=timesMin[peak.peakIndex]
    for peak in peaksL: peak.minRT=timesMin[peak.peakIndex]
    printObjectsAsTable(peaks, attrs=["minRT", "peakLeftFlank", "leftInflection", "peakIndex", "rightInflection", "peakRightFlank", "peakArea", "leftDelta", "leftInflectionDelta", "rightInflectionDelta", "rightDelta"])
    printObjectsAsTable(peaksL, attrs=["minRT", "peakLeftFlank", "leftInflection", "peakIndex", "rightInflection", "peakRightFlank", "peakArea", "leftDelta", "leftInflectionDelta", "rightInflectionDelta", "rightDelta"])

    for peak in peaks:
        lind=peak.peakIndex-peak.peakLeftFlank
        rind=peak.peakIndex+peak.peakRightFlank
        #print "startMin:",timesMin[lind], "startInflection:",timesMin[peak.peakIndex-peak.leftInflection], "centerMin:",timesMin[peak.peakIndex], "endInflection:",timesMin[peak.peakIndex+peak.rightInflection], "endMin:",timesMin[rind], "peaksCovered:",peak.peakRightFlank+peak.peakLeftFlank, peak.peakArea, ((eic[peak.peakIndex]-eic[peak.peakIndex+peak.peakRightFlank])/eic[peak.peakIndex+peak.peakRightFlank])

        plt.plot([timesMin[lind],timesMin[peak.peakIndex-peak.leftInflection],timesMin[peak.peakIndex],timesMin[peak.peakIndex+peak.rightInflection],timesMin[rind]],
                 [eic[lind],eic[peak.peakIndex-peak.leftInflection],eic[peak.peakIndex],eic[peak.peakIndex+peak.rightInflection],eic[rind]], 'yD')
        plt.fill_between(timesMin[lind:(rind+1)], 0, eic[lind:(rind+1)], color="red")
        plt.plot([timesMin[lind], timesMin[peak.peakIndex-peak.leftInflection], timesMin[peak.peakIndex], timesMin[peak.peakIndex+peak.rightInflection], timesMin[rind]],
                 [eic[lind], eic[peak.peakIndex-peak.leftInflection], eic[peak.peakIndex], eic[peak.peakIndex+peak.rightInflection], eic[rind]], "black")
    eicL=[-e for e in eicL]
    for peak in peaksL:
        lind=peak.peakIndex-peak.peakLeftFlank
        rind=peak.peakIndex+peak.peakRightFlank
        #print "startMin:",timesMin[lind], "startInflection:",timesMin[peak.peakIndex-peak.leftInflection], "centerMin:",timesMin[peak.peakIndex], "endInflection:",timesMin[peak.peakIndex+peak.rightInflection], "endMin:",timesMin[rind], "peaksCovered:",peak.peakRightFlank+peak.peakLeftFlank, peak.peakArea, ((eic[peak.peakIndex]-eic[peak.peakIndex+peak.peakRightFlank])/eic[peak.peakIndex+peak.peakRightFlank])

        plt.plot([timesMin[lind],timesMin[peak.peakIndex-peak.leftInflection],timesMin[peak.peakIndex],timesMin[peak.peakIndex+peak.rightInflection],timesMin[rind]],
                 [eicL[lind],eicL[peak.peakIndex-peak.leftInflection],eicL[peak.peakIndex],eicL[peak.peakIndex+peak.rightInflection],eicL[rind]], 'yD')
        plt.fill_between(timesMin[lind:(rind+1)], 0, eicL[lind:(rind+1)], color="red")
        plt.plot([timesMin[lind], timesMin[peak.peakIndex-peak.leftInflection], timesMin[peak.peakIndex], timesMin[peak.peakIndex+peak.rightInflection], timesMin[rind]],
                 [eicL[lind], eicL[peak.peakIndex-peak.leftInflection], eicL[peak.peakIndex], eicL[peak.peakIndex+peak.rightInflection], eicL[rind]], "black")

        #plt.plot(timesMin[lind:(rind+1)], eic[lind:(rind+1)]/max(eic[lind:(rind+1)]), color="black")

        #pdf=mlab.normpdf(np.array(timesMin[lind:rind]), timesMin[peak.peakIndex], (rind-lind)/500.)
        #pdf=pdf/max(pdf)
        #plt.plot(timesMin[lind:rind], pdf, color="yellow")

    import MExtract
    from MassSpecWavelet import MassSpecWavelet

    cp=MassSpecWavelet("./MassSpecWaveletIdentification.R")
    print "------------"

    peaks=cp.getPeaksFor(times, eic, scales=[3,19])

    for peak in peaks:
        peak.peakAtTime=times[peak.peakIndex] / 60.

    from utils import printObjectsAsTable
    print "MassSpecWavelet"
    printObjectsAsTable(peaks, ["peakAtTime", "peakLeftFlank", "peakIndex", "peakRightFlank"])
    for peak in peaks:
        lind=int(peak.peakIndex-peak.peakLeftFlank)
        rind=int(peak.peakIndex+peak.peakRightFlank)
        #print timesMin[lind], timesMin[peak[0]], timesMin[rind], peak[4]+peak[5]
        plt.fill_between(timesMin[lind:(rind+1)], 0, eic[lind:(rind+1)], color="green", alpha=.5)

    plt.show()



















