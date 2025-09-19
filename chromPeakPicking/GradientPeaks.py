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
    def __init__(self, minInt=1000, minIntFlanks=10, minIncreaseRatio=.05, expTime=[5,45], minDelta=1000, minInflectionDelta=2000, minBaseLineRatio=0.5):
        self.minInt = minInt
        self.minIntFlanks=minIntFlanks
        self.minIncreaseRatio = minIncreaseRatio
        self.expTime = expTime
        self.minDelta=minDelta
        self.minInflectionDelta=minInflectionDelta
        self.minBaseLineRatio=minBaseLineRatio

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
        ##print("peak at ", peak.peakIndex, x[peak.peakIndex])

        # find left border
        peak.peakLeftFlank=0
        goOn=True
        goOnInflection=True
        lastDelta=0
        while goOn and (peak.peakIndex-peak.peakLeftFlank-1)>0:
            rat=y[peak.peakIndex-peak.peakLeftFlank-1]/y[peak.peakIndex-peak.peakLeftFlank]
            delta=(y[peak.peakIndex-peak.peakLeftFlank]-y[peak.peakIndex-peak.peakLeftFlank-1])/(x[peak.peakIndex-peak.peakLeftFlank]-x[peak.peakIndex-peak.peakLeftFlank-1])
            #print("left", y[peak.peakIndex-peak.peakLeftFlank-1], y[peak.peakIndex-peak.peakLeftFlank], rat, delta, goOnInflection, goOn)

            if goOnInflection and delta<=lastDelta:
                peak.leftInflection=peak.peakLeftFlank
                goOnInflection=False
            else:
                lastDelta=delta


            if rat>self.minIncreaseRatio and rat<1. and y[peak.peakIndex-peak.peakLeftFlank-1]>=self.minIntFlanks:
                peak.peakLeftFlank+=1
            elif not goOnInflection:
                goOn=False

        #find right border
        peak.peakRightFlank=0
        goOn=True
        goOnInflection=True
        lastDelta=0
        while goOn and (peak.peakIndex+peak.peakRightFlank+1)<len(x):
            rat=y[peak.peakIndex+peak.peakRightFlank+1]/y[peak.peakIndex+peak.peakRightFlank]
            delta=(y[peak.peakIndex+peak.peakRightFlank]-y[peak.peakIndex+peak.peakRightFlank+1])/(x[peak.peakIndex+peak.peakRightFlank+1]-x[peak.peakIndex+peak.peakRightFlank])
            #print("right", y[peak.peakIndex-peak.peakLeftFlank-1], y[peak.peakIndex-peak.peakLeftFlank], rat, delta, goOnInflection, goOn)

            if goOnInflection and delta<=lastDelta:
                peak.rightInflection=peak.peakRightFlank
                goOnInflection=False
            else:
                lastDelta=delta


            if rat>self.minIncreaseRatio and rat<1. and y[peak.peakIndex+peak.peakRightFlank+1]>=self.minIntFlanks:
                peak.peakRightFlank+=1
            elif not goOnInflection:
                goOn=False

        peakWidthValid=self.expTime[0]<=peak.peakRightFlank<=self.expTime[1] or self.expTime[0]<=peak.peakLeftFlank<=self.expTime[1]

        return peakWidthValid

    def calcPeakProperties(self, x, y, peak):
        peak.leftDelta=(y[peak.peakIndex]-y[peak.peakIndex-peak.peakLeftFlank])
        peak.rightDelta=(y[peak.peakIndex]-y[peak.peakIndex+peak.peakRightFlank])
        peak.leftInflectionDelta=(y[peak.peakIndex]-y[peak.peakIndex-peak.leftInflection])
        peak.rightInflectionDelta=(y[peak.peakIndex]-y[peak.peakIndex+peak.rightInflection])

        peak.baseLineRatio=(y[peak.peakIndex]-(y[peak.peakIndex-peak.peakLeftFlank]+y[peak.peakIndex+peak.peakRightFlank])/2.)/y[peak.peakIndex]

        return (peak.leftDelta>=self.minDelta and peak.rightDelta>=self.minDelta) and (peak.leftInflectionDelta>=self.minInflectionDelta and peak.rightInflectionDelta>=self.minInflectionDelta) and peak.baseLineRatio>=self.minBaseLineRatio




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


    from utilities.Stopwatch import Stopwatch
    sw=Stopwatch()
    durationMessages=[]

    t = Chromatogram()

    if False:
        mzXMLFile="F:/160112_238_posneg_Labelled_wheat_experiment/160112_posneg_236_12C13C_fullyLab_Remus_DON_1.mzXML"
        mz = 375.14439
        ppm = 8.
        cn=20
        dmz=1.00335
        z=1

    if False:
        mzXMLFile="E:/Cambridge/Jan_2015/wetransfer-cd97f3/QE_PR493_PV_051014-AIF-02.mzXML"
        mzXMLFile="E:/Cambridge/Apr_2015/wetransfer-2d73b0/RNA-elegans-large_150328005031.mzXML"
        mzXMLFile="E:/Cambridge/Apr_2015/wetransfer-2d73b0/RNA-elegans-small.mzXML"
        mzXMLFile="E:/Cambridge/Apr_2015/wetransfer-2d73b0/RNA-elegans-large13C.mzXML"
        mz=268.103685
        mz=365.166
        cn=5
        ppm=5.
        dmz=1.00335
        z=1

    if False:
        mzXMLFile="F:/BenediktWarth/160303_final raw data/neg/GEN_24h_5uM_lysate_2_neg_066.mzXML"
        mz = 175.9574#349.002275
        ppm = 5.
        cn = 3
        z = 1
        dmz=1.00628

    if False:
        mzXMLFile="F:/Thermo_OribtrapTest/OrbitrapHF_Swiss/Vial1_pos_60k_06.mzXML"
        mz = 359.131328#527.155036#354.17606#359.131328
        ppm = 5.
        cn = 14#27#14
        z = 1
        dmz = 1.00335
    if False:
        mzXMLFile="F:/HrKoch/Dataset/151215_pos_236_13C_Trp/151215_pos_236_fullyLab_Remus.mzXML"

        mz = 295.22656
        ppm = 5.
        cn = 18
        z = 1
        scales = [5,19]
        dmz = 1.00335

    if True:
        mzXMLFile="F:/Waterhouse/Quinones/mzXML/a12131_neg.mzXML"

        mz = 1006.34489#912.31805
        ppm = 50.
        cn = 12
        z = 1
        scales = [5,19]
        dmz = 1.00335
        gp=GradientPeaks(minInt=1000, minIntFlanks=100, minIncreaseRatio=.5, minDelta=1, expTime=[3, 150])


    if False:
        mzXMLFile="F:/Bayer_BernhardKluger/HSS_pos/pos_12C13C_L19_Mock_Rep1_28.mzXML"

        mz=519.11316
        ppm=5.
        cn=24
        z=1
        scales=[3,6]
        dmz=1.00335
        gp=GradientPeaks(minInt=1000, minIntFlanks=10, minIncreaseRatio=.05, minDelta=10000, expTime=[5, 150]) ## Bernhard


    if False:
        mzXMLFile="F:/Sciex/IFA_Tulin_mzXML_A/IFA_Tulln_neg_CM_stem_12C13C_1_01_CM_stem_12C13C_Vial_1.mzXML"

        mz=583.21851
        ppm=5.
        cn=31
        z=1
        scales=[3,6]
        dmz=1.00335
        gp=GradientPeaks(minInt=50, minIntFlanks=10, minIncreaseRatio=.05, expTime=[15, 150], minDelta=10000, minInflectionDelta=2)


    if False:                            ##Bruker
        mzXMLFile="F:/10er Mix 100 ppb foc_RA3_01_1138.d/10er_Mix_100_ppb_foc_RA3_01_1138.mzXML"

        mz=716.45728
        ppm=5.
        cn=31
        z=1
        scales=[3,6]
        dmz=1.00335
        gp=GradientPeaks(minInt=50, minIntFlanks=10, minIncreaseRatio=.05, expTime=[15, 150], minDelta=10000, minInflectionDelta=2)

    sw.start()
    t.parse_file(mzXMLFile)
    sw.stop()
    durationMessages.append("Importing the file took %.1f seconds (1 time)"%sw.getDurationInSeconds())

    filterLine=[i for i in t.getFilterLines(includeMS2=False)][0]

    sw.start()
    for i in range(100):
        eicN_raw, times, scanIds, mzs = t.getEIC(mz, ppm, filterLine=filterLine)
        eicL_raw, times, scanIds, mzs = t.getEIC(mz+cn*dmz/z, ppm, filterLine=filterLine)
    sw.stop()
    durationMessages.append("Calculating the EICs took %.1f seconds (100 times)"%sw.getDurationInSeconds())




    timesMin=[t/60. for t in times]


    import matplotlib.pyplot as plt
    plt.plot(timesMin, eicN_raw)
    plt.plot(timesMin, [-e for e in eicL_raw])

    smoothWin="savitzkygolay"  ## "gaussian", "flat", "triangle", "gaussian", "flat", "savitzkygolay"

    sw.start()
    for i in range(1000):
        eicN = smoothDataSeries(times, eicN_raw, windowLen=9, window=smoothWin, polynom=3)
        eicL = smoothDataSeries(times, eicL_raw, windowLen=9, window=smoothWin, polynom=3)
    sw.stop()
    durationMessages.append("Smoothing the EICs took %.1f seconds (1000 times)"%sw.getDurationInSeconds())

    eic=eicN



    plt.plot(timesMin, eicN)
    plt.plot(timesMin, [-e for e in eicL])



    sw.start()
    for i in range(1000):
        peaks=gp.findPeaks(times, eicN)
        peaksL=gp.findPeaks(times, eicL)
    sw.stop()
    durationMessages.append("Detecting chromatographic peaks with GradientDescend took %.1f seconds (1000 times)"%sw.getDurationInSeconds())

    from utils import printObjectsAsTable

    "GradientPeaks"
    for peak in peaks: peak.minRT=timesMin[peak.peakIndex]
    for peak in peaksL: peak.minRT=timesMin[peak.peakIndex]
    printObjectsAsTable(peaks , attrs=["minRT", "peakLeftFlank", "leftInflection", "peakIndex", "rightInflection", "peakRightFlank", "peakArea", "leftDelta", "leftInflectionDelta", "rightInflectionDelta", "rightDelta", "baseLineRatio"])
    printObjectsAsTable(peaksL, attrs=["minRT", "peakLeftFlank", "leftInflection", "peakIndex", "rightInflection", "peakRightFlank", "peakArea", "leftDelta", "leftInflectionDelta", "rightInflectionDelta", "rightDelta", "baseLineRatio"])

    for peak in peaks:
        lind=peak.peakIndex-peak.peakLeftFlank
        rind=peak.peakIndex+peak.peakRightFlank

        plt.plot([timesMin[lind],timesMin[peak.peakIndex-peak.leftInflection],timesMin[peak.peakIndex],timesMin[peak.peakIndex+peak.rightInflection],timesMin[rind]],
                 [eic[lind],eic[peak.peakIndex-peak.leftInflection],eic[peak.peakIndex],eic[peak.peakIndex+peak.rightInflection],eic[rind]], 'yD')
        plt.fill_between(timesMin[lind:(rind+1)], 0, eic[lind:(rind+1)], color="red", alpha=0.25)
        plt.plot([timesMin[lind], timesMin[peak.peakIndex-peak.leftInflection], timesMin[peak.peakIndex], timesMin[peak.peakIndex+peak.rightInflection], timesMin[rind]],
                 [eic[lind], eic[peak.peakIndex-peak.leftInflection], eic[peak.peakIndex], eic[peak.peakIndex+peak.rightInflection], eic[rind]], "black")
    eicL=[-e for e in eicL]
    for peak in peaksL:
        lind=peak.peakIndex-peak.peakLeftFlank
        rind=peak.peakIndex+peak.peakRightFlank

        plt.plot([timesMin[lind],timesMin[peak.peakIndex-peak.leftInflection],timesMin[peak.peakIndex],timesMin[peak.peakIndex+peak.rightInflection],timesMin[rind]],
                 [eicL[lind],eicL[peak.peakIndex-peak.leftInflection],eicL[peak.peakIndex],eicL[peak.peakIndex+peak.rightInflection],eicL[rind]], 'yD')
        plt.fill_between(timesMin[lind:(rind+1)], 0, eicL[lind:(rind+1)], color="red", alpha=0.25)
        plt.plot([timesMin[lind], timesMin[peak.peakIndex-peak.leftInflection], timesMin[peak.peakIndex], timesMin[peak.peakIndex+peak.rightInflection], timesMin[rind]],
                 [eicL[lind], eicL[peak.peakIndex-peak.leftInflection], eicL[peak.peakIndex], eicL[peak.peakIndex+peak.rightInflection], eicL[rind]], "black")



    import MExtract
    from MassSpecWavelet import MassSpecWavelet

    cp=MassSpecWavelet("./MassSpecWaveletIdentification.R")
    print("------------")

    sw.start()
    for i in range(10):
        peaks=cp.getPeaksFor(times, eicN_raw, scales=[5,19])
        peaksL=cp.getPeaksFor(times, eicL_raw, scales=[5,19])
    sw.stop()
    durationMessages.append("Detecting chromatographic peaks with MassSpecWavelet took %.1f seconds (10 times)"%sw.getDurationInSeconds())

    for peak in peaks:
        peak.peakAtTime=times[peak.peakIndex] / 60.
    for peak in peaksL:
        peak.peakAtTime=times[peak.peakIndex] / 60.

    from utils import printObjectsAsTable
    print("MassSpecWavelet")
    printObjectsAsTable(peaks, ["peakAtTime", "peakLeftFlank", "peakIndex", "peakRightFlank"])
    printObjectsAsTable(peaksL, ["peakAtTime", "peakLeftFlank", "peakIndex", "peakRightFlank"])
    for peak in peaks:
        lind=int(peak.peakIndex-peak.peakLeftFlank)
        rind=int(peak.peakIndex+peak.peakRightFlank)
        plt.fill_between(timesMin[lind:(rind+1)], 0, eicN[lind:(rind+1)], color="green", alpha=.15)
    for peak in peaksL:
        lind=int(peak.peakIndex-peak.peakLeftFlank)
        rind=int(peak.peakIndex+peak.peakRightFlank)
        plt.fill_between(timesMin[lind:(rind+1)], 0, [e for e in eicL[lind:(rind+1)]], color="green", alpha=.15)





    print("\n".join(durationMessages))


    plt.show()












