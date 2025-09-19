if __name__=="__main__":
    import MExtract

try:
    import rpy2
    import rpy2.robjects as ro
    r = ro.r
    R_AVAILABLE = True
except ImportError:
    print("Warning: rpy2 not available. MassSpecWavelet will use fallback methods.")
    R_AVAILABLE = False
    r = None

from utils import Bunch

errorIndex = 2


class MassSpecWavelet:
    def __init__(self, scriptLocation=None, scales=None, snrTh=None, minScans=None):
        self.scriptLocation=scriptLocation
        if self.scriptLocation==None:
            from utils import get_main_dir
            self.scriptLocation=get_main_dir() + "/chromPeakPicking/MassSpecWaveletIdentification.r"

        r('source(\"' + self.scriptLocation + '\")')

        self.scales=scales
        if self.scales==None:
            self.scales=[11,66]

        self.snrTh=snrTh
        if self.snrTh==None:
            self.snrTh=0.1

        self.minScans=minScans
        if self.minScans==None:
            self.minScans=3


    def cleanUp(self):
        r('rm(list=ls());gc();')

    def getPeaksFor(self, timesi, eici, startIndex=None, endIndex=None):

        snrTh=3

        # precheck
        if sum(eici) == 0:
            return []

        if startIndex is None:
            startIndex=0

        if endIndex is None:
            endIndex=len(eici)

        # crop EIC
        eici=eici[startIndex:endIndex]
        timesi=timesi[startIndex:endIndex]

        # add scales at begin and end of eic
        eicAdd = int(self.scales[1]) *2
        eic = [0 for u in range(0, eicAdd)]
        times = [0 for u in range(0, eicAdd)]
        eic.extend(eici)
        times.extend(timesi)

        eic.extend([0 for u in range(0, eicAdd)])
        mt=max(timesi)
        times.extend([mt for u in range(0, eicAdd)])

        # create R-vector of EIC
        eicRC = str([str(v) for v in eic]).replace('[', '').replace(']', '').replace('\'', '')
        timesRC= str([str(v) for v in times]).replace('[', '').replace(']', '').replace('\'', '')

        # call R-MassSpecWavelet
        try:   #eic, times, scales=c(5,33), snrTh=1, eicSmoothing="None", minCorr=0.85, testCorr=TRUE
            ret = r('getMajorPeaks(eic=c(' + eicRC + '), times=c(' + timesRC + '), scales=c(' + str(self.scales[0]) + ', ' + str(self.scales[1]) + '), snrTh=' + str(self.snrTh) + ')')
        except Exception as e:
            ret = []

        # Convert results from R-Objects to Python Objects/Arrays
        retl = []
        if len(ret) > 1:
            # Convert linear R-Vector (4 consecutive elements represent one peak) to Python Array
            for i in range(0, len(ret), 4):
                cur = []
                for j in range(0, 4):
                    cur.append(ret[i + j])

                peak=Bunch(peakIndex=int(cur[0]), peakScale=float(cur[1]), peakSNR=float(cur[2]),
                           peakArea=float(cur[3]), peakLeftFlank=float(cur[1]), peakRightFlank=float(cur[1]))
                retl.append(peak)

            # detect and remove very close-by peaks
            todel = []
            for f in range(0, len(retl)):
                for s in range(0, len(retl)):
                    if f < s:
                        if abs(retl[f].peakIndex - retl[s].peakIndex) < errorIndex:
                            todel.append(f)

                retl[f].peakIndex = retl[f].peakIndex - 1
                peakCenter = int(retl[f].peakIndex)
                count = 0
                for i in range(max(0, int(peakCenter - retl[f].peakScale)), min(len(eic) - 1, int(peakCenter + retl[f].peakScale))):
                    if eic[i] > 0:
                        count = count + 1
                if count < self.minScans:
                    todel.append(f)

            # remove very close-by peaks
            if len(todel) > 0:
                todel = [x for x in set(todel)]
                todel.sort()
                todel.reverse()
                for i in todel:
                    retl.pop(i)

            # calculate peak flanks (left and right side)
            for i1 in range(len(retl)):
                for i2 in range(len(retl)):
                    p1 = retl[i1]
                    p2 = retl[i2]
                    if p1.peakIndex < p2.peakIndex and (p1.peakIndex + p1.peakScale) > (p2.peakIndex - p2.peakScale):
                        overlap = (p1.peakIndex + p1.peakScale) - (p2.peakIndex - p2.peakScale)
                        if p1.peakScale > p2.peakScale:
                            p1.peakRightFlank = p1.peakRightFlank - overlap
                        else:
                            p2.peakLeftFlank = p2.peakLeftFlank - overlap

            # remove added ms scans from peak center
            # and add startIndex
            for f in retl:
                f.peakIndex = f.peakIndex - eicAdd + startIndex

        return retl



if __name__ == '__main__':

    import matplotlib
    import matplotlib.pyplot as plt


    import Chromatogram
    from utils import getLastTimeBefore, smoothDataSeries

    chromatogram = Chromatogram.Chromatogram()

    if True:
        chromatogram.parse_file("H:/180924_463_12C13CDON_Celegans_short_polar/12C13C_DON_onRemus_Exp27.mzXML")
        mz = 297.1331
        ppm = 5.
        cn = 15
        z = 1
        scales = [2, 5]
        dmz = 1.00335484



    scanEvents = chromatogram.getFilterLines()
    for s in scanEvents:
        print(s)

    scanEvent = sorted(scanEvents)[0]
    print("selected scan event:", scanEvent)

    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(111)

    eic, times, timesI, mzs  = chromatogram.getEIC(mz, ppm, scanEvent)
    eicL, times, timesI, mzsL = chromatogram.getEIC(mz+cn*dmz/z, ppm, scanEvent)

    #eic=[e/6/1.6 for e in eic]
    #eicL=[e/6/1.6 for e in eicL]

    smoothWin="gaussian"  ## "triangle", "gaussian", 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
    #eic  = smoothDataSeries(times, eic , windowLen=1, window=smoothWin)
    #eicL = smoothDataSeries(times, eicL, windowLen=1, window=smoothWin)


    startIndex=getLastTimeBefore(times, refTime=0*60.)
    endIndex=getLastTimeBefore(times, refTime=10000*60.)
    for i in range(getLastTimeBefore(times, refTime=600.4*60), len(eic)):
        eic[i]=0
    
    from utils import printObjectsAsTable
    print("Startindex", startIndex, "EndIndex", endIndex)


    CP=MassSpecWavelet("./MassSpecWaveletIdentification.R", scales=scales, minScans=1)
    ret = CP.getPeaksFor(times, eic, startIndex=startIndex, endIndex=endIndex)
    for peak in ret:
        peak.peakAtTime=times[peak.peakIndex] / 60.
    print("Native")
    printObjectsAsTable(ret, ["peakAtTime", "peakIndex", "peakScale", "peakArea", "peakLeftFlank", "peakIndex", "peakRightFlank", "peakSNR"])

    ret = CP.getPeaksFor(times, eicL, startIndex=startIndex, endIndex=endIndex)
    for peak in ret:
        peak.peakAtTime=times[peak.peakIndex] / 60.
    print("Labeled")
    printObjectsAsTable(ret, ["peakAtTime", "peakIndex", "peakScale", "peakArea", "peakLeftFlank", "peakIndex", "peakRightFlank", "peakSNR"])
    

    ax.plot([t / 60. for t in times], [e for e in eic])
    ax.plot([t / 60. for t in times], [-e for e in eicL])

    plt.xlabel("RT [minutes]")
    plt.ylabel("Intensity")
    plt.title("EIC")

    plt.show()