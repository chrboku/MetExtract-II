import rpy2
import rpy2.robjects as ro

from utils import Bunch


r = ro.r

errorIndex = 2


class MassSpecWavelet:
    def __init__(self, scriptLocation):
        r('source(\"' + scriptLocation + '\")')

    def cleanUp(self):
        r('rm(list=ls());gc();')

    def getPeaksFor(self, times, eici, scales=[11, 66], snrTh=0.1, eicSmoothing="None", minScans=3, startIndex=None, endIndex=None):

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
        #times=times[startIndex:endIndex]

        # add scales at begin and end of eic
        eicAdd = int(scales[1]) *2
        eic = [0 for u in range(0, eicAdd)]
        eic.extend(eici)
        eic.extend([0 for u in range(0, eicAdd)])

        # create R-vector of EIC
        eicRC = str([str(v) for v in eic]).replace('[', '').replace(']', '').replace('\'', '')

        # call R-MassSpecWavelet
        try:
            ret = r('getMajorPeaks(c(' + eicRC + '), c(' + str(scales[0]) + ', ' + str(scales[1]) + '), snrTh=' + str(
                snrTh) + ', eicSmoothing=\"' + str(eicSmoothing) + '\")')
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
                if count < minScans:
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

    CP = MassSpecWavelet("./MassSpecWaveletIdentification.r")

    from Chromatogram import Chromatogram
    from utils import getLastTimeBefore


    chromatogram = Chromatogram()
    chromatogram.parse_file("/Volumes/Storage/IFA/implTest/LTQ_Orbitrap_XL/Remus_DON_1_Base.mzXML")
    mz = 297.132
    ppm = 5.
    scales = [3, 19]

    scanEvents = chromatogram.getFilterLines()
    for s in scanEvents:
        print s
    scanEvent = sorted(scanEvents)[0]
    print "selected scan event:", scanEvent

    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(111)

    eic, times, timesI = chromatogram.getEIC(mz, ppm, scanEvent)
    startIndex=getLastTimeBefore(times, refTime=6.5*60)
    endIndex=getLastTimeBefore(times, refTime=13.5*60)
    startIndex=258
    endIndex=643
    
    from utils import printObjectsAsTable
    ret = CP.getPeaksFor(times, eic, scales=scales, startIndex=startIndex, endIndex=endIndex)
    for peak in ret:
        peak.peakAtTime=times[peak.peakIndex] / 60.
    printObjectsAsTable(ret, ["peakAtTime", "peakArea", "peakLeftFlank", "peakIndex", "peakRightFlank"])
    

    ax.plot([t / 60. for t in times], [e for e in eic])

    plt.xlabel("RT [minutes]")
    plt.ylabel("Intensity")
    plt.title("EIC")

    plt.show()