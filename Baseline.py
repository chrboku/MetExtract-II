
if __name__=="__main__":
    import MExtract

try:
    import rpy2
    R_AVAILABLE = True
except ImportError:
    R_AVAILABLE = False
    rpy2 = None
try:
    import rpy2.robjects as ro
    R_AVAILABLE = True
except ImportError:
    R_AVAILABLE = False
    ro = None

from utils import Bunch


r = ro.r if R_AVAILABLE else None

class Baseline:

    def __init__(self):
        r("suppressWarnings(suppressMessages(library(baseline))); getBaseline<-function(eic){  return(baseline(eic, method='medianWindow', hwm=20)@baseline[1,])  };")

    def getBaseline(self, eic, times):

        eicR="c(" + ",".join(str(e) for e in eic) + ")"

        try:
            ret = r("getBaseline(t("+eicR+"));")

            ret=[float(i) for i in ret]
            return ret

        except Exception as e:
            raise Exception("Baseline could not be calculated successfully")





if __name__=="__main__":
    import matplotlib
    import matplotlib.pyplot as plt

    BL = Baseline()

    import Chromatogram
    from utils import getLastTimeBefore, smoothDataSeries
    from utils import printObjectsAsTable

    import chromPeakPicking.MassSpecWavelet

    CP = chromPeakPicking.MassSpecWavelet.MassSpecWavelet("chromPeakPicking/MassSpecWaveletIdentification.r")


    chromatogram = Chromatogram.Chromatogram()

    if True:
        chromatogram.parse_file("C:/PyMetExtract/implTest/D-GEN_Warth/GEN_24h_5uM_supern_3_neg_079.mzXML")
        mz = 370.90112
        #mz = 283.06069
        ppm = 10.
        cn = 4
        z = 1
        scales = [2, 7]
        dmz = 1.00628

    scanEvents = chromatogram.getFilterLines()
    for s in scanEvents:
        print(s)

    scanEvent = sorted(scanEvents)[0]
    print("selected scan event:", scanEvent)

    fig = plt.figure(facecolor='white')
    ax = fig.add_subplot(111)

    eic, times, timesI, mzs = chromatogram.getEIC(mz, ppm, scanEvent)
    eicL, times, timesI, mzsL = chromatogram.getEIC(mz + cn * dmz / z, ppm, scanEvent)


    ax.plot([t / 60. for t in times], [e for e in eic])

    eicBL=BL.getBaseline(eic, times)

    ax.plot([t / 60. for t in times], [e for e in eicBL])

    ax.plot([t / 60. for t in times], [-max(0, eic[i]-eicBL[i]) for i in range(len(eic))])

    ret = CP.getPeaksFor(times, eic, scales=scales, minScans=1)
    for peak in ret:
        peak.peakAtTime = times[peak.peakIndex] / 60.
    print("Native")
    printObjectsAsTable(ret, ["peakAtTime", "peakIndex", "peakScale", "peakArea", "peakLeftFlank", "peakIndex",
                              "peakRightFlank", "peakSNR"])


    ret = CP.getPeaksFor(times, [max(0, eic[i]-eicBL[i]) for i in range(len(eic))], scales=scales, minScans=1)
    for peak in ret:
        peak.peakAtTime = times[peak.peakIndex] / 60.
    print("Native without baseline")
    printObjectsAsTable(ret, ["peakAtTime", "peakIndex", "peakScale", "peakArea", "peakLeftFlank", "peakIndex",
                              "peakRightFlank", "peakSNR"])

    plt.show()