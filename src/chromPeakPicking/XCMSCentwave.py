if __name__ == "__main__":
    from . import MExtract
    import os

    os.environ["R_HOME"] = "C:/development/R-3.3.2"

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

from ..utils import Bunch


r = ro.r if R_AVAILABLE else None

errorIndex = 2


class XCMSCentwave:
    def __init__(self, scriptLocation):
        r('source("' + scriptLocation + '")')

    def cleanUp(self):
        r("rm(list=ls());gc();")

    def getPeaksFor(self, times, eic, scales=[11, 66]):
        # precheck
        if sum(eic) == 0:
            return []

        # create R-vector of EIC
        eicRC = str([str(v) for v in eic]).replace("[", "").replace("]", "").replace("'", "")
        timesRC = str([str(v) for v in times]).replace("[", "").replace("]", "").replace("'", "")

        # call R-MassSpecWavelet
        try:  # eic, times, scales=c(5,33), snrTh=1, eicSmoothing="None", minCorr=0.85, testCorr=TRUE
            ret = r("getCentwavePeaks(scanTime=c(" + timesRC + "), intensity=c(" + eicRC + "), peakwidth=c(" + str(scales[0]) + ", " + str(scales[1]) + "));")
        except Exception as e:
            ret = []

        # Convert results from R-Objects to Python Objects/Arrays
        retl = []
        if len(ret) > 1:
            for i in range(0, ret.nrow):
                peak = Bunch(
                    peakIndex=int(ret.rx(i + 1, 12)[0]) - 1,
                    peakScale=int(ret.rx(i + 1, 11)[0]),
                    peakSNR=ret.rx(i + 1, 7)[0],
                    peakArea=ret.rx(i + 1, 5)[0],
                    peakLeftFlank=int(ret.rx(i + 1, 13)[0]) - 1,
                    peakRightFlank=int(ret.rx(i + 1, 14)[0]) - 1,
                )
                retl.append(peak)

        return retl


if __name__ == "__main__":
    import matplotlib
    import matplotlib.pyplot as plt

    CP = XCMSCentwave("./XCMSCentwave.r")

    from . import Chromatogram
    from ..utils import getLastTimeBefore, smoothDataSeries

    chromatogram = Chromatogram.Chromatogram()

    if True:
        chromatogram.parse_file("H:/180924_463_12C13CDON_Celegans_short_polar/12C13C_DON_onRemus_Exp27.mzXML")
        mz = 297.1331
        ppm = 5.0
        cn = 15
        z = 1
        scales = [15, 90]
        dmz = 1.00335484

    scanEvents = chromatogram.getFilterLines()
    for s in scanEvents:
        print(s)

    scanEvent = sorted(scanEvents)[0]
    print("selected scan event:", scanEvent)

    fig = plt.figure(facecolor="white")
    ax = fig.add_subplot(111)

    eic, times, timesI, mzs = chromatogram.getEIC(mz, ppm, scanEvent)
    eicL, times, timesI, mzsL = chromatogram.getEIC(mz + cn * dmz / z, ppm, scanEvent)

    ax.plot(times, [e for e in eic], color="slategrey")
    ax.plot(times, [-e for e in eicL], color="slategrey")

    from ..utils import printObjectsAsTable

    ret = CP.getPeaksFor(times, eic, scales=scales)
    for peak in ret:
        peak.peakAtTime = times[peak.peakIndex]
        ax.plot([peak.peakAtTime, peak.peakAtTime], [0, max(eic)])
        ax.plot(
            times[peak.peakLeftFlank : peak.peakRightFlank],
            eic[peak.peakLeftFlank : peak.peakRightFlank],
            linewidth=3,
        )
    print("Native")
    printObjectsAsTable(
        ret,
        [
            "peakAtTime",
            "peakIndex",
            "peakScale",
            "peakArea",
            "peakLeftFlank",
            "peakIndex",
            "peakRightFlank",
            "peakSNR",
        ],
    )

    ret = CP.getPeaksFor(times, eicL, scales=scales)
    for peak in ret:
        peak.peakAtTime = times[peak.peakIndex]
        ax.plot([peak.peakAtTime, peak.peakAtTime], [-max(eicL), 0])
        ax.plot(
            times[peak.peakLeftFlank : peak.peakRightFlank],
            [-e for e in eicL[peak.peakLeftFlank : peak.peakRightFlank]],
            linewidth=3,
        )
    print("Labeled")
    printObjectsAsTable(
        ret,
        [
            "peakAtTime",
            "peakIndex",
            "peakScale",
            "peakArea",
            "peakLeftFlank",
            "peakIndex",
            "peakRightFlank",
            "peakSNR",
        ],
    )

    plt.xlabel("RT [minutes]")
    plt.ylabel("Intensity")
    plt.title("EIC")

    plt.show()
