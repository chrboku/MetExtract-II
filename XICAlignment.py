try:
    import rpy2.robjects as ro
    R_AVAILABLE = True
except ImportError:
    R_AVAILABLE = False
    ro = None
r = ro.r if R_AVAILABLE else None

from utils import mapArrayToRefTimes
from copy import deepcopy


# HELPER METHOD for adding a constant part before and after the actual EIC
def preAppendEICs(eics, scanTimes, pretend=25, scanDuration=1, mIntMult=1, addFront=True, addBack=False):
    assert len(eics)==len(scanTimes)
    eics=deepcopy(eics)
    scanTimes=deepcopy(scanTimes)

    if pretend==0:
        return eics, scanTimes

    for j in range(len(eics)):
        pInt=max(eics[j])*mIntMult

        scanTimes[j]=[scanTimes[j][i]+scanDuration*pretend for i in range(len(scanTimes[j]))]
        mScanTime=max(scanTimes[j])

        for i in range(pretend):
            if addFront:
                eics[j].insert(0, pInt)
                scanTimes[j].insert(i, i*scanDuration)
            if addBack:
                eics[j].append(pInt)
                scanTimes[j].append(mScanTime+(1+i)*scanDuration)


    return eics, scanTimes

# HELPER METHOD for removing a constant part before and after the actual EIC.
# Used im combination with preAppendEICs
def undoPreAppendEICs(eics, refTimes, pretend=25, scanDuration=1, removeFront=True, removeBack=False):

    eics=deepcopy(eics)
    refTimes=deepcopy(refTimes)

    if pretend==0:
        return eics, refTimes

    remTime=scanDuration*(pretend)
    mStart=min(enumerate(refTimes), key=lambda x: abs(x[1]-remTime))[0]

    if removeFront:
        refTimes=refTimes[mStart:]
        refTimes=[refTimes[i]-remTime for i in range(len(refTimes))]

    if removeBack:
        refTimes=refTimes[:-mStart]


    for j in range(len(eics)):
        if removeFront:
            eics[j]=eics[j][mStart:]
        if removeBack:
            eics[j]=eics[j][:-mStart]

    return eics, refTimes


# Class for chromatographic alignment using PTW (http://cran.fhcrc.org/web/packages/ptw/index.html) and
# bracketing of different feature pairs
class XICAlignment:
    def __init__(self, scriptLocation):
        r('source(\"' + scriptLocation + '\")')

    # aligns the EICs and matches the chromatographic peaks detected earlier
    def alignXIC(self, eics, peakss, scantimes, align=True, nPolynom=3, maxTimeDiff=0.36 * 60, pretend=100, scanDuration=1):
        #assert (len(eics) == len(peakss) == len(scantimes))


        nPolynom = max(nPolynom, 1)

        if len(peakss) == 0:
            return []

        if len(peakss) == 1:
            ret = []
            for i in range(len(peakss[0])):
                ret.append([[i, peakss[0][i].NPeakCenterMin]])
            return ret

        refTimes = ""
        eicsR = ""
        peakssR = ""
        if align:
            # add a constant part before and after the actual EIC as an anchor for the alignment
            eics, scantimesEnlarged= preAppendEICs(eics, scantimes, pretend=pretend, scanDuration=scanDuration)

            maxTime = max((max(scanTime) for scanTime in scantimes))
            refTimes = [i*scanDuration for i in range(int(maxTime/scanDuration))]

            for i in range(len(eics)):
                eics[i] = mapArrayToRefTimes(eics[i], scantimesEnlarged[i], refTimes)
                # peak translation to ref Time
                for j in range(len(peakss[i])):
                    peakTime = scantimes[i][peakss[i][j].NPeakCenter]+pretend*scanDuration
                    minindex, minvalue = min(enumerate((abs(z - peakTime) for z in refTimes)), key=lambda x:x[1])
                    peakss[i][j].NPeakCenter = refTimes[minindex]


            # convert the EICs to R-function parameters
            maxLen = max([len(eic) for eic in eics])
            append = False
            for eic in eics:
                if append:
                    eicsR = eicsR + ","
                else:
                    append = True
                eicsR = eicsR + ",".join((str(d) for d in eic))
                if len(eic) < maxLen:
                    eicsR = eicsR + "," + ",".join((str(0) for p in range(len(eic), maxLen)))

        # convert the detected chromatographic peaks to R-function parameters
        maxPeaks = max([len(peaks) for peaks in peakss])
        append = False
        for peaks in peakss:
            if append:
                peakssR = peakssR + ","
            else:
                append = True
            peakssR = peakssR + ",".join((str(peak.NPeakCenter) if align else str(peak.NPeakCenterMin) for peak in peaks))
            if len(peaks) < maxPeaks:
                peakssR = peakssR + "," + ",".join((str(0) for p in range(len(peaks), maxPeaks)))

        # align the EICs and bracket the chromatographic peaks
        ret = r('alignPeaks(c(' + eicsR + '), c(' + ",".join([str(f) for f in refTimes]) + '), c(' + peakssR + '), ' + str(len(peakss)) + ', align=' + ("TRUE" if align else "FALSE") + ', npoly=' + str(nPolynom) + ', maxGroupDist=' + str(int(maxTimeDiff)) + ', scanDuration=' + str(scanDuration) + ')')

        # convert R-results object to python arrays
        retl = []
        if len(ret) > 1:
            for i in range(0, len(ret), 2):
                cur = []
                for j in range(0, 2):
                    cur.append(ret[i + j])
                retl.append(cur)

        # match aligned results to the supplied chromatographic peaks
        retld = []
        i = 0
        for peaks in peakss:
            j = 0
            dd = []
            for peak in peaks:
                dd.append((retl[(i * maxPeaks) + j]))
                j = j + 1
            retld.append(dd)
            i = i + 1

        return retld

    # returns the algined EICs
    def getAligendXICs(self, eics, scantimes, align=True, nPolynom=3, pretend=100, scanDuration=1):
        assert (len(eics) == len(scantimes))
        if len(eics) == 0:
            return [[]]

        eics, scantimes= preAppendEICs(eics, scantimes, pretend=pretend, scanDuration=scanDuration)

        nPolynom = max(nPolynom, 0)
        maxTime = max([max(scanTime) for scanTime in scantimes])
        refTimes = range(int(maxTime))

        if len(eics) == 1:
            retl, refTimes=undoPreAppendEICs([mapArrayToRefTimes(eics[0], scantimes[0], refTimes)], [u for u in refTimes], pretend=pretend, scanDuration=scanDuration)
            undoPreAppendEICs(retl, [u for u in refTimes], pretend=pretend)

        maxTime = max([max(scanTime) for scanTime in scantimes])
        refTimes = range(int(maxTime))

        eicds = []
        for i in range(len(eics)):
            eicds.append(mapArrayToRefTimes(eics[i], scantimes[i], refTimes))
        eics = eicds

        maxLen = max([len(eic) for eic in eics])
        eicsR = ""
        append = False
        for eic in eics:
            if append:
                eicsR = eicsR + ","
            else:
                append = True
            eicsR = eicsR + ",".join([str(d) for d in eic])
            if len(eic) < maxLen:
                eicsR = eicsR + "," + ",".join([str(0) for p in range(len(eic), maxLen)])

        ret = r('getAlignXICs(c(' + eicsR + '),  ' + str(len(eics)) + ', align=' + ("TRUE" if align else "FALSE") + ', npoly=' + str(nPolynom) + ')')
        retl = []
        for i in range(len(eics)):
            d = []
            for j in range(maxLen):
                d.append(ret[(i * maxLen) + j])
            retl.append(d)

        retl, refTimes=undoPreAppendEICs(retl, [u for u in refTimes], pretend=pretend, scanDuration=scanDuration)
        return retl, refTimes









