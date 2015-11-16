#    MetExtract II
#    Copyright (C) 2015
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


from utils import getNormRatio

maxSub = 5


def getSubstitutionArray(purity, xMax, maxSub):
    ret = []
    ret.append([-1 for n in range(0, maxSub + 1)])
    for i in range(1, xMax + 1):
        cur = []
        for j in range(0, maxSub + 1):
            cur.append(getNormRatio(purity, i, j))
        ret.append(cur)
    return ret


# results object (ion signal pairs)
class mzFeature:
    def __init__(self, mz, lmz, xCount, scanIndex, loading, nIntensity, lIntensity, ionMode):
        self.mz = mz
        self.lmz = lmz
        self.xCount = xCount
        self.scanIndex = scanIndex
        self.loading = loading
        self.nIntensity = nIntensity
        self.lIntensity = lIntensity
        self.ionMode = ionMode

    def __str__(self):
        return "mz: %.4f, lmz: %.4f, xCount: %d, charge: %d, NInt: %.1f, LInt: %.1f, ionMode: %s, scanIndex: %d"%(
                    self.mz, self.lmz, self.xCount, self.loading, self.nIntensity, self.lIntensity, self.ionMode, self.scanIndex)

# detects in each recorded MS scan (lvl. 1) isotope patterns originating from a native and a (partially) labelled
# metabolite / biotransformation product. It is important, that all atoms that can be labelled are actually
# labelled and have the same chance to be labelled. Thus, fluxomics applications or such isotope patterns
# are in general not supported. They can, however, be detected if not isotopolog verification step is
# used (not recommended)
def matchPartners(mzXMLData, labellingElement, useCValidation, intensityThres, isotopologIntensityThres, maxLoading, xMin, xMax, xOffset, ppm,
                  intensityErrorN, intensityErrorL, purityN, purityL, startTime, stopTime, filterLine, ionMode,
                  peakCountLeft, peakCountRight, lowAbundanceIsotopeCutoff, metabolisationExperiment, checkRatio,
                  minRatio, maxRatio, reportFunction=None):
    scanRTRange = stopTime - startTime

    cValidationOffset = 1.00335   # mass difference between 12C and 13C

    detectedIonPairs = []

    oriXOffset = xOffset
    oriCValidationOffset = cValidationOffset

    if labellingElement == 'C':
        useCValidation = 0
        # Carbon validation is the only possible method
        # when the labelling element is carbon (to some extend that's logical)

    # substitution arrays for checking the number of carbon atoms
    purityNArray = getSubstitutionArray(purityN, xMax + 3, maxSub)   # native metabolite
    purityLArray = getSubstitutionArray(purityL, xMax + 3, maxSub)   # labelled metabolite

    # iterate over all MS scans (lvl. 1)
    curScanIndex=0
    for j in range(0, len(mzXMLData.MS1_list)):

        try:
            curScan = mzXMLData.MS1_list[j]
            curScanDetectedIonPairs = []

            # check for correct filterline and scan time
            if curScan.filter_line == filterLine:
                if startTime <= (curScan.retention_time / 60.) <= stopTime:

                    if reportFunction is not None:
                        reportFunction((curScan.retention_time / 60. - startTime) / scanRTRange,
                                       "RT %.2f" % (curScan.retention_time / 60.))

                    dontUsePeakIndices = []

                    # assume each peak to be a valid M (monoisotopic, native metabolite ion)
                    # and verify this assumption (search for (partially) labelled pendant)
                    for currentPeakIndex in range(0, len(curScan.mz_list)):
                        if not (currentPeakIndex in dontUsePeakIndices):
                            curPeakmz = curScan.mz_list[currentPeakIndex]
                            curPeakIntensity = curScan.intensity_list[currentPeakIndex]

                            curPeakDetectedIonPairs = []

                            # only consider peaks above the threshold
                            if curPeakIntensity >= intensityThres:
                                skipOtherLoadings = False


                                possibleLoadings=[]
                                ## figure out possible loadings
                                for l in range(1, maxLoading+1, 1):
                                    iso = curScan.findMZ(curPeakmz + oriXOffset / l, ppm, start=currentPeakIndex)
                                    iso = curScan.getMostIntensePeak(iso[0], iso[1])

                                    if iso != -1:
                                        possibleLoadings.append(l)

                                if len(possibleLoadings)==0:
                                    possibleLoadings=[1]


                                for curLoading in possibleLoadings:
                                    if not skipOtherLoadings:
                                        xOffset = oriXOffset / float(curLoading)
                                        cValidationOffset = oriCValidationOffset / float(curLoading)


                                        # C-isotope distribution validation for labelling with N, S, ... (useCValidation == 2)
                                        # region
                                        # The carbon distribution of both isotopologs is checked for equality
                                        # checks if the isotope patterns of M, M+1.. and M', M'+1.. are approximately the same
                                        # E.g. 15N-labelling
                                        # |    <--   Nn   -->    |
                                        # ||                     ||
                                        # |||                   ||||
                                        # Required for some labelling applications (e.g. S, N, Cl)
                                        # Requires: - a high resolution and separation of different isotoplogs (especially carbon)
                                        # EXPERIMENTAL: has not been tested with real data (not N or S labelled sample material
                                        #               was available)
                                        if useCValidation == 2:

                                            # search for M+1 isotoplog (if required it needs to be present; slight speed gain)
                                            isoM_p1 = curScan.findMZ(curPeakmz + cValidationOffset, ppm, start=currentPeakIndex)
                                            isoM_p1 = curScan.getMostIntensePeak(isoM_p1[0], isoM_p1[1])
                                            if isoM_p1 != -1:
                                                dontUseXCount = []

                                                # test certain number of labelled atoms
                                                for xCount in range(xMax, xMin - 1, -1):
                                                    if not (xCount in dontUseXCount):

                                                        # search for M'
                                                        isoM_pX = curScan.findMZ(curPeakmz + xCount * xOffset, ppm,
                                                                                 start=currentPeakIndex)
                                                        isoM_pX = curScan.getMostIntensePeak(isoM_pX[0], isoM_pX[1],
                                                                                             intensityThres)
                                                        if isoM_pX != -1:
                                                            # check if ratio is okay
                                                            rat=curPeakIntensity/curScan.intensity_list[isoM_pX]
                                                            if not(checkRatio) or minRatio<=rat<=maxRatio:

                                                                # search for M'+1
                                                                isoM_pXp1 = curScan.findMZ(
                                                                    curPeakmz + xCount * xOffset + cValidationOffset, ppm,
                                                                    start=currentPeakIndex)
                                                                isoM_pXp1 = curScan.getMostIntensePeak(isoM_pXp1[0],
                                                                                                       isoM_pXp1[1])
                                                                if isoM_pXp1 != -1:
                                                                    acceptIsoPeak = []

                                                                    if not (isoM_p1 in dontUsePeakIndices):
                                                                        acceptIsoPeak.append(
                                                                            curScan.intensity_list[isoM_p1] / curPeakIntensity)
                                                                    acceptIsoPeakO = []

                                                                    # check if ratios of M+1/M and M'+1/M' are approximately
                                                                    # the same (<=intensityErrorL)
                                                                    if len(acceptIsoPeak) > 0:

                                                                        if not (isoM_pX in dontUsePeakIndices):
                                                                            if not (isoM_pXp1 in dontUsePeakIndices):
                                                                                for isoRatio in acceptIsoPeak:
                                                                                    if abs(curScan.intensity_list[isoM_pXp1] /
                                                                                                   curScan.intensity_list[
                                                                                                       isoM_pX] - isoRatio) < intensityErrorL:
                                                                                        acceptIsoPeakO.append(
                                                                                            (isoM_pX, isoM_pXp1))
                                                                    if len(acceptIsoPeak) > 0 and len(acceptIsoPeakO) > 0:
                                                                        #mzs.append(
                                                                        #    [mz, xCount, curScanNum, curLoading, intensity])
                                                                        curScanDetectedIonPairs.append(
                                                                            mzFeature(mz=curPeakmz,
                                                                                      lmz=curScan.mz_list[isoM_pX],
                                                                                      xCount=xCount,
                                                                                      scanIndex=curScanIndex,
                                                                                      loading=curLoading,
                                                                                      nIntensity=curPeakIntensity,
                                                                                      lIntensity=curScan.intensity_list[acceptIsoPeakO[0][0]],
                                                                                      ionMode=ionMode))
                                                                        skipOtherLoadings = True
                                        # endregion
                                        # Mixed isotope validation (useCValidation == 1)
                                        # region
                                        # First, the distinct isotope patterns of the native and
                                        # x-labelled metabolite forms are tested. If these do not successfully validate
                                        # a ion signal pairing, C-validation (similar carbon isotope pattern) is tested
                                        # Required for some labelling applications (e.g. S, N, Cl)
                                        # Requires: - a high resolution and separation of different isotoplogs (especially carbon)
                                        # EXPERIMENTAL: has not been tested with real data (not N or S labelled sample material
                                        #               was available)
                                        # CAUTION: Not tested, partially implemented. Use with caution
                                        if useCValidation == 1:
                                            #Mixed Isotope Validation
                                            isoM_p1 = curScan.findMZ(curPeakmz + xOffset, ppm, start=currentPeakIndex)
                                            isoM_p1 = curScan.getMostIntensePeak(isoM_p1[0], isoM_p1[1])
                                            if peakCountLeft == 1 or isoM_p1 != -1:
                                                dontUseXCount = []
                                                for xCount in range(xMax, xMin - 1, -1):
                                                    if not (xCount in dontUseXCount):
                                                        isoM_pX = curScan.findMZ(curPeakmz + xCount * xOffset, ppm,
                                                                                 start=currentPeakIndex)
                                                        isoM_pX = curScan.getMostIntensePeak(isoM_pX[0], isoM_pX[1],
                                                                                             intensityThres)
                                                        if isoM_pX != -1:
                                                            isoM_pXm1 = curScan.findMZ(curPeakmz + (xCount - 1) * xOffset, ppm,
                                                                                       start=currentPeakIndex)
                                                            isoM_pXm1 = curScan.getMostIntensePeak(isoM_pXm1[0],
                                                                                                   isoM_pXm1[1])
                                                            if peakCountRight == 1 or isoM_pXm1[0] != -1:
                                                                acceptIsoPeak = []
                                                                normRatio = purityNArray[xCount][1]
                                                                if peakCountLeft > 1:
                                                                    if metabolisationExperiment:
                                                                        isoM_pXp1 = curScan.findMZ(
                                                                            curPeakmz + (xCount + 1) * xOffset, ppm,
                                                                            start=currentPeakIndex)
                                                                        isoM_pXp1 = curScan.getMostIntensePeak(
                                                                            isoM_pXp1[0], isoM_pXp1[1])


                                                                        if not (isoM_p1 in dontUsePeakIndices):
                                                                            if abs(curScan.intensity_list[
                                                                                       isoM_p1] / curPeakIntensity - normRatio -
                                                                                           curScan.intensity_list[
                                                                                               isoM_pXp1]) < intensityErrorN:
                                                                                acceptIsoPeak.append(isoM_p1)
                                                                    else:

                                                                        if not (isoM_p1 in dontUsePeakIndices):
                                                                            if abs(curScan.intensity_list[
                                                                                       isoM_p1] / curPeakIntensity - normRatio) < intensityErrorN:
                                                                                acceptIsoPeak.append(isoM_p1)
                                                                else:
                                                                    acceptIsoPeak.append(1)
                                                                acceptIsoPeakO = []
                                                                normRatio = purityLArray[xCount][1]

                                                                if peakCountRight > 1:

                                                                    if not (isoM_pX in dontUsePeakIndices):
                                                                        if not (isoM_pXm1 in dontUsePeakIndices):
                                                                            if abs(curScan.intensity_list[isoM_pXm1] /
                                                                                           curScan.intensity_list[
                                                                                               isoM_pX] - normRatio) < intensityErrorL:
                                                                                acceptIsoPeakO.append((isoM_pX, isoM_pXm1))
                                                                else:
                                                                    acceptIsoPeakO.append(1)
                                                                if (peakCountLeft == 1 or len(acceptIsoPeak) > 0) and (
                                                                                peakCountRight == 1 or len(
                                                                            acceptIsoPeakO) > 0):
                                                                    #mzs.append(
                                                                    #    [mz, xCount, curScanNum, curLoading, intensity])

                                                                    curScanDetectedIonPairs.append(
                                                                        mzFeature(mz=curPeakmz,
                                                                                  xCount=xCount,
                                                                                  lmz=curScan.mz_list[isoM_pX],
                                                                                  scanIndex=curScanIndex,
                                                                                  loading=curLoading,
                                                                                  nIntensity=curPeakIntensity,
                                                                                  lIntensity=curScan.intensity_list[acceptIsoPeakO[0][0]],
                                                                                  ionMode=ionMode))
                                                                    skipOtherLoadings = True
                                                                    #dontUseXCount.append(xCount-1)

                                                            else:
                                                                isoM_pXm1 = curScan.findMZ(
                                                                    curPeakmz + (xCount + 1) * cValidationOffset, ppm,
                                                                    start=currentPeakIndex)
                                                                isoM_pXm1 = curScan.getMostIntensePeak(isoM_pXm1[0],
                                                                                                       isoM_pXm1[1])
                                                                if isoM_pXm1 != -1:
                                                                    acceptIsoPeak = []

                                                                    if not (isoM_p1 in dontUsePeakIndices):
                                                                        acceptIsoPeak.append(
                                                                            curScan.intensity_list[isoM_p1] / curPeakIntensity)
                                                                    acceptIsoPeakO = []

                                                                    if len(acceptIsoPeak) > 0:
                                                                        if not (isoM_pX in dontUsePeakIndices):
                                                                            if not (isoM_pXm1 in dontUsePeakIndices):
                                                                                for isoRatio in acceptIsoPeak:
                                                                                    if abs(curScan.intensity_list[
                                                                                               isoM_pXm1] /
                                                                                                   curScan.intensity_list[
                                                                                                       isoM_pX] - isoRatio) < intensityErrorL:
                                                                                        acceptIsoPeakO.append(
                                                                                            (isoM_pX, isoM_pXm1))
                                                                    if len(acceptIsoPeak) > 0 and len(
                                                                            acceptIsoPeakO) > 0:
                                                                        #mzs.append([mz, xCount, curScanNum, curLoading,
                                                                        #            intensity])

                                                                        curScanDetectedIonPairs.append(
                                                                            mzFeature(mz=curPeakmz,
                                                                                      lmz=curScan.mz_list[isoM_pX],
                                                                                      xCount=xCount,
                                                                                      scanIndex=curScanIndex,
                                                                                      loading=curLoading,
                                                                                      nIntensity=curPeakIntensity,
                                                                                      lIntensity=curScan.intensity_list[acceptIsoPeakO[0][0]],
                                                                                      ionMode=ionMode))
                                                                        skipOtherLoadings = True
                                        # endregion
                                        # Isotope pattern validation (useCValidation == 0)
                                        # region
                                        # It is tested, if the expected isotope patterns
                                        # separate for the native and labelled metabolite follow a theoretical pattern
                                        # E.g. native 12C and uniformly / partially 13C-labelled metabolites
                                        # |    <--   Cn   -->    |
                                        # ||                    ||
                                        # |||                  ||||
                                        # Necessary mainly for 13C-labelling with mirror-symmetric isotope patterns
                                        # NOTE: - Approach is mainly used for 13C-labelling
                                        if useCValidation==0:
                                            # find M+1 peak
                                            isoM_p1 = curScan.findMZ(curPeakmz + xOffset, ppm, start=currentPeakIndex)
                                            isoM_p1 = curScan.getMostIntensePeak(isoM_p1[0], isoM_p1[1])
                                            if isoM_p1 != -1 or peakCountLeft == 1 or lowAbundanceIsotopeCutoff:
                                                # test certain number of labelled carbon atoms

                                                for xCount in range(xMax, xMin - 1, -1):
                                                    # find corresponding M' peak
                                                    isoM_pX = curScan.findMZ(curPeakmz + xCount * xOffset, ppm, start=currentPeakIndex)
                                                    isoM_pX = curScan.getMostIntensePeak(isoM_pX[0], isoM_pX[1], intensityThres)
                                                    if isoM_pX != -1:

                                                        labPeakmz=curScan.mz_list[isoM_pX]
                                                        labPeakIntensity=curScan.intensity_list[isoM_pX]

                                                        # (1.) check if M' and M ratio are as expected
                                                        if checkRatio:
                                                            rat = curPeakIntensity / labPeakIntensity
                                                            if minRatio <= rat <= maxRatio:
                                                                pass     ## ratio check passed
                                                            else:
                                                                continue ## ratio check not passed

                                                        ## no isotopolog verification needs to be performed
                                                        if peakCountLeft == 1 and peakCountRight == 1:
                                                            curPeakDetectedIonPairs.append(
                                                                mzFeature(mz=curPeakmz,
                                                                          lmz=curScan.mz_list[isoM_pX],
                                                                          xCount=xCount,
                                                                          scanIndex=curScanIndex,
                                                                          loading=curLoading,
                                                                          nIntensity=curPeakIntensity,
                                                                          lIntensity=labPeakIntensity,
                                                                          ionMode=ionMode))

                                                            skipOtherLoadings = True
                                                            continue

                                                        # find M'-1 peak
                                                        isoM_pXm1 = curScan.findMZ(curPeakmz + (xCount - 1) * xOffset, ppm, start=currentPeakIndex)
                                                        isoM_pXm1 = curScan.getMostIntensePeak(isoM_pXm1[0], isoM_pXm1[1])
                                                        normRatioL = purityLArray[xCount][1]
                                                        normRatioN = purityNArray[xCount][1]

                                                        # 2. check if the observed M'-1/M' ratio fits the theoretical one
                                                        if peakCountRight > 1:
                                                            if isoM_pXm1 != -1:
                                                                observedRatioMp=curScan.intensity_list[isoM_pXm1] / labPeakIntensity
                                                                if abs(normRatioL-observedRatioMp) <= intensityErrorL:
                                                                    pass     ## accept peak
                                                                else:
                                                                    continue ## discard current peak
                                                            elif lowAbundanceIsotopeCutoff:
                                                                if labPeakIntensity*normRatioL <= isotopologIntensityThres:
                                                                    pass     ## accept peak
                                                                else:
                                                                    continue ## discard current peak
                                                            else:
                                                                continue     ## discard current peak
                                                        else:
                                                            pass ## accept peak

                                                        # 3. check if the observed M+1/M ratoi fiths the theoretical one
                                                        if peakCountLeft > 1:
                                                            observedRatioM = curScan.intensity_list[isoM_p1] / curPeakIntensity
                                                            adjRatio=0

                                                            if metabolisationExperiment:
                                                                # if experiment is a tracer-fate study, assume a conjugated moiety and correct M+1/M ratio for it
                                                                isoM_pXp1 = curScan.findMZ( curPeakmz + (xCount + 1) * xOffset, ppm, start=currentPeakIndex)
                                                                isoM_pXp1 = curScan.getMostIntensePeak( isoM_pXp1[0], isoM_pXp1[1])
                                                                if isoM_pXp1 != -1:
                                                                    adjRatio=curScan.intensity_list[isoM_pXp1] / labPeakIntensity

                                                            if abs(abs(observedRatioM-adjRatio)-normRatioN) <= intensityErrorN:
                                                                pass         ## acceptPeak
                                                            elif lowAbundanceIsotopeCutoff:
                                                                if curPeakIntensity*normRatioN <= isotopologIntensityThres:
                                                                    pass     ## accept peak
                                                                else:
                                                                    continue ## discard peak
                                                            else:
                                                                continue     ## discard current peak




                                                        # All verification criteria are passed, store the ion signal pair
                                                        # for further processing
                                                        curPeakDetectedIonPairs.append(
                                                            mzFeature(mz=curPeakmz,
                                                                      lmz=curScan.mz_list[isoM_pX],
                                                                      xCount=xCount,
                                                                      scanIndex=curScanIndex,
                                                                      loading=curLoading,
                                                                      nIntensity=curPeakIntensity,
                                                                      lIntensity=labPeakIntensity,
                                                                      ionMode=ionMode))

                                                        skipOtherLoadings = True
                                                        # endregion

                            if True:  ## select best fit
                                if len(curPeakDetectedIonPairs)>0:
                                    bestFit=None
                                    bestFitPPMDiff=1000000000

                                    for mt in curPeakDetectedIonPairs:
                                        if abs(mt.lmz-mt.mz-mt.xCount*1.00335)*1000000./mt.mz < bestFitPPMDiff:
                                            bestFit=mt
                                            bestFitPPMDiff=abs(mt.lmz-mt.mz-mt.xCount*1.00335)*1000000./mt.mz

                                    curScanDetectedIonPairs.append(bestFit)
                            else:     ## use all peak pairs
                                curScanDetectedIonPairs.extend(curPeakDetectedIonPairs)

                detectedIonPairs.extend(curScanDetectedIonPairs)
                curScanIndex += 1

        except Exception as e:
            import traceback
            traceback.print_exc()

    return detectedIonPairs


















