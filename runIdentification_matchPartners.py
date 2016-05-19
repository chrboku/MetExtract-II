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
    def __init__(self, mz, lmz, tmz, xCount, scanIndex, loading, nIntensity, lIntensity, ionMode):
        self.mz = mz
        self.lmz = lmz
        self.tmz = tmz
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
                                for l in range(maxLoading, 0, -1):
                                    iso = curScan.findMZ(curPeakmz + oriCValidationOffset / l, ppm, start=currentPeakIndex)
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
                                        # region
                                        if useCValidation == 2:
                                            # find M+1 peak
                                            isoM_p1 = curScan.findMZ(curPeakmz + cValidationOffset, ppm, start=currentPeakIndex)
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
                                                                          lmz=labPeakmz,
                                                                          tmz=xCount * xOffset,
                                                                          xCount=xCount,
                                                                          scanIndex=curScanIndex,
                                                                          loading=curLoading,
                                                                          nIntensity=curPeakIntensity,
                                                                          lIntensity=labPeakIntensity,
                                                                          ionMode=ionMode))

                                                            skipOtherLoadings = True
                                                            continue

                                                        # find M'+1 peak
                                                        isoM_pXp1 = curScan.findMZ(curPeakmz + xCount * xOffset + cValidationOffset, ppm, start=currentPeakIndex)
                                                        isoM_pXp1 = curScan.getMostIntensePeak(isoM_pXp1[0], isoM_pXp1[1])

                                                        # calculate the ratio of M+1/M
                                                        isoPeakIntensity=curScan.intensity_list[isoM_p1]
                                                        if peakCountLeft == 1:
                                                            ratioN=None
                                                        elif peakCountLeft > 1 and lowAbundanceIsotopeCutoff and isoPeakIntensity <= isotopologIntensityThres:
                                                            ratioN=None
                                                        else:
                                                            ratioN=isoPeakIntensity/curPeakIntensity

                                                        # calculate the ratio of M'+1/M'
                                                        isoLabPeakIntensity=curScan.intensity_list[isoM_pXp1]
                                                        if peakCountRight == 1:
                                                            ratioL=None
                                                        elif peakCountRight > 1 and lowAbundanceIsotopeCutoff and isoLabPeakIntensity <= isotopologIntensityThres:
                                                            ratioL=None
                                                        else:
                                                            ratioL=isoLabPeakIntensity/labPeakIntensity

                                                        # 2. check if the observed M'+1/M' ratio and the M+1/M ratio are approximately equal
                                                        if (ratioN!=None and ratioL!=None) and abs(ratioN-ratioL)<=intensityErrorL:
                                                            curPeakDetectedIonPairs.append(
                                                                mzFeature(mz=curPeakmz,
                                                                          lmz=labPeakmz,
                                                                          tmz=xCount * xOffset,
                                                                          xCount=xCount,
                                                                          scanIndex=curScanIndex,
                                                                          loading=curLoading,
                                                                          nIntensity=curPeakIntensity,
                                                                          lIntensity=labPeakIntensity,
                                                                          ionMode=ionMode))

                                                            skipOtherLoadings = True
                                                            continue
                                                        elif lowAbundanceIsotopeCutoff and (ratioN==None or ratioL==None):
                                                            curPeakDetectedIonPairs.append(
                                                                mzFeature(mz=curPeakmz,
                                                                          lmz=labPeakmz,
                                                                          tmz=xCount * xOffset,
                                                                          xCount=xCount,
                                                                          scanIndex=curScanIndex,
                                                                          loading=curLoading,
                                                                          nIntensity=curPeakIntensity,
                                                                          lIntensity=labPeakIntensity,
                                                                          ionMode=ionMode))

                                                            skipOtherLoadings = True
                                                            continue
                                        # endregion




                                        # Isotope pattern validation (useCValidation == 0)
                                        # It is tested, if the expected isotope patterns
                                        # separate for the native and labelled metabolite follow a theoretical pattern
                                        # E.g. native 12C and uniformly / partially 13C-labelled metabolites
                                        # |    <--   Cn   -->    |
                                        # ||                    ||
                                        # |||                  ||||
                                        # Necessary mainly for 13C-labelling with mirror-symmetric isotope patterns
                                        # NOTE: - Approach is mainly used for 13C-labelling
                                        # region
                                        if useCValidation==0:
                                            # find M+1 peak
                                            isoM_p1 = curScan.findMZ(curPeakmz + cValidationOffset, ppm, start=currentPeakIndex)
                                            isoM_p1 = curScan.getMostIntensePeak(isoM_p1[0], isoM_p1[1])
                                            if isoM_p1 != -1 or peakCountLeft == 1 or lowAbundanceIsotopeCutoff:
                                                # test certain number of labelled carbon atoms

                                                for xCount in range(xMax, xMin - 1, -1):
                                                    # find corresponding M' peak
                                                    isoM_pX = curScan.findMZ(curPeakmz + xCount * cValidationOffset, ppm, start=currentPeakIndex)
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
                                                                          lmz=labPeakmz,
                                                                          tmz=xCount * xOffset,
                                                                          xCount=xCount,
                                                                          scanIndex=curScanIndex,
                                                                          loading=curLoading,
                                                                          nIntensity=curPeakIntensity,
                                                                          lIntensity=labPeakIntensity,
                                                                          ionMode=ionMode))

                                                            skipOtherLoadings = True
                                                            continue

                                                        # find M'-1 peak
                                                        isoM_pXm1 = curScan.findMZ(curPeakmz + (xCount - 1) * cValidationOffset, ppm, start=currentPeakIndex)
                                                        isoM_pXm1 = curScan.getMostIntensePeak(isoM_pXm1[0], isoM_pXm1[1])
                                                        normRatioL = purityLArray[xCount][1]
                                                        normRatioN = purityNArray[xCount][1]


                                                        # 2. check if the observed M'-1/M' ratio fits the theoretical one
                                                        if peakCountRight > 1 or (lowAbundanceIsotopeCutoff and labPeakIntensity*normRatioL <= isotopologIntensityThres):
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

                                                        # 3. check if the observed M+1/M ratio fits the theoretical one
                                                        if peakCountLeft > 1 or (lowAbundanceIsotopeCutoff and curPeakIntensity*normRatioN <= isotopologIntensityThres):
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
                                                                      lmz=labPeakmz,
                                                                      tmz=xCount * xOffset,
                                                                      xCount=xCount,
                                                                      scanIndex=curScanIndex,
                                                                      loading=curLoading,
                                                                      nIntensity=curPeakIntensity,
                                                                      lIntensity=labPeakIntensity,
                                                                      ionMode=ionMode))

                                                        skipOtherLoadings = True
                                        # endregion

                            if False:  ## select best fit
                                if len(curPeakDetectedIonPairs)>0:
                                    bestFit=None
                                    bestFitPPMDiff=1000000000

                                    ## TODO select best fit based on isotopic pattern (e.g. intensity)


                                    for mt in curPeakDetectedIonPairs:
                                        if abs(mt.lmz-mt.mz-mt.xCount*1.00335)*1000000./mt.mz < bestFitPPMDiff:
                                            bestFit=mt
                                            bestFitPPMDiff=abs(mt.lmz-mt.mz-mt.xCount*1.00335)*1000000./mt.mz

                                    curScanDetectedIonPairs.append(bestFit)
                            else:     ## use all peak pairs
                                curScanDetectedIonPairs.extend(curPeakDetectedIonPairs)

                if len(curScanDetectedIonPairs)>0 and False:
                    from utils import printObjectsAsTable
                    print "\n"
                    print curScan.retention_time/60.
                    printObjectsAsTable(curScanDetectedIonPairs, attrs=["mz", "xCount", "loading", "nIntensity", "lIntensity", "ionMode"])

                detectedIonPairs.extend(curScanDetectedIonPairs)
                curScanIndex += 1

        except Exception as e:
            import traceback
            traceback.print_exc()

    return detectedIonPairs


















