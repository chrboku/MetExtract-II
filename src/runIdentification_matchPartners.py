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


from .utils import getNormRatio
from .utils import Bunch
from copy import deepcopy
from .formulaTools import formulaTools

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
        return "mz: %.4f, lmz: %.4f, xCount: %d, charge: %d, NInt: %.1f, LInt: %.1f, ionMode: %s, scanIndex: %d" % (
            self.mz,
            self.lmz,
            self.xCount,
            self.loading,
            self.nIntensity,
            self.lIntensity,
            self.ionMode,
            self.scanIndex,
        )


## calculates all possible combinations of labeling elements
## required for double labeling
def getCombinationsOfLabel(
    useElems,
    labelingElements,
    minLabelingAtoms,
    maxLabelingAtoms,
    used=None,
    startAt=0,
    ind=2,
):
    combs = []
    if used is None:
        used = {}

    x = sum(used.values())
    if minLabelingAtoms <= x <= maxLabelingAtoms:
        b = Bunch(
            atoms=deepcopy(used),
            atomsCount=sum([used[e] for e in used.keys()]),
            mzdelta=sum([(labelingElements[e].massLabeled - labelingElements[e].massNative) * used[e] for e in used.keys()]),
        )
        combs.append(b)

    if startAt < len(useElems):
        for cStartAt in range(startAt, len(useElems)):
            e = useElems[cStartAt]
            for i in range(labelingElements[e].minXn, labelingElements[e].maxXn + 1):
                c = deepcopy(used)
                c[e] = i

                combs.extend(
                    getCombinationsOfLabel(
                        useElems,
                        labelingElements,
                        minLabelingAtoms,
                        maxLabelingAtoms,
                        c,
                        cStartAt + 1,
                        ind=ind + 2,
                    )
                )

    return combs


# detects in each recorded MS scan (lvl. 1) isotope patterns originating from a native and a (partially) labelled
# metabolite / biotransformation product. It is important, that all atoms that can be labelled are actually
# labelled and have the same chance to be labelled. Thus, fluxomics applications or such isotope patterns
# are in general not supported. They can, however, be detected if not isotopolog verification step is
# used (not recommended)
def matchPartners(
    mzXMLData,
    forFile,
    labellingElement,
    useCIsotopePatternValidation,
    intensityThres,
    isotopologIntensityThres,
    maxLoading,
    xCounts,
    xOffset,
    ppm,
    intensityErrorN,
    intensityErrorL,
    purityN,
    purityL,
    startTime,
    stopTime,
    filterLine,
    ionMode,
    peakCountLeft,
    peakCountRight,
    lowAbundanceIsotopeCutoff,
    metabolisationExperiment,
    checkRatio,
    minRatio,
    maxRatio,
    reportFunction=None,
    writeExtendedDiagnostics=True,
):
    # debug: print all input parameters
    try:
        params = {k: v for k, v in locals().items()}
        print("matchPartners parameters:\n" + str(params))
    except Exception as _e:
        try:
            for k, v in locals().items():
                try:
                    print("%s: %r" % (k, v))
                except Exception:
                    print("%s: <unprintable>" % k)
        except Exception:
            print("Failed to print parameters")

    scanRTRange = stopTime - startTime

    cValidationOffset = 1.00335484  # mass difference between 12C and 13C

    detectedIonPairs = []

    oriXOffset = xOffset
    oriCValidationOffset = cValidationOffset

    if labellingElement == "C":
        useCIsotopePatternValidation = 0
        # Carbon validation is the only possible method
        # when the labelling element is carbon (to some extend that's logical)

    # substitution arrays for checking the number of carbon atoms
    purityNArray = getSubstitutionArray(purityN, max(xCounts) + 1, maxSub)  # native metabolite
    purityLArray = getSubstitutionArray(purityL, max(xCounts) + 1, maxSub)  # labelled metabolite

    labelingElements = {}
    labelingElements["C"] = Bunch(
        nativeIsotope="12C",
        labelingIsotope="13C",
        massNative=12.0,
        massLabeled=13.00335,
        isotopicEnrichmentNative=0.9893,
        isotopicEnrichmentLabeled=0.995,
        minXn=1,
        maxXn=3,
    )
    labelingElements["H"] = Bunch(
        nativeIsotope="1H",
        labelingIsotope="2H",
        massNative=1.0078250,
        massLabeled=2.0141018,
        isotopicEnrichmentNative=0.9999,
        isotopicEnrichmentLabeled=0.96,
        minXn=0,
        maxXn=9,
    )

    minLabelingAtoms = 5
    maxLabelingAtoms = 5
    ## combinations of labeling elements
    tempCombs = getCombinationsOfLabel(["C", "H"], labelingElements, 2, 12)
    combs = []
    for comb in tempCombs:
        ## remove implausible labeling combinations
        if True:
            if "C" in comb.atoms.keys() and comb.atoms["C"] == 1 and ("H" not in comb.atoms.keys() or ("H" in comb.atoms.keys() and comb.atoms["H"] in [1, 2, 3])):
                combs.append(comb)
            if "C" in comb.atoms.keys() and comb.atoms["C"] == 2 and ("H" not in comb.atoms.keys() or ("H" in comb.atoms.keys() and comb.atoms["H"] in [1, 2, 3, 4, 5, 6])):
                combs.append(comb)
            if "C" in comb.atoms.keys() and comb.atoms["C"] == 3 and ("H" not in comb.atoms.keys() or ("H" in comb.atoms.keys() and comb.atoms["H"] in [1, 2, 3, 4, 5, 6, 7, 8, 9])):
                combs.append(comb)

    fT = formulaTools()

    useDoubleLabelingCombinations = False

    if useDoubleLabelingCombinations:
        print("The following labeling configurations are used:")
        for comb in combs:
            print(comb)

    # iterate over all MS scans (lvl. 1)
    curScanIndex = 0
    for j in range(0, len(mzXMLData.MS1_list)):
        try:
            curScan = mzXMLData.MS1_list[j]
            curScanDetectedIonPairs = []

            # check for correct filterline and scan time
            if curScan.filter_line == filterLine:
                if startTime <= (curScan.retention_time / 60.0) <= stopTime:
                    if reportFunction is not None:
                        reportFunction(
                            (curScan.retention_time / 60.0 - startTime) / scanRTRange,
                            "RT %.2f" % (curScan.retention_time / 60.0),
                        )

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

                                ## do not process peaks that are likely isotopologs
                                backIsos = []
                                for l in range(maxLoading, 0, -1):
                                    iso = curScan.findMZ(curPeakmz - oriCValidationOffset / l, ppm)
                                    iso = curScan.getMostIntensePeak(iso[0], iso[1])

                                    if iso != -1 and curScan.intensity_list[iso] > curPeakIntensity:
                                        backIsos.append(l)
                                if len(backIsos) > 0:
                                    continue

                                possibleLoadings = []
                                ## figure out possible loadings
                                for l in range(maxLoading, 0, -1):
                                    iso = curScan.findMZ(
                                        curPeakmz + oriCValidationOffset / l,
                                        ppm,
                                        start=currentPeakIndex,
                                    )
                                    iso = curScan.getMostIntensePeak(iso[0], iso[1])

                                    if iso != -1:
                                        possibleLoadings.append(l)
                                        break  ## skip other loadings

                                if len(possibleLoadings) == 0:
                                    possibleLoadings = [1]

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
                                        # |||                   .|||
                                        # Required for some labelling applications (e.g. S, N, Cl)
                                        # Requires: - a high resolution and separation of different isotoplogs (especially carbon)
                                        # EXPERIMENTAL: has not been tested with real data (not N or S labelled sample material
                                        #               was available)
                                        if not useDoubleLabelingCombinations and useCIsotopePatternValidation != 0:
                                            isoM_m1 = curScan.findMZ(curPeakmz - cValidationOffset, ppm * 2)
                                            isoM_m1 = curScan.getMostIntensePeak(isoM_m1[0], isoM_m1[1])
                                            if isoM_m1 != -1:
                                                obRatio = curScan.intensity_list[isoM_m1] / curPeakIntensity
                                                if obRatio >= 0.5:
                                                    continue

                                            # find M+1 peak
                                            isoM_p1 = curScan.findMZ(
                                                curPeakmz + cValidationOffset,
                                                ppm,
                                                start=currentPeakIndex,
                                            )
                                            isoM_p1 = curScan.getMostIntensePeak(isoM_p1[0], isoM_p1[1])
                                            if isoM_p1 != -1 or peakCountLeft == 1 or lowAbundanceIsotopeCutoff:
                                                # test certain number of labelled carbon atoms

                                                for xCount in sorted(xCounts, reverse=True):
                                                    # find corresponding M' peak
                                                    isoM_pX = curScan.findMZ(
                                                        curPeakmz + xCount * xOffset,
                                                        ppm,
                                                        start=currentPeakIndex,
                                                    )
                                                    isoM_pX = curScan.getMostIntensePeak(
                                                        isoM_pX[0],
                                                        isoM_pX[1],
                                                        intensityThres,
                                                    )
                                                    if isoM_pX != -1:
                                                        labPeakmz = curScan.mz_list[isoM_pX]
                                                        labPeakIntensity = curScan.intensity_list[isoM_pX]

                                                        # find M'-1 peak
                                                        isoM_pXm1 = curScan.findMZ(
                                                            curPeakmz + xCount * xOffset - cValidationOffset,
                                                            ppm * 2,
                                                            start=currentPeakIndex,
                                                        )
                                                        isoM_pXm1 = curScan.getMostIntensePeak(
                                                            isoM_pXm1[0],
                                                            isoM_pXm1[1],
                                                        )
                                                        if isoM_pXm1 != -1:
                                                            obRatio = curScan.intensity_list[isoM_pXm1] / curScan.intensity_list[isoM_pX]
                                                            if obRatio >= 0.5:
                                                                continue

                                                        # (1.) check if M' and M ratio are as expected
                                                        if checkRatio:
                                                            rat = curPeakIntensity / labPeakIntensity
                                                            if minRatio <= rat <= maxRatio:
                                                                pass  ## ratio check passed
                                                            else:
                                                                continue  ## ratio check not passed

                                                        ## no isotopolog verification needs to be performed
                                                        if peakCountLeft == 1 and peakCountRight == 1:
                                                            curPeakDetectedIonPairs.append(
                                                                mzFeature(
                                                                    mz=curPeakmz,
                                                                    lmz=labPeakmz,
                                                                    tmz=xCount * xOffset,
                                                                    xCount=xCount,
                                                                    scanIndex=curScanIndex,
                                                                    loading=curLoading,
                                                                    nIntensity=curPeakIntensity,
                                                                    lIntensity=labPeakIntensity,
                                                                    ionMode=ionMode,
                                                                )
                                                            )

                                                            skipOtherLoadings = True
                                                            continue

                                                        # find M'+1 peak
                                                        isoM_pXp1 = curScan.findMZ(
                                                            curPeakmz + xCount * xOffset + cValidationOffset,
                                                            ppm,
                                                            start=currentPeakIndex,
                                                        )
                                                        isoM_pXp1 = curScan.getMostIntensePeak(
                                                            isoM_pXp1[0],
                                                            isoM_pXp1[1],
                                                        )

                                                        # calculate the ratio of M+1/M
                                                        isoPeakIntensity = curScan.intensity_list[isoM_p1]
                                                        if peakCountLeft == 1:
                                                            ratioN = None
                                                        elif peakCountLeft > 1 and lowAbundanceIsotopeCutoff and isoPeakIntensity <= isotopologIntensityThres:
                                                            ratioN = None
                                                        else:
                                                            ratioN = isoPeakIntensity / curPeakIntensity

                                                        # calculate the ratio of M'+1/M'
                                                        isoLabPeakIntensity = curScan.intensity_list[isoM_pXp1]
                                                        if peakCountRight == 1:
                                                            ratioL = None
                                                        elif peakCountRight > 1 and lowAbundanceIsotopeCutoff and isoLabPeakIntensity <= isotopologIntensityThres:
                                                            ratioL = None
                                                        else:
                                                            ratioL = isoLabPeakIntensity / labPeakIntensity
                                                        # 2. check if the observed M'+1/M' ratio and the M+1/M ratio are approximately equal
                                                        if (ratioN != None and ratioL != None) and abs(ratioN - ratioL) <= intensityErrorL:
                                                            curPeakDetectedIonPairs.append(
                                                                mzFeature(
                                                                    mz=curPeakmz,
                                                                    lmz=labPeakmz,
                                                                    tmz=xCount * xOffset,
                                                                    xCount=xCount,
                                                                    scanIndex=curScanIndex,
                                                                    loading=curLoading,
                                                                    nIntensity=curPeakIntensity,
                                                                    lIntensity=labPeakIntensity,
                                                                    ionMode=ionMode,
                                                                )
                                                            )

                                                            skipOtherLoadings = True
                                                            continue
                                                        elif lowAbundanceIsotopeCutoff and (ratioN == None or ratioL == None):
                                                            curPeakDetectedIonPairs.append(
                                                                mzFeature(
                                                                    mz=curPeakmz,
                                                                    lmz=labPeakmz,
                                                                    tmz=xCount * xOffset,
                                                                    xCount=xCount,
                                                                    scanIndex=curScanIndex,
                                                                    loading=curLoading,
                                                                    nIntensity=curPeakIntensity,
                                                                    lIntensity=labPeakIntensity,
                                                                    ionMode=ionMode,
                                                                )
                                                            )

                                                            skipOtherLoadings = True
                                                            continue
                                            # endregion

                                        # Isotope pattern validation (useCValidation == 0)
                                        # It is tested, if the expected isotope patterns
                                        # separate for the native and labelled metabolite follow a theoretical pattern
                                        # E.g. native 12C and uniformly / partially 13C-labelled metabolites
                                        # |    <--   Cn   -->    |
                                        # ||                    ||
                                        # |||                  |||
                                        # Necessary mainly for 13C-labelling with mirror-symmetric isotope patterns
                                        # NOTE: - Approach is mainly used for 13C-labelling
                                        if not useDoubleLabelingCombinations and useCIsotopePatternValidation == 0:
                                            # region
                                            # (0.) verify that peak is M and not something else (e.g. M+1, M+1...)
                                            ## TODO improve me. Use seven golden rules or the number of carbon atoms
                                            isoM_m1 = curScan.findMZ(curPeakmz - cValidationOffset, ppm)
                                            isoM_m1 = curScan.getMostIntensePeak(isoM_m1[0], isoM_m1[1])
                                            if isoM_m1 != -1:
                                                obRatio = curScan.intensity_list[isoM_m1] / curPeakIntensity
                                                if obRatio >= 0.5:
                                                    continue

                                            # find M+1 peak
                                            isoM_p1 = curScan.findMZ(
                                                curPeakmz + cValidationOffset,
                                                ppm,
                                                start=currentPeakIndex,
                                            )
                                            isoM_p1 = curScan.getMostIntensePeak(isoM_p1[0], isoM_p1[1])
                                            if isoM_p1 != -1 or peakCountLeft == 1 or lowAbundanceIsotopeCutoff:
                                                # test certain number of labelled carbon atoms

                                                for xCount in sorted(xCounts, reverse=True):
                                                    # stop for impossible carbon atom number
                                                    if xCount > curPeakmz * curLoading / 12:
                                                        continue

                                                    # find corresponding M' peak
                                                    isoM_pX = curScan.findMZ(
                                                        curPeakmz + xCount * cValidationOffset,
                                                        ppm,
                                                        start=currentPeakIndex,
                                                    )
                                                    isoM_pX = curScan.getMostIntensePeak(
                                                        isoM_pX[0],
                                                        isoM_pX[1],
                                                        intensityThres,
                                                    )
                                                    if isoM_pX != -1:
                                                        labPeakmz = curScan.mz_list[isoM_pX]
                                                        labPeakIntensity = curScan.intensity_list[isoM_pX]

                                                        # (0.) verify that peak is M' and not something else (e.g. M'-1, M'-2)
                                                        # only for AllExtract experiments
                                                        adjRatio = 0
                                                        isoM_pXp1 = curScan.findMZ(
                                                            curPeakmz + (xCount + 1) * cValidationOffset,
                                                            ppm,
                                                            start=currentPeakIndex,
                                                        )
                                                        isoM_pXp1 = curScan.getMostIntensePeak(
                                                            isoM_pXp1[0],
                                                            isoM_pXp1[1],
                                                        )
                                                        if isoM_pXp1 != -1:
                                                            adjRatio = curScan.intensity_list[isoM_pXp1] / labPeakIntensity

                                                        if not metabolisationExperiment:
                                                            if adjRatio >= 0.5:
                                                                continue

                                                        # (1.) check if M' and M ratio are as expected
                                                        if checkRatio:
                                                            rat = curPeakIntensity / labPeakIntensity
                                                            if minRatio <= rat <= maxRatio:
                                                                pass  ## ratio check passed
                                                            else:
                                                                continue  ## ratio check not passed

                                                        ## no isotopolog verification needs to be performed
                                                        if peakCountLeft == 1 and peakCountRight == 1:
                                                            curPeakDetectedIonPairs.append(
                                                                mzFeature(
                                                                    mz=curPeakmz,
                                                                    lmz=labPeakmz,
                                                                    tmz=xCount * cValidationOffset,
                                                                    xCount=xCount,
                                                                    scanIndex=curScanIndex,
                                                                    loading=curLoading,
                                                                    nIntensity=curPeakIntensity,
                                                                    lIntensity=labPeakIntensity,
                                                                    ionMode=ionMode,
                                                                )
                                                            )

                                                            skipOtherLoadings = True
                                                            continue

                                                        # find M'-1 peak
                                                        isoM_pXm1 = curScan.findMZ(
                                                            curPeakmz + (xCount - 1) * cValidationOffset,
                                                            ppm,
                                                            start=currentPeakIndex,
                                                        )
                                                        isoM_pXm1 = curScan.getMostIntensePeak(
                                                            isoM_pXm1[0],
                                                            isoM_pXm1[1],
                                                        )
                                                        normRatioL = purityLArray[xCount][1]
                                                        normRatioN = purityNArray[xCount][1]

                                                        # 2. check if the observed M'-1/M' ratio fits the theoretical one
                                                        if peakCountRight == 1:
                                                            pass
                                                        elif isoM_pXm1 == -1:
                                                            if lowAbundanceIsotopeCutoff and labPeakIntensity * normRatioL <= isotopologIntensityThres:
                                                                pass
                                                            else:
                                                                continue
                                                        elif isoM_pXm1 != -1:
                                                            isoM_pXm1_Intensity = curScan.intensity_list[isoM_pXm1]
                                                            observedRatioMp = isoM_pXm1_Intensity / labPeakIntensity
                                                            if abs(normRatioL - observedRatioMp) <= intensityErrorL:
                                                                pass
                                                            elif lowAbundanceIsotopeCutoff and isoM_pXm1_Intensity <= isotopologIntensityThres:
                                                                pass
                                                            else:
                                                                continue
                                                        else:
                                                            continue

                                                        # 3. check if the observed M+1/M ratio fits the theoretical one
                                                        if peakCountLeft == 1:
                                                            pass
                                                        elif isoM_p1 == -1:
                                                            if lowAbundanceIsotopeCutoff and curPeakIntensity * (normRatioN + adjRatio) <= isotopologIntensityThres:
                                                                pass
                                                            else:
                                                                continue
                                                        elif isoM_p1 != -1:
                                                            isoM_p1_Intensity = curScan.intensity_list[isoM_p1]
                                                            observedRatioM = isoM_p1_Intensity / curPeakIntensity
                                                            if abs((normRatioN + adjRatio) - observedRatioM) <= intensityErrorN:
                                                                pass
                                                            elif lowAbundanceIsotopeCutoff and isoM_p1_Intensity <= isotopologIntensityThres:
                                                                pass
                                                            else:
                                                                continue
                                                        else:
                                                            continue

                                                        # All verification criteria are passed, store the ion signal pair
                                                        # for further processing
                                                        curPeakDetectedIonPairs.append(
                                                            mzFeature(
                                                                mz=curPeakmz,
                                                                lmz=labPeakmz,
                                                                tmz=xCount * cValidationOffset,
                                                                xCount=xCount,
                                                                scanIndex=curScanIndex,
                                                                loading=curLoading,
                                                                nIntensity=curPeakIntensity,
                                                                lIntensity=labPeakIntensity,
                                                                ionMode=ionMode,
                                                            )
                                                        )
                                                        skipOtherLoadings = True
                                            # endregion

                                        # labeling patterns derived from one or more labeling-elements (e.g. 13C and D)
                                        # Currently, onle the m/z delta is checked but no isotopolog distribution
                                        # E.g. 13CD3
                                        # |  <-- 13CD3 -->  |
                                        # ||                ||
                                        # |||              ||||
                                        # NOTE: currently, the labeling elements must be defined in the code
                                        # NOTE: The option must be activated and the other two options must be deactivated
                                        if useDoubleLabelingCombinations:
                                            # find M+1 peak
                                            isoM_p1 = curScan.findMZ(
                                                curPeakmz + cValidationOffset / curLoading,
                                                ppm,
                                                start=currentPeakIndex,
                                            )
                                            isoM_p1 = curScan.getMostIntensePeak(isoM_p1[0], isoM_p1[1])
                                            intIsoM_p1 = 0
                                            if isoM_p1 != -1:
                                                intIsoM_p1 = curScan.intensity_list[isoM_p1]

                                            isoM_m1 = curScan.findMZ(
                                                curPeakmz - cValidationOffset / curLoading,
                                                ppm,
                                                start=currentPeakIndex,
                                            )
                                            isoM_m1 = curScan.getMostIntensePeak(isoM_m1[0], isoM_m1[1])
                                            intIsoM_m1 = 0
                                            if isoM_m1 != -1:
                                                intIsoM_m1 = curScan.intensity_list[isoM_m1]

                                            if intIsoM_p1 < curPeakIntensity and intIsoM_m1 < curPeakIntensity and (isoM_p1 != -1 or peakCountLeft == 1 or lowAbundanceIsotopeCutoff):
                                                # test certain number of labelled carbon atoms

                                                for comb in combs:
                                                    # find corresponding M' peak
                                                    isoM_pX = curScan.findMZ(
                                                        curPeakmz + comb.mzdelta / curLoading,
                                                        ppm,
                                                        start=currentPeakIndex,
                                                    )
                                                    isoM_pX = curScan.getMostIntensePeak(
                                                        isoM_pX[0],
                                                        isoM_pX[1],
                                                        intensityThres,
                                                    )
                                                    if isoM_pX != -1:
                                                        labPeakmz = curScan.mz_list[isoM_pX]
                                                        labPeakIntensity = curScan.intensity_list[isoM_pX]

                                                        # (1.) check if M' and M ratio are as expected
                                                        if False:
                                                            rat = curPeakIntensity / labPeakIntensity
                                                            if minRatio <= rat <= maxRatio:
                                                                pass  ## ratio check passed
                                                            else:
                                                                continue  ## ratio check not passed

                                                        # find M'-1 peak
                                                        isoM_pXm1 = curScan._findMZGeneric(
                                                            curPeakmz + (comb.mzdelta - 1.00628 * (1.0 + curPeakmz * ppm / 1000000)) / curLoading,
                                                            curPeakmz + (comb.mzdelta - 1.00335 * (1.0 - curPeakmz * ppm / 1000000)) / curLoading,
                                                        )
                                                        isoM_pXm1 = curScan.getMostIntensePeak(
                                                            isoM_pXm1[0],
                                                            isoM_pXm1[1],
                                                        )

                                                        if isoM_pXm1 != -1:
                                                            isoPeakIntensity = curScan.intensity_list[isoM_pXm1]
                                                            rat = isoPeakIntensity / labPeakIntensity

                                                            if rat <= 0.75:
                                                                pass
                                                            else:
                                                                continue

                                                        # find M'+1 peak
                                                        isoM_pXp1 = curScan._findMZGeneric(
                                                            curPeakmz + (comb.mzdelta + 1.00335 * (1.0 - curPeakmz * ppm / 1000000)) / curLoading,
                                                            curPeakmz + (comb.mzdelta + 1.00628 * (1.0 + curPeakmz * ppm / 1000000)) / curLoading,
                                                        )
                                                        isoM_pXp1 = curScan.getMostIntensePeak(
                                                            isoM_pXp1[0],
                                                            isoM_pXp1[1],
                                                        )

                                                        if isoM_pXp1 != -1:
                                                            isoPeakIntensity = curScan.intensity_list[isoM_pXp1]
                                                            rat = isoPeakIntensity / labPeakIntensity

                                                            if rat <= 0.9:
                                                                pass
                                                            else:
                                                                continue

                                                        # All verification criteria are passed, store the ion signal pair
                                                        # for further processing
                                                        curPeakDetectedIonPairs.append(
                                                            mzFeature(
                                                                mz=curPeakmz,
                                                                lmz=curScan.mz_list[isoM_pX],
                                                                tmz=comb.mzdelta / curLoading,
                                                                xCount=fT.flatToString(comb.atoms),
                                                                scanIndex=curScanIndex,
                                                                loading=curLoading,
                                                                nIntensity=curPeakIntensity,
                                                                lIntensity=labPeakIntensity,
                                                                ionMode=ionMode,
                                                            )
                                                        )

                                                        skipOtherLoadings = True
                                                        # endregion

                            if False:  ## select best fit
                                if len(curPeakDetectedIonPairs) > 0:
                                    bestFit = None
                                    bestFitPPMDiff = 1000000000

                                    ## TODO select best fit based on isotopic pattern (e.g. intensity)

                                    for mt in curPeakDetectedIonPairs:
                                        if abs(mt.lmz - mt.mz - mt.xCount * 1.00335) * 1000000.0 / mt.mz < bestFitPPMDiff:
                                            bestFit = mt
                                            bestFitPPMDiff = abs(mt.lmz - mt.mz - mt.xCount * 1.00335) * 1000000.0 / mt.mz

                                    curScanDetectedIonPairs.append(bestFit)
                            else:  ## use all peak pairs
                                curScanDetectedIonPairs.extend(curPeakDetectedIonPairs)

                if len(curScanDetectedIonPairs) > 0 and False:
                    from .utils import printObjectsAsTable

                    print("\n")
                    print(curScan.retention_time / 60.0)
                    printObjectsAsTable(
                        curScanDetectedIonPairs,
                        attrs=[
                            "mz",
                            "xCount",
                            "loading",
                            "nIntensity",
                            "lIntensity",
                            "ionMode",
                        ],
                    )

                detectedIonPairs.extend(curScanDetectedIonPairs)
                curScanIndex += 1

        except Exception as e:
            import traceback

            traceback.print_exc()

    return detectedIonPairs
