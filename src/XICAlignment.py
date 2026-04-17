"""XICAlignment – Simple peak grouping without R/PTW dependency.

The former PTW-based chromatographic alignment has been removed. This module
now provides a pure-Python fallback that groups peaks by retention-time
proximity using scipy hierarchical clustering, without performing any warping.
"""

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from copy import deepcopy

from .utils import mapArrayToRefTimes


# HELPER METHOD for adding a constant part before and after the actual EIC
def preAppendEICs(
    eics,
    scanTimes,
    pretend=25,
    scanDuration=1,
    mIntMult=1,
    addFront=True,
    addBack=False,
):
    assert len(eics) == len(scanTimes)
    eics = deepcopy(eics)
    scanTimes = deepcopy(scanTimes)

    if pretend == 0:
        return eics, scanTimes

    for j in range(len(eics)):
        eic_arr = np.asarray(eics[j])
        scanTime_arr = np.asarray(scanTimes[j])

        pInt = np.max(eic_arr) * mIntMult

        scanTime_arr = scanTime_arr + scanDuration * pretend
        mScanTime = np.max(scanTime_arr)

        if addFront:
            front_eic = np.full(pretend, pInt)
            front_times = np.arange(pretend) * scanDuration
            eic_arr = np.concatenate([front_eic, eic_arr])
            scanTime_arr = np.concatenate([front_times, scanTime_arr])

        if addBack:
            back_eic = np.full(pretend, pInt)
            back_times = mScanTime + (np.arange(1, pretend + 1) * scanDuration)
            eic_arr = np.concatenate([eic_arr, back_eic])
            scanTime_arr = np.concatenate([scanTime_arr, back_times])

        eics[j] = eic_arr
        scanTimes[j] = scanTime_arr

    return eics, scanTimes


# HELPER METHOD for removing a constant part before and after the actual EIC.
def undoPreAppendEICs(eics, refTimes, pretend=25, scanDuration=1, removeFront=True, removeBack=False):
    eics = deepcopy(eics)
    refTimes = np.asarray(refTimes, dtype=np.float64).copy()

    if pretend == 0:
        return eics, refTimes

    remTime = scanDuration * (pretend)
    mStart = np.argmin(np.abs(refTimes - remTime))

    if removeFront:
        refTimes = refTimes[mStart:]
        refTimes = refTimes - remTime

    if removeBack:
        refTimes = refTimes[:-mStart]

    for j in range(len(eics)):
        eic_arr = np.asarray(eics[j])
        if removeFront:
            eic_arr = eic_arr[mStart:]
        if removeBack:
            eic_arr = eic_arr[:-mStart]
        eics[j] = eic_arr

    return eics, refTimes


class XICAlignment:
    """Peak grouping by retention-time proximity (alignment removed).

    This class replaces the former PTW-based alignment. It groups peaks
    solely by hierarchical clustering on retention times without performing
    any chromatographic warping.
    """

    def __init__(self, scriptLocation=None):
        # scriptLocation kept for API compatibility (no longer used)
        pass

    def alignXIC(
        self,
        eics,
        peakss,
        scantimes,
        align=True,
        nPolynom=3,
        maxTimeDiff=0.36 * 60,
        pretend=100,
        scanDuration=1,
    ):
        """Group chromatographic peaks across samples by RT proximity.

        Returns list of lists; each inner list contains [rt_value, group_id].
        The outer list order matches the flattened order of *peakss*.
        """
        if len(peakss) == 0:
            return []

        if len(peakss) == 1:
            ret = []
            for i in range(len(peakss[0])):
                ret.append([[i, peakss[0][i].NPeakCenterMin]])
            return ret

        # Collect all peak RT values (in minutes) from all samples
        all_rts = []
        sample_indices = []
        for si, peaks in enumerate(peakss):
            for peak in peaks:
                all_rts.append(peak.NPeakCenterMin)
                sample_indices.append(si)

        if len(all_rts) == 0:
            return []

        all_rts = np.array(all_rts, dtype=np.float64)

        if len(all_rts) == 1:
            return [[[all_rts[0], 1]]]

        # Hierarchical clustering on RT values
        rt_matrix = all_rts.reshape(-1, 1)
        Z = linkage(rt_matrix, method="complete", metric="euclidean")
        labels = fcluster(Z, t=maxTimeDiff, criterion="distance")

        # Build result in the same format as the R version
        retld = []
        idx = 0
        for peaks in peakss:
            dd = []
            for peak in peaks:
                dd.append([all_rts[idx], int(labels[idx])])
                idx += 1
            retld.append(dd)

        return retld

    def getAligendXICs(self, eics, scantimes, align=True, nPolynom=3, pretend=100, scanDuration=1):
        """Return EICs mapped to a common time grid (no warping applied).

        Returns (aligned_eics, refTimes).
        """
        if len(eics) == 0:
            return [[]], np.array([])

        eics_copy, scantimes_copy = preAppendEICs(eics, scantimes, pretend=pretend, scanDuration=scanDuration)

        maxTime = np.max([np.max(st) for st in scantimes_copy])
        refTimes = np.arange(int(maxTime))

        aligned = []
        for i in range(len(eics_copy)):
            aligned.append(mapArrayToRefTimes(eics_copy[i], scantimes_copy[i], refTimes))

        aligned, refTimes = undoPreAppendEICs(aligned, refTimes, pretend=pretend, scanDuration=scanDuration)
        return aligned, refTimes
