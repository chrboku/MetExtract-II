if __name__ == "__main__":
    import MExtract

import numpy as np
from scipy.ndimage import median_filter

from .utils import Bunch


class Baseline:
    """Median-window baseline estimation (pure Python replacement for R baseline package).

    Uses a median filter with a configurable half-window size to estimate the
    baseline of a chromatographic signal.
    """

    def __init__(self, hwm: int = 20):
        """
        Args:
            hwm: Half-window size for the median filter.
        """
        self.hwm = hwm

    def getBaseline(self, eic, times) -> np.ndarray:
        """Estimate baseline using a rolling median window.

        Args:
            eic: 1-D intensity array.
            times: 1-D retention time array (unused, kept for API compat).

        Returns:
            1-D numpy array of estimated baseline values.
        """
        eic_arr = np.asarray(eic, dtype=np.float64)
        if eic_arr.size == 0:
            return eic_arr.copy()
        window_size = 2 * self.hwm + 1
        # scipy median_filter handles edge padding internally
        baseline = median_filter(eic_arr, size=window_size, mode="reflect")
        return np.asarray(baseline, dtype=np.float64)


if __name__ == "__main__":
    import matplotlib
    import matplotlib.pyplot as plt

    BL = Baseline()

    from . import Chromatogram
    from .utils import getLastTimeBefore, smoothDataSeries
    from .utils import printObjectsAsTable

    from .chromPeakPicking.GradientPeaks import GradientPeaks

    from .utils import get_main_dir

    from .chromPeakPicking.peakpickers import WaveletTransformPeakPicker, PeakPickerAdapter
    CP = PeakPickerAdapter(WaveletTransformPeakPicker())

    chromatogram = Chromatogram.Chromatogram()

    if False:  # Disabled hardcoded test
        chromatogram.parse_file("C:/PyMetExtract/implTest/D-GEN_Warth/GEN_24h_5uM_supern_3_neg_079.mzXML")
        mz = 370.90112
        ppm = 10.0
        cn = 4
        z = 1
        scales = [2, 7]
        dmz = 1.00628

        scanEvents = chromatogram.getFilterLines()
        for s in scanEvents:
            print(s)

        scanEvent = sorted(scanEvents)[0]
        print("selected scan event:", scanEvent)

        fig = plt.figure(facecolor="white")
        ax = fig.add_subplot(111)

        eic, times, timesI, mzs = chromatogram.getEIC(mz, ppm, scanEvent)
        eicL, times, timesI, mzsL = chromatogram.getEIC(mz + cn * dmz / z, ppm, scanEvent)

        ax.plot(times / 60.0, eic)

        eicBL = BL.getBaseline(eic, times)

        ax.plot(times / 60.0, eicBL)

        ax.plot(times / 60.0, -np.maximum(0, eic - eicBL))

        ret = CP.getPeaksFor(times, eic)
        for peak in ret:
            peak.peakAtTime = times[peak.peakIndex] / 60.0
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

        ret = CP.getPeaksFor(
            times,
            [max(0, eic[i] - eicBL[i]) for i in range(len(eic))],
        )
        for peak in ret:
            peak.peakAtTime = times[peak.peakIndex] / 60.0
        print("Native without baseline")
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

        plt.show()
