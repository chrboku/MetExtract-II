"""MassSpecWavelet-compatible interface using pure-Python CWT peak picking.

This module provides the MassSpecWavelet peak picking interface using
the native WaveletTransformPeakPicker while keeping the same public API.
"""

from .peakpickers import PeakPickerAdapter, WaveletTransformPeakPicker


class MassSpecWavelet:
    """Drop-in replacement for the R-based MassSpecWavelet class.

    Uses :class:`WaveletTransformPeakPicker` under the hood but exposes
    the same ``getPeaksFor()`` interface expected by the rest of MetExtract.
    """

    def __init__(self, scriptLocation=None, scales=None, snrTh=None, minScans=None):
        # scriptLocation kept for API compatibility but no longer used
        self.scriptLocation = scriptLocation

        self.scales = scales if scales is not None else [11, 66]
        self.snrTh = snrTh if snrTh is not None else 0.1
        self.minScans = minScans if minScans is not None else 3

        self._picker = WaveletTransformPeakPicker(
            min_scale=max(1, self.scales[0]),
            max_scale=max(self.scales[1], self.scales[0] + 1),
            min_ridge_length=max(self.minScans, 1),
        )
        self._adapter = PeakPickerAdapter(self._picker)

    def cleanUp(self):
        """No-op for backward compatibility."""
        pass

    def getPeaksFor(self, timesi, eici, scales=None, snrTh=None, startIndex=None, endIndex=None):
        """Detect peaks using CWT – returns list of Bunch objects."""
        return self._adapter.getPeaksFor(
            timesi,
            eici,
            scales=scales,
            snrTh=snrTh,
            startIndex=startIndex,
            endIndex=endIndex,
        )
