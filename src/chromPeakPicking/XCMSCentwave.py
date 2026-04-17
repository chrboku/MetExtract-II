"""XCMSCentwave-compatible interface using native Python peak picking.

This module replaces the former rpy2/R-based XCMS centWave implementation
with the native MatchedFilterPeakPicker while keeping the same public API.
"""

from ..utils import Bunch
from .peakpickers import MatchedFilterPeakPicker, PeakPickerAdapter
import numpy as np


class XCMSCentwave:
    """Drop-in replacement for the R-based XCMSCentwave class."""

    def __init__(self, scriptLocation=None):
        # scriptLocation kept for API compatibility but no longer used
        self.scriptLocation = scriptLocation
        self._picker = MatchedFilterPeakPicker(expected_peak_width=15)
        self._adapter = PeakPickerAdapter(self._picker)

    def cleanUp(self):
        """No-op for backward compatibility."""
        pass

    def getPeaksFor(self, times, eic, scales=None):
        """Detect peaks – returns list of Bunch objects."""
        if scales is not None and len(scales) >= 2:
            self._picker.expected_peak_width = max(3, (scales[0] + scales[1]) // 2)
        return self._adapter.getPeaksFor(times, eic, scales=scales)
