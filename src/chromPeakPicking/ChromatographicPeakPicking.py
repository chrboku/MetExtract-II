from .GradientPeaks import GradientPeaks
from .peakpickers import (
    BasePeakPicker,
    ChromatographicPeak,
    GradientDescentPeakPicker,
    GMMPeakPicker,
    WaveletTransformPeakPicker,
    SecondDerivativePeakPicker,
    MatchedFilterPeakPicker,
    PeakPickerAdapter,
)

# Supported algorithm names (case-insensitive)
_ALGORITHMS: dict[str, type[BasePeakPicker]] = {
    "gradientdescent": GradientDescentPeakPicker,
    "gmm": GMMPeakPicker,
    "wavelettransform": WaveletTransformPeakPicker,
    "secondderivative": SecondDerivativePeakPicker,
    "matchedfilter": MatchedFilterPeakPicker,
}

# Legacy aliases
_ALIASES: dict[str, str] = {
    "massspecwavelet": "wavelettransform",
    "gradientpeaks": "gradientdescent",
}


class ChromPeakPicking:
    """Factory for chromatographic peak picking algorithms.

    Supports both the new BasePeakPicker implementations and the legacy
    GradientPeaks class.
    """

    method = None

    def __init__(self, algorithm: str = "GradientDescent", **kwargs):
        self.algorithm = algorithm
        key = algorithm.lower().replace("_", "").replace(" ", "")
        key = _ALIASES.get(key, key)

        if key == "gradientpeaks_legacy":
            self._picker = GradientPeaks(**kwargs)
        elif key in _ALGORITHMS:
            self._picker = PeakPickerAdapter(_ALGORITHMS[key](**kwargs))
        else:
            raise ValueError(
                f"Unknown peak picking algorithm '{algorithm}'. "
                f"Available: {', '.join(list(_ALGORITHMS.keys()) + ['GradientPeaks_legacy'])}"
            )

    def getPeaksFor(self, times, eic, **kwargs):
        """Legacy-compatible interface."""
        return self._picker.getPeaksFor(times, eic, **kwargs)

    @staticmethod
    def getCurrentPeakPicking(*argv, **kwargs):
        if ChromPeakPicking.method is None:
            ChromPeakPicking.method = ChromPeakPicking(*argv, **kwargs)
        return ChromPeakPicking.method
