"""
Peak Picking Architecture for LC-HRMS Extracted Ion Chromatograms (EICs).

This module provides a robust, object-oriented peak picking framework with
five distinct algorithms, replacing the previous R/rpy2-based implementation.

Classes:
    ChromatographicPeak - Dataclass for peak results
    BasePeakPicker - Abstract base class with shared utilities
    GradientDescentPeakPicker - Gradient-based peak boundary detection
    GMMPeakPicker - Gaussian Mixture Model / curve fitting approach
    WaveletTransformPeakPicker - CWT-based (MassSpecWavelet-inspired) detection
    SecondDerivativePeakPicker - Savitzky-Golay second derivative approach
    MatchedFilterPeakPicker - Cross-correlation with template peak shape
"""

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Optional, Sequence

import numpy as np
from scipy import signal, optimize, stats
from scipy.ndimage import gaussian_filter1d, uniform_filter1d


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ChromatographicPeak:
    """Stores the output of chromatographic peak detection.

    Attributes:
        start_rt: Retention time at the left boundary of the peak (seconds).
        end_rt: Retention time at the right boundary of the peak (seconds).
        apex_rt: Retention time of the peak apex (seconds).
        fwhm: Full Width at Half Maximum (seconds).
        snr: Signal-to-Noise Ratio.
        area: Integrated peak area (trapezoidal rule).
        apex_intensity: Intensity at the peak apex.
        start_index: Array index of the left boundary.
        end_index: Array index of the right boundary.
        apex_index: Array index of the peak apex.
        baseline: Estimated local baseline intensity.
    """

    start_rt: float
    end_rt: float
    apex_rt: float
    fwhm: float
    snr: float
    area: float
    apex_intensity: float = 0.0
    start_index: int = 0
    end_index: int = 0
    apex_index: int = 0
    baseline: float = 0.0


class SmoothingMethod(str, Enum):
    """Supported smoothing kernels."""

    MOVING_AVERAGE = "moving_average"
    GAUSSIAN = "gaussian"
    SAVITZKY_GOLAY = "savitzky_golay"


class DistributionType(str, Enum):
    """Distribution types for curve fitting."""

    GAUSSIAN = "gaussian"
    EMG = "emg"  # Exponentially Modified Gaussian


# ---------------------------------------------------------------------------
# Base class
# ---------------------------------------------------------------------------

class BasePeakPicker(ABC):
    """Abstract base class for all peak picking algorithms."""

    @abstractmethod
    def pick_peaks(
        self,
        rt: np.ndarray,
        intensity: np.ndarray,
        **kwargs,
    ) -> list[ChromatographicPeak]:
        """Detect chromatographic peaks in an EIC.

        Args:
            rt: 1-D array of retention times (seconds).
            intensity: 1-D array of corresponding intensities.

        Returns:
            List of detected ChromatographicPeak objects.
        """
        ...

    # -- Shared utilities ---------------------------------------------------

    @staticmethod
    def _validate_inputs(rt: np.ndarray, intensity: np.ndarray) -> bool:
        """Return *True* if the inputs are usable, *False* otherwise."""
        if rt is None or intensity is None:
            return False
        rt = np.asarray(rt, dtype=np.float64)
        intensity = np.asarray(intensity, dtype=np.float64)
        if rt.size == 0 or intensity.size == 0:
            return False
        if rt.shape != intensity.shape:
            return False
        if np.all(intensity == intensity[0]):
            return False  # flat signal
        return True

    @staticmethod
    def compute_area(rt: np.ndarray, intensity: np.ndarray, start: int, end: int) -> float:
        """Trapezoidal integration between *start* and *end* indices."""
        if end <= start:
            return 0.0
        _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz", None)
        if _trapz is None:
            # Manual trapezoidal rule fallback
            x = rt[start : end + 1]
            y = intensity[start : end + 1]
            return float(np.sum((y[:-1] + y[1:]) * np.diff(x) / 2.0))
        return float(_trapz(intensity[start : end + 1], rt[start : end + 1]))

    @staticmethod
    def compute_snr(
        intensity: np.ndarray,
        apex_index: int,
        start_index: int,
        end_index: int,
    ) -> float:
        """Compute signal-to-noise ratio.

        Noise is estimated as the median of the intensities outside the peak
        region falling back to the minimum of the flanks.
        """
        apex_val = intensity[apex_index]
        flank_min = min(intensity[start_index], intensity[end_index])
        # noise from outside peak
        mask = np.ones(len(intensity), dtype=bool)
        mask[start_index : end_index + 1] = False
        outside = intensity[mask]
        if outside.size > 0:
            noise = float(np.median(outside))
        else:
            noise = flank_min
        if noise <= 0:
            noise = flank_min if flank_min > 0 else 1.0
        return float((apex_val - flank_min) / noise)

    @staticmethod
    def compute_fwhm(
        rt: np.ndarray,
        intensity: np.ndarray,
        apex_index: int,
        start_index: int,
        end_index: int,
    ) -> float:
        """Estimate FWHM by linear interpolation."""
        apex_val = intensity[apex_index]
        baseline = min(intensity[start_index], intensity[end_index])
        half_max = (apex_val + baseline) / 2.0

        # Left side
        left_rt = rt[start_index]
        for i in range(apex_index, start_index - 1, -1):
            if intensity[i] <= half_max:
                if i < apex_index:
                    frac = (half_max - intensity[i]) / max(intensity[i + 1] - intensity[i], 1e-12)
                    left_rt = rt[i] + frac * (rt[i + 1] - rt[i])
                break

        # Right side
        right_rt = rt[end_index]
        for i in range(apex_index, end_index + 1):
            if intensity[i] <= half_max:
                if i > apex_index:
                    frac = (half_max - intensity[i]) / max(intensity[i - 1] - intensity[i], 1e-12)
                    right_rt = rt[i] - frac * (rt[i] - rt[i - 1])
                break

        return max(float(right_rt - left_rt), 0.0)

    @staticmethod
    def smooth_eic(
        intensity: np.ndarray,
        method: SmoothingMethod = SmoothingMethod.GAUSSIAN,
        window_size: int = 5,
        polynom_order: int = 3,
        iterations: int = 1,
    ) -> np.ndarray:
        """Apply *iterations* rounds of smoothing to *intensity*.

        Args:
            intensity: Raw intensity array.
            method: Smoothing kernel type.
            window_size: Kernel width (must be odd for Savitzky-Golay).
            polynom_order: Polynomial order (Savitzky-Golay only).
            iterations: Number of sequential smoothing passes.

        Returns:
            Smoothed intensity array.
        """
        smoothed = np.array(intensity, dtype=np.float64)
        for _ in range(max(1, iterations)):
            if method == SmoothingMethod.MOVING_AVERAGE:
                smoothed = uniform_filter1d(smoothed, size=max(window_size, 1))
            elif method == SmoothingMethod.GAUSSIAN:
                sigma = max(window_size / 4.0, 0.5)
                smoothed = gaussian_filter1d(smoothed, sigma=sigma)
            elif method == SmoothingMethod.SAVITZKY_GOLAY:
                win = max(window_size, polynom_order + 2)
                if win % 2 == 0:
                    win += 1
                smoothed = signal.savgol_filter(smoothed, window_length=win, polyorder=polynom_order)
            smoothed = np.maximum(smoothed, 0.0)
        return smoothed

    def _build_peak(
        self,
        rt: np.ndarray,
        intensity: np.ndarray,
        apex_index: int,
        start_index: int,
        end_index: int,
    ) -> ChromatographicPeak:
        """Convenience factory that populates all common metrics."""
        area = self.compute_area(rt, intensity, start_index, end_index)
        fwhm = self.compute_fwhm(rt, intensity, apex_index, start_index, end_index)
        snr = self.compute_snr(intensity, apex_index, start_index, end_index)
        baseline = float(min(intensity[start_index], intensity[end_index]))
        return ChromatographicPeak(
            start_rt=float(rt[start_index]),
            end_rt=float(rt[end_index]),
            apex_rt=float(rt[apex_index]),
            fwhm=fwhm,
            snr=snr,
            area=area,
            apex_intensity=float(intensity[apex_index]),
            start_index=int(start_index),
            end_index=int(end_index),
            apex_index=int(apex_index),
            baseline=baseline,
        )


# ---------------------------------------------------------------------------
# Method 1: GradientDescentPeakPicker
# ---------------------------------------------------------------------------

class GradientDescentPeakPicker(BasePeakPicker):
    """Smoothes the EIC, finds local maxima and descends slopes to define bounds.

    Parameters:
        smoothing_method: Kernel used to smooth the EIC before detection.
        smoothing_window: Kernel size for smoothing.
        smoothing_iterations: Number of sequential smoothing passes.
        smoothing_polynom: Polynomial order (SG only).
        min_intensity: Minimum apex intensity for a peak to be kept.
        min_flank_intensity: Minimum intensity at the flank to continue descending.
        min_increase_ratio: Ratio threshold to stop descending.
        consecutive_scans: Number of consecutive scans to check for intensity change.
        intensity_change_threshold: Stop when change over *consecutive_scans* < this.
        min_width: Minimum peak width in scans.
        max_width: Maximum peak width in scans.
    """

    def __init__(
        self,
        smoothing_method: SmoothingMethod = SmoothingMethod.GAUSSIAN,
        smoothing_window: int = 5,
        smoothing_iterations: int = 1,
        smoothing_polynom: int = 3,
        min_intensity: float = 1000.0,
        min_flank_intensity: float = 10.0,
        min_increase_ratio: float = 0.05,
        consecutive_scans: int = 3,
        intensity_change_threshold: float = 0.0,
        min_width: int = 3,
        max_width: int = 300,
    ):
        self.smoothing_method = smoothing_method
        self.smoothing_window = smoothing_window
        self.smoothing_iterations = smoothing_iterations
        self.smoothing_polynom = smoothing_polynom
        self.min_intensity = min_intensity
        self.min_flank_intensity = min_flank_intensity
        self.min_increase_ratio = min_increase_ratio
        self.consecutive_scans = max(consecutive_scans, 1)
        self.intensity_change_threshold = intensity_change_threshold
        self.min_width = min_width
        self.max_width = max_width

    # -- helpers ------------------------------------------------------------

    def _find_local_maxima(self, y: np.ndarray) -> list[int]:
        """Return indices of local maxima above *min_intensity*."""
        if y.size < 3:
            return []
        maxima: list[int] = []
        if y[0] >= self.min_intensity and y[0] > y[1]:
            maxima.append(0)
        mid = (y[1:-1] >= self.min_intensity) & (y[:-2] < y[1:-1]) & (y[1:-1] > y[2:])
        maxima.extend((np.where(mid)[0] + 1).tolist())
        if y[-1] >= self.min_intensity and y[-2] <= y[-1]:
            maxima.append(len(y) - 1)
        return maxima

    def _descend_left(self, y: np.ndarray, apex: int) -> int:
        """Descend from *apex* to the left until the signal rises or stalls."""
        pos = apex
        rising_count = 0
        while pos > 0:
            ratio = y[pos - 1] / max(y[pos], 1e-12)
            # signal is rising again
            if ratio > 1.0:
                rising_count += 1
                if rising_count >= self.consecutive_scans:
                    break
            else:
                rising_count = 0

            # intensity change stall check
            window_start = max(0, pos - self.consecutive_scans)
            if (y[pos] - y[window_start]) < self.intensity_change_threshold and pos < apex:
                break

            if y[pos - 1] < self.min_flank_intensity:
                break
            if ratio < self.min_increase_ratio:
                break
            pos -= 1
        return pos

    def _descend_right(self, y: np.ndarray, apex: int) -> int:
        """Descend from *apex* to the right."""
        pos = apex
        rising_count = 0
        n = len(y)
        while pos < n - 1:
            ratio = y[pos + 1] / max(y[pos], 1e-12)
            if ratio > 1.0:
                rising_count += 1
                if rising_count >= self.consecutive_scans:
                    break
            else:
                rising_count = 0

            window_end = min(n - 1, pos + self.consecutive_scans)
            if (y[pos] - y[window_end]) < self.intensity_change_threshold and pos > apex:
                break

            if y[pos + 1] < self.min_flank_intensity:
                break
            if ratio < self.min_increase_ratio:
                break
            pos += 1
        return pos

    # -- public API ---------------------------------------------------------

    def pick_peaks(
        self,
        rt: np.ndarray,
        intensity: np.ndarray,
        **kwargs,
    ) -> list[ChromatographicPeak]:
        if not self._validate_inputs(rt, intensity):
            return []
        rt = np.asarray(rt, dtype=np.float64)
        intensity = np.asarray(intensity, dtype=np.float64)

        smoothed = self.smooth_eic(
            intensity,
            method=self.smoothing_method,
            window_size=self.smoothing_window,
            polynom_order=self.smoothing_polynom,
            iterations=self.smoothing_iterations,
        )

        maxima = self._find_local_maxima(smoothed)
        peaks: list[ChromatographicPeak] = []
        for apex in maxima:
            left = self._descend_left(smoothed, apex)
            right = self._descend_right(smoothed, apex)
            width = right - left
            if width < self.min_width or width > self.max_width:
                continue
            peaks.append(self._build_peak(rt, intensity, apex, left, right))
        return peaks


# ---------------------------------------------------------------------------
# Method 2: GMMPeakPicker
# ---------------------------------------------------------------------------

def _gaussian(x: np.ndarray, amplitude: float, mean: float, sigma: float) -> np.ndarray:
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def _emg(x: np.ndarray, amplitude: float, mean: float, sigma: float, lam: float) -> np.ndarray:
    """Exponentially Modified Gaussian."""
    lam = max(lam, 1e-12)
    part1 = (lam / 2.0) * np.exp((lam / 2.0) * (2.0 * mean + lam * sigma**2 - 2.0 * x))
    with np.errstate(over="ignore", invalid="ignore"):
        arg = (mean + lam * sigma**2 - x) / (np.sqrt(2.0) * sigma)
        from scipy.special import erfc
        part2 = erfc(arg)
    result = amplitude * part1 * part2
    return np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)


class GMMPeakPicker(BasePeakPicker):
    """Chunks the EIC into regions of interest and fits overlapping distributions.

    Parameters:
        baseline_fraction: Fraction of median intensity used as baseline threshold.
        min_chunk_size: Minimum number of scans in a region of interest.
        max_peaks_per_chunk: Maximum number of peaks fitted per chunk.
        distribution: Distribution type for curve fitting.
        fit_max_iterations: Maximum iterations for the curve fitter.
    """

    def __init__(
        self,
        baseline_fraction: float = 0.1,
        min_chunk_size: int = 5,
        max_peaks_per_chunk: int = 10,
        distribution: DistributionType = DistributionType.GAUSSIAN,
        fit_max_iterations: int = 5000,
    ):
        self.baseline_fraction = baseline_fraction
        self.min_chunk_size = min_chunk_size
        self.max_peaks_per_chunk = max_peaks_per_chunk
        self.distribution = distribution
        self.fit_max_iterations = fit_max_iterations

    # -- helpers ------------------------------------------------------------

    @staticmethod
    def _find_chunks(mask: np.ndarray, min_size: int) -> list[tuple[int, int]]:
        """Return (start, end) inclusive index pairs for contiguous True regions."""
        diffs = np.diff(mask.astype(np.int8))
        starts = np.where(diffs == 1)[0] + 1
        ends = np.where(diffs == -1)[0]
        if mask[0]:
            starts = np.concatenate(([0], starts))
        if mask[-1]:
            ends = np.concatenate((ends, [len(mask) - 1]))
        chunks = [(int(s), int(e)) for s, e in zip(starts, ends) if (e - s + 1) >= min_size]
        return chunks

    def _estimate_num_peaks(self, y: np.ndarray) -> int:
        """Rough local-maxima count."""
        if y.size < 3:
            return 1
        mid = (y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])
        return max(int(np.sum(mid)), 1)

    def _fit_chunk(
        self,
        rt_chunk: np.ndarray,
        int_chunk: np.ndarray,
    ) -> list[tuple[float, float, float]]:
        """Fit overlapping distributions, return list of (mean, sigma, amplitude)."""
        n_peaks = min(self._estimate_num_peaks(int_chunk), self.max_peaks_per_chunk)

        # Build initial guesses by finding n_peaks highest local maxima
        if int_chunk.size < 3:
            return [(float(rt_chunk[np.argmax(int_chunk)]), float(rt_chunk[-1] - rt_chunk[0]) / 4, float(np.max(int_chunk)))]

        local_max_mask = np.zeros(len(int_chunk), dtype=bool)
        local_max_mask[1:-1] = (int_chunk[1:-1] > int_chunk[:-2]) & (int_chunk[1:-1] > int_chunk[2:])
        local_max_idx = np.where(local_max_mask)[0]
        if len(local_max_idx) == 0:
            local_max_idx = np.array([np.argmax(int_chunk)])
        sorted_idx = local_max_idx[np.argsort(int_chunk[local_max_idx])[::-1]]
        chosen = sorted_idx[: n_peaks]

        rt_range = rt_chunk[-1] - rt_chunk[0]
        sigma_guess = max(rt_range / (2.0 * n_peaks), 0.1)

        p0: list[float] = []
        lower: list[float] = []
        upper: list[float] = []
        for idx in chosen:
            p0.extend([float(int_chunk[idx]), float(rt_chunk[idx]), sigma_guess])
            lower.extend([0, float(rt_chunk[0]), 1e-3])
            upper.extend([float(np.max(int_chunk)) * 3, float(rt_chunk[-1]), rt_range])
            if self.distribution == DistributionType.EMG:
                p0.append(1.0)
                lower.append(1e-6)
                upper.append(100.0)

        if self.distribution == DistributionType.GAUSSIAN:
            n_params = 3

            def model(x: np.ndarray, *params: float) -> np.ndarray:
                result = np.zeros_like(x)
                for i in range(0, len(params), n_params):
                    result += _gaussian(x, params[i], params[i + 1], params[i + 2])
                return result
        else:
            n_params = 4

            def model(x: np.ndarray, *params: float) -> np.ndarray:
                result = np.zeros_like(x)
                for i in range(0, len(params), n_params):
                    result += _emg(x, params[i], params[i + 1], params[i + 2], params[i + 3])
                return result

        try:
            popt, _ = optimize.curve_fit(
                model,
                rt_chunk,
                int_chunk,
                p0=p0,
                bounds=(lower, upper),
                maxfev=self.fit_max_iterations,
            )
        except (RuntimeError, ValueError):
            # Fallback: return single peak at maximum
            amax = int(np.argmax(int_chunk))
            return [(float(rt_chunk[amax]), sigma_guess, float(int_chunk[amax]))]

        results: list[tuple[float, float, float]] = []
        for i in range(0, len(popt), n_params):
            amp, mu, sig = popt[i], popt[i + 1], popt[i + 2]
            results.append((mu, sig, amp))
        return results

    # -- public API ---------------------------------------------------------

    def pick_peaks(
        self,
        rt: np.ndarray,
        intensity: np.ndarray,
        **kwargs,
    ) -> list[ChromatographicPeak]:
        if not self._validate_inputs(rt, intensity):
            return []
        rt = np.asarray(rt, dtype=np.float64)
        intensity = np.asarray(intensity, dtype=np.float64)

        # Determine baseline threshold
        med = float(np.median(intensity))
        threshold = med + self.baseline_fraction * med

        above = intensity > threshold
        chunks = self._find_chunks(above, self.min_chunk_size)

        peaks: list[ChromatographicPeak] = []
        for s, e in chunks:
            rt_c = rt[s : e + 1]
            int_c = intensity[s : e + 1]
            fits = self._fit_chunk(rt_c, int_c)
            for mu, sig, amp in fits:
                # Define bounds at +/- 3 sigma clipped to chunk
                left_rt = max(mu - 3 * sig, rt_c[0])
                right_rt = min(mu + 3 * sig, rt_c[-1])
                left_idx = s + int(np.searchsorted(rt_c, left_rt))
                right_idx = s + int(np.searchsorted(rt_c, right_rt, side="right") - 1)
                left_idx = max(left_idx, s)
                right_idx = min(right_idx, e)
                apex_idx = left_idx + int(np.argmax(intensity[left_idx : right_idx + 1]))
                peaks.append(self._build_peak(rt, intensity, apex_idx, left_idx, right_idx))
        return peaks


# ---------------------------------------------------------------------------
# Method 3: WaveletTransformPeakPicker
# ---------------------------------------------------------------------------

class WaveletTransformPeakPicker(BasePeakPicker):
    """CWT-based peak detection inspired by MassSpecWavelet.

    Uses a Ricker (Mexican hat) wavelet and 2-D ridge finding across scales.

    Parameters:
        min_scale: Minimum wavelet scale in number of scans.
        max_scale: Maximum wavelet scale in number of scans.
        num_scales: Number of scales to evaluate.
        snr_threshold: Minimum SNR for accepting a ridge.
        min_ridge_length: Minimum number of scales a ridge must span.
        gap_threshold: Max consecutive scale gaps allowed in a ridge.
    """

    def __init__(
        self,
        min_scale: int = 1,
        max_scale: int = 64,
        num_scales: int = 32,
        snr_threshold: float = 1.0,
        min_ridge_length: int = 4,
        gap_threshold: int = 2,
    ):
        self.min_scale = max(min_scale, 1)
        self.max_scale = max_scale
        self.num_scales = num_scales
        self.snr_threshold = snr_threshold
        self.min_ridge_length = max(min_ridge_length, 1)
        self.gap_threshold = gap_threshold

    # -- helpers ------------------------------------------------------------

    @staticmethod
    def _ricker_wavelet(points: int, width: float) -> np.ndarray:
        """Generate a Ricker (Mexican hat) wavelet."""
        x = np.arange(points) - (points - 1) / 2
        sigma = width
        norm = 2.0 / (np.sqrt(3.0 * sigma) * (np.pi ** 0.25))
        wsq = sigma ** 2
        vec = norm * (1.0 - (x ** 2) / wsq) * np.exp(-(x ** 2) / (2.0 * wsq))
        return vec

    def _cwt_matrix(self, intensity: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Compute CWT coefficient matrix and corresponding widths."""
        widths = np.linspace(self.min_scale, self.max_scale, self.num_scales)
        n = len(intensity)
        cwt_mat = np.empty((len(widths), n), dtype=np.float64)
        for i, w in enumerate(widths):
            wavelet_len = min(10 * int(np.ceil(w)), n)
            if wavelet_len % 2 == 0:
                wavelet_len += 1
            wavelet = self._ricker_wavelet(wavelet_len, w)
            cwt_mat[i] = np.convolve(intensity, wavelet, mode="same")
        return cwt_mat, widths

    def _find_ridges(
        self,
        cwt_mat: np.ndarray,
        widths: np.ndarray,
    ) -> list[list[tuple[int, int]]]:
        """Simple greedy ridge-line tracing across scales.

        Returns list of ridges; each ridge is a list of (scale_index, position) tuples.
        """
        n_scales, n_points = cwt_mat.shape
        # Start from the finest scale and trace upward
        # Active ridges: list of (list of (scale_idx, pos), last_pos, gap_count)
        active: list[tuple[list[tuple[int, int]], int, int]] = []
        finished: list[list[tuple[int, int]]] = []

        for si in range(n_scales):
            row = cwt_mat[si]
            # Find local maxima in this scale
            max_indices: list[int] = []
            if n_points >= 3:
                mid = (row[1:-1] > row[:-2]) & (row[1:-1] > row[2:]) & (row[1:-1] > 0)
                max_indices = (np.where(mid)[0] + 1).tolist()
            if n_points >= 1:
                if row[0] > 0 and (n_points < 2 or row[0] > row[1]):
                    max_indices.insert(0, 0)
                if n_points >= 2 and row[-1] > 0 and row[-1] > row[-2]:
                    max_indices.append(n_points - 1)

            used = set()
            new_active: list[tuple[list[tuple[int, int]], int, int]] = []
            tolerance = max(int(widths[si]), 2)

            for ridge, last_pos, gap in active:
                best_idx = None
                best_dist = tolerance + 1
                for mi in max_indices:
                    if mi in used:
                        continue
                    d = abs(mi - last_pos)
                    if d < best_dist:
                        best_dist = d
                        best_idx = mi
                if best_idx is not None and best_dist <= tolerance:
                    ridge.append((si, best_idx))
                    used.add(best_idx)
                    new_active.append((ridge, best_idx, 0))
                else:
                    if gap + 1 <= self.gap_threshold:
                        new_active.append((ridge, last_pos, gap + 1))
                    else:
                        finished.append(ridge)

            # Start new ridges from unused maxima
            for mi in max_indices:
                if mi not in used:
                    new_active.append(([(si, mi)], mi, 0))

            active = new_active

        finished.extend([r for r, _, _ in active])
        return finished

    def _peak_bounds_from_cwt(
        self,
        intensity: np.ndarray,
        cwt_mat: np.ndarray,
        scale_idx: int,
        position: int,
    ) -> tuple[int, int]:
        """Find zero-crossings at the given scale to set peak boundaries."""
        row = cwt_mat[scale_idx]
        n = len(row)
        left = position
        while left > 0 and row[left] > 0:
            left -= 1
        right = position
        while right < n - 1 and row[right] > 0:
            right += 1
        return left, right

    # -- public API ---------------------------------------------------------

    def pick_peaks(
        self,
        rt: np.ndarray,
        intensity: np.ndarray,
        **kwargs,
    ) -> list[ChromatographicPeak]:
        if not self._validate_inputs(rt, intensity):
            return []
        rt = np.asarray(rt, dtype=np.float64)
        intensity = np.asarray(intensity, dtype=np.float64)

        cwt_mat, widths = self._cwt_matrix(intensity)
        ridges = self._find_ridges(cwt_mat, widths)

        peaks: list[ChromatographicPeak] = []
        used_positions: set[int] = set()
        for ridge in ridges:
            if len(ridge) < self.min_ridge_length:
                continue
            # Use the best (highest CWT coeff) point in the ridge as the apex
            best_scale_idx, best_pos = max(ridge, key=lambda sp: cwt_mat[sp[0], sp[1]])
            if best_pos in used_positions:
                continue
            # Check SNR at the best scale
            row = cwt_mat[best_scale_idx]
            noise_est = float(np.median(np.abs(row)))
            if noise_est > 0 and cwt_mat[best_scale_idx, best_pos] / noise_est < self.snr_threshold:
                continue
            left, right = self._peak_bounds_from_cwt(intensity, cwt_mat, best_scale_idx, best_pos)
            # Refine apex as the argmax in the original intensity within bounds
            apex = left + int(np.argmax(intensity[left : right + 1]))
            used_positions.add(apex)
            peaks.append(self._build_peak(rt, intensity, apex, left, right))
        return peaks


# ---------------------------------------------------------------------------
# Method 4: SecondDerivativePeakPicker
# ---------------------------------------------------------------------------

class SecondDerivativePeakPicker(BasePeakPicker):
    """Uses the Savitzky-Golay second derivative to find peaks.

    Local minima of the second derivative correspond to peak apexes;
    zero-crossings (negative→positive) mark inflection points → peak bounds.

    Parameters:
        sg_window: Window length for the Savitzky-Golay filter.
        sg_polyorder: Polynomial order for the SG filter.
        pre_smoothing: Optional pre-smoothing methods/iterations.
        pre_smoothing_window: Window size for pre-smoothing.
        pre_smoothing_iterations: Number of pre-smoothing passes.
        min_apex_intensity: Minimum raw intensity at peak apex.
    """

    def __init__(
        self,
        sg_window: int = 11,
        sg_polyorder: int = 3,
        pre_smoothing: Optional[SmoothingMethod] = None,
        pre_smoothing_window: int = 5,
        pre_smoothing_iterations: int = 1,
        min_apex_intensity: float = 0.0,
    ):
        self.sg_window = sg_window if sg_window % 2 == 1 else sg_window + 1
        self.sg_polyorder = sg_polyorder
        self.pre_smoothing = pre_smoothing
        self.pre_smoothing_window = pre_smoothing_window
        self.pre_smoothing_iterations = pre_smoothing_iterations
        self.min_apex_intensity = min_apex_intensity

    # -- public API ---------------------------------------------------------

    def pick_peaks(
        self,
        rt: np.ndarray,
        intensity: np.ndarray,
        **kwargs,
    ) -> list[ChromatographicPeak]:
        if not self._validate_inputs(rt, intensity):
            return []
        rt = np.asarray(rt, dtype=np.float64)
        intensity = np.asarray(intensity, dtype=np.float64)

        # Optional pre-smoothing
        smoothed = intensity.copy()
        if self.pre_smoothing is not None:
            smoothed = self.smooth_eic(
                smoothed,
                method=self.pre_smoothing,
                window_size=self.pre_smoothing_window,
                iterations=self.pre_smoothing_iterations,
            )

        n = len(smoothed)
        if n < self.sg_window:
            return []

        # Compute second derivative via Savitzky-Golay
        d2 = signal.savgol_filter(smoothed, window_length=self.sg_window, polyorder=self.sg_polyorder, deriv=2)

        # Find minima of d2 (= peak apexes in original signal)
        apex_candidates: list[int] = []
        if n >= 3:
            mid = (d2[1:-1] < d2[:-2]) & (d2[1:-1] < d2[2:]) & (d2[1:-1] < 0)
            apex_candidates = (np.where(mid)[0] + 1).tolist()

        # For each apex, find zero-crossings of d2 (neg→pos) as boundaries
        # Zero crossings where d2 goes from negative to positive
        zero_crossings: list[int] = []
        for i in range(n - 1):
            if d2[i] < 0 and d2[i + 1] >= 0:
                zero_crossings.append(i + 1)

        peaks: list[ChromatographicPeak] = []
        for apex in apex_candidates:
            if intensity[apex] < self.min_apex_intensity:
                continue
            # Find nearest zero crossing to the left
            left = 0
            for zc in reversed(zero_crossings):
                if zc < apex:
                    left = zc
                    break
            # Find nearest zero crossing to the right
            right = n - 1
            for zc in zero_crossings:
                if zc > apex:
                    right = zc
                    break
            if right <= left:
                continue
            # Refine apex to actual intensity maximum within bounds
            real_apex = left + int(np.argmax(intensity[left : right + 1]))
            peaks.append(self._build_peak(rt, intensity, real_apex, left, right))
        return peaks


# ---------------------------------------------------------------------------
# Method 5: MatchedFilterPeakPicker
# ---------------------------------------------------------------------------

class MatchedFilterPeakPicker(BasePeakPicker):
    """Cross-correlates the EIC with an ideal peak template.

    Parameters:
        expected_peak_width: Expected peak width in scans.
        correlation_threshold: Minimum normalised cross-correlation score.
        template_type: Shape of the ideal template ('gaussian').
        min_distance: Minimum distance between detected peaks in scans.
    """

    def __init__(
        self,
        expected_peak_width: int = 15,
        correlation_threshold: float = 0.5,
        template_type: str = "gaussian",
        min_distance: int = 5,
    ):
        self.expected_peak_width = max(expected_peak_width, 3)
        self.correlation_threshold = correlation_threshold
        self.template_type = template_type
        self.min_distance = max(min_distance, 1)

    def _make_template(self) -> np.ndarray:
        """Generate the ideal peak template."""
        n = self.expected_peak_width * 3
        if n % 2 == 0:
            n += 1
        x = np.arange(n) - n // 2
        sigma = self.expected_peak_width / 2.355  # FWHM → sigma
        template = np.exp(-0.5 * (x / sigma) ** 2)
        template /= np.linalg.norm(template)
        return template

    # -- public API ---------------------------------------------------------

    def pick_peaks(
        self,
        rt: np.ndarray,
        intensity: np.ndarray,
        **kwargs,
    ) -> list[ChromatographicPeak]:
        if not self._validate_inputs(rt, intensity):
            return []
        rt = np.asarray(rt, dtype=np.float64)
        intensity = np.asarray(intensity, dtype=np.float64)

        template = self._make_template()
        half_tpl = len(template) // 2

        # Normalised cross-correlation via scipy.signal.correlate
        norm_int = intensity / max(np.linalg.norm(intensity), 1e-12)
        cc = signal.correlate(norm_int, template, mode="same")

        # Find peaks in the cross-correlation signal
        peak_indices, properties = signal.find_peaks(
            cc,
            height=self.correlation_threshold,
            distance=self.min_distance,
        )

        peaks: list[ChromatographicPeak] = []
        for pi in peak_indices:
            left = max(0, pi - half_tpl)
            right = min(len(intensity) - 1, pi + half_tpl)
            apex = left + int(np.argmax(intensity[left : right + 1]))
            peaks.append(self._build_peak(rt, intensity, apex, left, right))
        return peaks


# ---------------------------------------------------------------------------
# Adapter to maintain backward compatibility with legacy getPeaksFor() API
# ---------------------------------------------------------------------------

class PeakPickerAdapter:
    """Wraps any BasePeakPicker to expose the legacy ``getPeaksFor()`` interface.

    This allows drop-in replacement in existing MetExtract code that expects
    Bunch objects with attributes: peakIndex, peakScale, peakSNR, peakArea,
    peakLeftFlank, peakRightFlank.
    """

    def __init__(self, picker: BasePeakPicker):
        self.picker = picker

    def getPeaksFor(
        self,
        times: np.ndarray,
        eic: np.ndarray,
        scales: Optional[Sequence[float]] = None,
        snrTh: float = 0.1,
        startIndex: Optional[int] = None,
        endIndex: Optional[int] = None,
    ) -> list:
        """Legacy-compatible peak detection interface.

        Returns list of Bunch-like objects with peakIndex, peakScale, peakSNR,
        peakArea, peakLeftFlank, peakRightFlank.
        """
        from ..utils import Bunch

        times = np.asarray(times, dtype=np.float64)
        eic = np.asarray(eic, dtype=np.float64)

        if startIndex is not None or endIndex is not None:
            si = startIndex if startIndex is not None else 0
            ei = endIndex if endIndex is not None else len(eic) - 1
            sub_times = times[si : ei + 1]
            sub_eic = eic[si : ei + 1]
            offset = si
        else:
            sub_times = times
            sub_eic = eic
            offset = 0

        detected = self.picker.pick_peaks(sub_times, sub_eic)

        results = []
        for pk in detected:
            idx = pk.apex_index + offset
            left_flank = pk.apex_index - pk.start_index
            right_flank = pk.end_index - pk.apex_index
            results.append(
                Bunch(
                    peakIndex=idx,
                    peakScale=(left_flank + right_flank) / 2.0,
                    peakSNR=pk.snr,
                    peakArea=pk.area,
                    peakLeftFlank=left_flank,
                    peakRightFlank=right_flank,
                )
            )
        return results
