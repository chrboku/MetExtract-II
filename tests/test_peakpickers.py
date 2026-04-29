"""Tests for the peak picking architecture in src/chromPeakPicking/peakpickers.py."""

import numpy as np
import pytest

from src.chromPeakPicking.peakpickers import (
    BasePeakPicker,
    ChromatographicPeak,
    DistributionType,
    GMMPeakPicker,
    GradientDescentPeakPicker,
    MatchedFilterPeakPicker,
    PeakPickerAdapter,
    SecondDerivativePeakPicker,
    SmoothingMethod,
    WaveletTransformPeakPicker,
)

# ---------------------------------------------------------------------------
# Fixtures – synthetic EICs
# ---------------------------------------------------------------------------


@pytest.fixture
def simple_gaussian_eic():
    """Single Gaussian peak at RT=300 s, sigma=10 s."""
    rt = np.linspace(0, 600, 601)
    intensity = 10000 * np.exp(-0.5 * ((rt - 300) / 10) ** 2)
    return rt, intensity


@pytest.fixture
def two_peak_eic():
    """Two well-separated Gaussian peaks at RT=200 and RT=400."""
    rt = np.linspace(0, 600, 601)
    peak1 = 8000 * np.exp(-0.5 * ((rt - 200) / 8) ** 2)
    peak2 = 12000 * np.exp(-0.5 * ((rt - 400) / 12) ** 2)
    intensity = peak1 + peak2
    return rt, intensity


@pytest.fixture
def noisy_eic():
    """Gaussian peak with added noise."""
    rng = np.random.default_rng(42)
    rt = np.linspace(0, 600, 601)
    intensity = 5000 * np.exp(-0.5 * ((rt - 300) / 15) ** 2) + rng.normal(0, 100, 601)
    intensity = np.maximum(intensity, 0)
    return rt, intensity


@pytest.fixture
def flat_eic():
    """Completely flat signal."""
    rt = np.linspace(0, 600, 601)
    intensity = np.full_like(rt, 500.0)
    return rt, intensity


@pytest.fixture
def empty_eic():
    """Empty arrays."""
    return np.array([]), np.array([])


# ---------------------------------------------------------------------------
# ChromatographicPeak dataclass
# ---------------------------------------------------------------------------


class TestChromatographicPeak:
    def test_creation(self):
        pk = ChromatographicPeak(
            start_rt=10.0,
            end_rt=20.0,
            apex_rt=15.0,
            fwhm=5.0,
            snr=10.0,
            area=1234.0,
        )
        assert pk.apex_rt == 15.0
        assert pk.snr == 10.0

    def test_defaults(self):
        pk = ChromatographicPeak(
            start_rt=0,
            end_rt=1,
            apex_rt=0.5,
            fwhm=0.5,
            snr=1.0,
            area=0.5,
        )
        assert pk.apex_intensity == 0.0
        assert pk.baseline == 0.0


# ---------------------------------------------------------------------------
# BasePeakPicker utility methods
# ---------------------------------------------------------------------------


class TestBasePeakPickerUtilities:
    def test_validate_inputs_empty(self, empty_eic):
        assert BasePeakPicker._validate_inputs(*empty_eic) is False

    def test_validate_inputs_flat(self, flat_eic):
        assert BasePeakPicker._validate_inputs(*flat_eic) is False

    def test_validate_inputs_valid(self, simple_gaussian_eic):
        assert BasePeakPicker._validate_inputs(*simple_gaussian_eic) is True

    def test_validate_inputs_mismatched(self):
        rt = np.arange(10)
        intensity = np.arange(5)
        assert BasePeakPicker._validate_inputs(rt, intensity) is False

    def test_validate_inputs_none(self):
        assert BasePeakPicker._validate_inputs(None, None) is False

    def test_compute_area(self, simple_gaussian_eic):
        rt, intensity = simple_gaussian_eic
        area = BasePeakPicker.compute_area(rt, intensity, 250, 350)
        assert area > 0

    def test_compute_snr(self, simple_gaussian_eic):
        rt, intensity = simple_gaussian_eic
        apex = int(np.argmax(intensity))
        snr = BasePeakPicker.compute_snr(intensity, apex, 250, 350)
        assert snr > 0

    def test_compute_fwhm(self, simple_gaussian_eic):
        rt, intensity = simple_gaussian_eic
        apex = int(np.argmax(intensity))
        fwhm = BasePeakPicker.compute_fwhm(rt, intensity, apex, 250, 350)
        # FWHM ≈ 2.355 * sigma = 23.55 for sigma=10
        assert 15 < fwhm < 35

    def test_smooth_eic_gaussian(self, simple_gaussian_eic):
        _, intensity = simple_gaussian_eic
        smoothed = BasePeakPicker.smooth_eic(intensity, method=SmoothingMethod.GAUSSIAN)
        assert smoothed.shape == intensity.shape

    def test_smooth_eic_moving_average(self, simple_gaussian_eic):
        _, intensity = simple_gaussian_eic
        smoothed = BasePeakPicker.smooth_eic(intensity, method=SmoothingMethod.MOVING_AVERAGE)
        assert smoothed.shape == intensity.shape

    def test_smooth_eic_savgol(self, simple_gaussian_eic):
        _, intensity = simple_gaussian_eic
        smoothed = BasePeakPicker.smooth_eic(intensity, method=SmoothingMethod.SAVITZKY_GOLAY)
        assert smoothed.shape == intensity.shape

    def test_smooth_eic_iterations(self, simple_gaussian_eic):
        _, intensity = simple_gaussian_eic
        s1 = BasePeakPicker.smooth_eic(intensity, iterations=1)
        s3 = BasePeakPicker.smooth_eic(intensity, iterations=3)
        # More iterations → smoother → lower max (for a sharp peak)
        assert np.max(s3) <= np.max(s1) + 1e-6


# ---------------------------------------------------------------------------
# GradientDescentPeakPicker
# ---------------------------------------------------------------------------


class TestGradientDescentPeakPicker:
    def test_single_peak(self, simple_gaussian_eic):
        picker = GradientDescentPeakPicker(min_intensity=100)
        peaks = picker.pick_peaks(*simple_gaussian_eic)
        assert len(peaks) >= 1
        apex_rts = [p.apex_rt for p in peaks]
        assert any(280 < rt < 320 for rt in apex_rts)

    def test_two_peaks(self, two_peak_eic):
        picker = GradientDescentPeakPicker(min_intensity=100)
        peaks = picker.pick_peaks(*two_peak_eic)
        assert len(peaks) >= 2

    def test_empty_input(self, empty_eic):
        picker = GradientDescentPeakPicker()
        assert picker.pick_peaks(*empty_eic) == []

    def test_flat_input(self, flat_eic):
        picker = GradientDescentPeakPicker()
        assert picker.pick_peaks(*flat_eic) == []

    def test_width_filter(self, simple_gaussian_eic):
        picker = GradientDescentPeakPicker(min_intensity=100, min_width=500)
        peaks = picker.pick_peaks(*simple_gaussian_eic)
        assert len(peaks) == 0  # peak too narrow

    def test_peak_has_valid_metrics(self, simple_gaussian_eic):
        picker = GradientDescentPeakPicker(min_intensity=100)
        peaks = picker.pick_peaks(*simple_gaussian_eic)
        for p in peaks:
            assert p.area > 0
            assert p.fwhm >= 0
            assert p.start_rt < p.apex_rt < p.end_rt


# ---------------------------------------------------------------------------
# GMMPeakPicker
# ---------------------------------------------------------------------------


class TestGMMPeakPicker:
    def test_single_peak(self, simple_gaussian_eic):
        picker = GMMPeakPicker(baseline_fraction=0.01)
        peaks = picker.pick_peaks(*simple_gaussian_eic)
        assert len(peaks) >= 1
        apex_rts = [p.apex_rt for p in peaks]
        assert any(280 < rt < 320 for rt in apex_rts)

    def test_two_peaks(self, two_peak_eic):
        picker = GMMPeakPicker(baseline_fraction=0.01)
        peaks = picker.pick_peaks(*two_peak_eic)
        assert len(peaks) >= 2

    def test_empty_input(self, empty_eic):
        picker = GMMPeakPicker()
        assert picker.pick_peaks(*empty_eic) == []

    def test_flat_input(self, flat_eic):
        picker = GMMPeakPicker()
        assert picker.pick_peaks(*flat_eic) == []

    def test_emg_distribution(self, simple_gaussian_eic):
        picker = GMMPeakPicker(baseline_fraction=0.01, distribution=DistributionType.EMG)
        peaks = picker.pick_peaks(*simple_gaussian_eic)
        assert len(peaks) >= 1


# ---------------------------------------------------------------------------
# WaveletTransformPeakPicker
# ---------------------------------------------------------------------------


class TestWaveletTransformPeakPicker:
    def test_single_peak(self, simple_gaussian_eic):
        picker = WaveletTransformPeakPicker(min_scale=1, max_scale=30, snr_threshold=0.5, min_ridge_length=2)
        peaks = picker.pick_peaks(*simple_gaussian_eic)
        assert len(peaks) >= 1
        apex_rts = [p.apex_rt for p in peaks]
        assert any(280 < rt < 320 for rt in apex_rts)

    def test_two_peaks(self, two_peak_eic):
        picker = WaveletTransformPeakPicker(min_scale=1, max_scale=30, snr_threshold=0.5, min_ridge_length=2)
        peaks = picker.pick_peaks(*two_peak_eic)
        assert len(peaks) >= 2

    def test_empty_input(self, empty_eic):
        picker = WaveletTransformPeakPicker()
        assert picker.pick_peaks(*empty_eic) == []

    def test_flat_input(self, flat_eic):
        picker = WaveletTransformPeakPicker()
        assert picker.pick_peaks(*flat_eic) == []


# ---------------------------------------------------------------------------
# SecondDerivativePeakPicker
# ---------------------------------------------------------------------------


class TestSecondDerivativePeakPicker:
    def test_single_peak(self, simple_gaussian_eic):
        picker = SecondDerivativePeakPicker(sg_window=11, sg_polyorder=3)
        peaks = picker.pick_peaks(*simple_gaussian_eic)
        assert len(peaks) >= 1
        apex_rts = [p.apex_rt for p in peaks]
        assert any(280 < rt < 320 for rt in apex_rts)

    def test_two_peaks(self, two_peak_eic):
        picker = SecondDerivativePeakPicker(sg_window=11, sg_polyorder=3)
        peaks = picker.pick_peaks(*two_peak_eic)
        assert len(peaks) >= 2

    def test_empty_input(self, empty_eic):
        picker = SecondDerivativePeakPicker()
        assert picker.pick_peaks(*empty_eic) == []

    def test_flat_input(self, flat_eic):
        picker = SecondDerivativePeakPicker()
        assert picker.pick_peaks(*flat_eic) == []

    def test_with_pre_smoothing(self, noisy_eic):
        picker = SecondDerivativePeakPicker(
            sg_window=15,
            sg_polyorder=3,
            pre_smoothing=SmoothingMethod.GAUSSIAN,
            pre_smoothing_window=5,
            pre_smoothing_iterations=2,
        )
        peaks = picker.pick_peaks(*noisy_eic)
        assert len(peaks) >= 1


# ---------------------------------------------------------------------------
# MatchedFilterPeakPicker
# ---------------------------------------------------------------------------


class TestMatchedFilterPeakPicker:
    def test_single_peak(self, simple_gaussian_eic):
        picker = MatchedFilterPeakPicker(expected_peak_width=20, correlation_threshold=0.3)
        peaks = picker.pick_peaks(*simple_gaussian_eic)
        assert len(peaks) >= 1
        apex_rts = [p.apex_rt for p in peaks]
        assert any(280 < rt < 320 for rt in apex_rts)

    def test_two_peaks(self, two_peak_eic):
        picker = MatchedFilterPeakPicker(expected_peak_width=20, correlation_threshold=0.3)
        peaks = picker.pick_peaks(*two_peak_eic)
        assert len(peaks) >= 2

    def test_empty_input(self, empty_eic):
        picker = MatchedFilterPeakPicker()
        assert picker.pick_peaks(*empty_eic) == []

    def test_flat_input(self, flat_eic):
        picker = MatchedFilterPeakPicker()
        assert picker.pick_peaks(*flat_eic) == []


# ---------------------------------------------------------------------------
# PeakPickerAdapter (legacy interface)
# ---------------------------------------------------------------------------


class TestPeakPickerAdapter:
    def test_adapter_returns_bunch_objects(self, simple_gaussian_eic):
        picker = WaveletTransformPeakPicker(min_scale=1, max_scale=30, snr_threshold=0.5, min_ridge_length=2)
        adapter = PeakPickerAdapter(picker)
        rt, intensity = simple_gaussian_eic
        results = adapter.getPeaksFor(rt, intensity)
        assert len(results) >= 1
        for r in results:
            assert hasattr(r, "peakIndex")
            assert hasattr(r, "peakScale")
            assert hasattr(r, "peakSNR")
            assert hasattr(r, "peakArea")
            assert hasattr(r, "peakLeftFlank")
            assert hasattr(r, "peakRightFlank")

    def test_adapter_with_index_range(self, simple_gaussian_eic):
        picker = WaveletTransformPeakPicker(min_scale=1, max_scale=30, snr_threshold=0.5, min_ridge_length=2)
        adapter = PeakPickerAdapter(picker)
        rt, intensity = simple_gaussian_eic
        results = adapter.getPeaksFor(rt, intensity, startIndex=200, endIndex=400)
        for r in results:
            assert 200 <= r.peakIndex <= 400


# ---------------------------------------------------------------------------
# MassSpecWavelet drop-in replacement
# ---------------------------------------------------------------------------


class TestMassSpecWaveletDropIn:
    def test_legacy_interface(self, simple_gaussian_eic):
        from src.chromPeakPicking.MassSpecWavelet import MassSpecWavelet

        cp = MassSpecWavelet(scales=[2, 30], snrTh=0.5, minScans=1)
        rt, intensity = simple_gaussian_eic
        results = cp.getPeaksFor(rt, intensity)
        assert isinstance(results, list)
        # Should find the peak
        assert len(results) >= 1
        for r in results:
            assert hasattr(r, "peakIndex")


# ---------------------------------------------------------------------------
# Baseline replacement
# ---------------------------------------------------------------------------


class TestBaseline:
    def test_median_baseline(self, simple_gaussian_eic):
        from src.Baseline import Baseline

        bl = Baseline()
        rt, intensity = simple_gaussian_eic
        baseline = bl.getBaseline(intensity, rt)
        assert baseline.shape == intensity.shape
        # Baseline should be below the peak apex
        apex_idx = np.argmax(intensity)
        assert baseline[apex_idx] < intensity[apex_idx]

    def test_empty_eic(self):
        from src.Baseline import Baseline

        bl = Baseline()
        result = bl.getBaseline(np.array([]), np.array([]))
        assert len(result) == 0
