"""
Peak Picking Settings Dialog for MetExtract II.

Provides a tabbed dialog where each tab corresponds to one of the five
peak picking algorithms. The active (selected) tab determines which
algorithm will be used for chromatographic peak picking.
"""

from PySide6 import QtCore, QtGui, QtWidgets


class PeakPickingSettingsDialog(QtWidgets.QDialog):
    """Dialog for configuring chromatographic peak picking parameters.

    Each algorithm has its own tab. The currently selected tab is the
    algorithm that will be used for processing.
    """

    # Signal emitted when the user accepts the dialog
    settingsAccepted = QtCore.Signal(dict)

    def __init__(self, parent=None, current_settings: dict | None = None):
        super().__init__(parent)
        self.setWindowTitle("Chromatographic Peak Picking Settings")
        self.setMinimumSize(520, 480)

        if current_settings is None:
            current_settings = {}

        self._settings = dict(current_settings)

        self._build_ui()
        self._load_settings(current_settings)

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        layout = QtWidgets.QVBoxLayout(self)

        # Algorithm tab widget
        self.tabWidget = QtWidgets.QTabWidget()
        layout.addWidget(self.tabWidget)

        # Tab 1: Wavelet Transform (default)
        self._build_wavelet_tab()
        # Tab 2: Gradient Descent
        self._build_gradient_tab()
        # Tab 3: GMM (Curve Fitting)
        self._build_gmm_tab()
        # Tab 4: Second Derivative
        self._build_second_derivative_tab()
        # Tab 5: Matched Filter
        self._build_matched_filter_tab()

        # Buttons
        btn_layout = QtWidgets.QHBoxLayout()
        self.btnOk = QtWidgets.QPushButton("OK")
        self.btnCancel = QtWidgets.QPushButton("Cancel")
        self.btnDefaults = QtWidgets.QPushButton("Restore Defaults")
        btn_layout.addWidget(self.btnDefaults)
        btn_layout.addStretch()
        btn_layout.addWidget(self.btnOk)
        btn_layout.addWidget(self.btnCancel)
        layout.addLayout(btn_layout)

        self.btnOk.clicked.connect(self._accept)
        self.btnCancel.clicked.connect(self.reject)
        self.btnDefaults.clicked.connect(self._restore_defaults)

    # -- Tab builders ---------------------------------------------------

    def _build_wavelet_tab(self):
        tab = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(tab)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)

        self.wt_min_scale = QtWidgets.QSpinBox()
        self.wt_min_scale.setRange(1, 500)
        self.wt_min_scale.setValue(1)
        form.addRow("Min wavelet scale (scans):", self.wt_min_scale)

        self.wt_max_scale = QtWidgets.QSpinBox()
        self.wt_max_scale.setRange(2, 1000)
        self.wt_max_scale.setValue(64)
        form.addRow("Max wavelet scale (scans):", self.wt_max_scale)

        self.wt_num_scales = QtWidgets.QSpinBox()
        self.wt_num_scales.setRange(4, 256)
        self.wt_num_scales.setValue(32)
        form.addRow("Number of scales:", self.wt_num_scales)

        self.wt_snr_threshold = QtWidgets.QDoubleSpinBox()
        self.wt_snr_threshold.setRange(0.0, 1000.0)
        self.wt_snr_threshold.setValue(1.0)
        self.wt_snr_threshold.setDecimals(2)
        form.addRow("SNR threshold:", self.wt_snr_threshold)

        self.wt_min_ridge_length = QtWidgets.QSpinBox()
        self.wt_min_ridge_length.setRange(1, 100)
        self.wt_min_ridge_length.setValue(4)
        form.addRow("Min ridge length:", self.wt_min_ridge_length)

        self.wt_gap_threshold = QtWidgets.QSpinBox()
        self.wt_gap_threshold.setRange(0, 50)
        self.wt_gap_threshold.setValue(2)
        form.addRow("Gap threshold:", self.wt_gap_threshold)

        self.tabWidget.addTab(tab, "Wavelet Transform")

    def _build_gradient_tab(self):
        tab = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(tab)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)

        self.gd_smoothing_method = QtWidgets.QComboBox()
        self.gd_smoothing_method.addItems(["gaussian", "moving_average", "savitzky_golay"])
        form.addRow("Smoothing method:", self.gd_smoothing_method)

        self.gd_smoothing_window = QtWidgets.QSpinBox()
        self.gd_smoothing_window.setRange(1, 201)
        self.gd_smoothing_window.setValue(5)
        form.addRow("Smoothing window size:", self.gd_smoothing_window)

        self.gd_smoothing_iterations = QtWidgets.QSpinBox()
        self.gd_smoothing_iterations.setRange(1, 20)
        self.gd_smoothing_iterations.setValue(1)
        form.addRow("Smoothing iterations:", self.gd_smoothing_iterations)

        self.gd_min_intensity = QtWidgets.QDoubleSpinBox()
        self.gd_min_intensity.setRange(0, 1e12)
        self.gd_min_intensity.setValue(1000.0)
        self.gd_min_intensity.setDecimals(1)
        form.addRow("Min apex intensity:", self.gd_min_intensity)

        self.gd_min_flank_intensity = QtWidgets.QDoubleSpinBox()
        self.gd_min_flank_intensity.setRange(0, 1e12)
        self.gd_min_flank_intensity.setValue(10.0)
        self.gd_min_flank_intensity.setDecimals(1)
        form.addRow("Min flank intensity:", self.gd_min_flank_intensity)

        self.gd_min_increase_ratio = QtWidgets.QDoubleSpinBox()
        self.gd_min_increase_ratio.setRange(0.0, 1.0)
        self.gd_min_increase_ratio.setValue(0.05)
        self.gd_min_increase_ratio.setDecimals(3)
        self.gd_min_increase_ratio.setSingleStep(0.01)
        form.addRow("Min increase ratio:", self.gd_min_increase_ratio)

        self.gd_consecutive_scans = QtWidgets.QSpinBox()
        self.gd_consecutive_scans.setRange(1, 50)
        self.gd_consecutive_scans.setValue(3)
        form.addRow("Consecutive scans (stall check):", self.gd_consecutive_scans)

        self.gd_min_width = QtWidgets.QSpinBox()
        self.gd_min_width.setRange(1, 500)
        self.gd_min_width.setValue(3)
        form.addRow("Min peak width (scans):", self.gd_min_width)

        self.gd_max_width = QtWidgets.QSpinBox()
        self.gd_max_width.setRange(2, 5000)
        self.gd_max_width.setValue(300)
        form.addRow("Max peak width (scans):", self.gd_max_width)

        self.tabWidget.addTab(tab, "Gradient Descent")

    def _build_gmm_tab(self):
        tab = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(tab)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)

        self.gmm_distribution = QtWidgets.QComboBox()
        self.gmm_distribution.addItems(["gaussian", "emg"])
        form.addRow("Distribution type:", self.gmm_distribution)

        self.gmm_baseline_fraction = QtWidgets.QDoubleSpinBox()
        self.gmm_baseline_fraction.setRange(0.0, 10.0)
        self.gmm_baseline_fraction.setValue(0.1)
        self.gmm_baseline_fraction.setDecimals(3)
        self.gmm_baseline_fraction.setSingleStep(0.01)
        form.addRow("Baseline fraction:", self.gmm_baseline_fraction)

        self.gmm_min_chunk_size = QtWidgets.QSpinBox()
        self.gmm_min_chunk_size.setRange(2, 500)
        self.gmm_min_chunk_size.setValue(5)
        form.addRow("Min chunk size (scans):", self.gmm_min_chunk_size)

        self.gmm_max_peaks_per_chunk = QtWidgets.QSpinBox()
        self.gmm_max_peaks_per_chunk.setRange(1, 50)
        self.gmm_max_peaks_per_chunk.setValue(10)
        form.addRow("Max peaks per chunk:", self.gmm_max_peaks_per_chunk)

        self.gmm_fit_max_iterations = QtWidgets.QSpinBox()
        self.gmm_fit_max_iterations.setRange(100, 100000)
        self.gmm_fit_max_iterations.setValue(5000)
        self.gmm_fit_max_iterations.setSingleStep(500)
        form.addRow("Fit max iterations:", self.gmm_fit_max_iterations)

        self.tabWidget.addTab(tab, "GMM (Curve Fitting)")

    def _build_second_derivative_tab(self):
        tab = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(tab)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)

        self.sd_sg_window = QtWidgets.QSpinBox()
        self.sd_sg_window.setRange(5, 201)
        self.sd_sg_window.setValue(11)
        self.sd_sg_window.setSingleStep(2)
        form.addRow("Savitzky-Golay window:", self.sd_sg_window)

        self.sd_sg_polyorder = QtWidgets.QSpinBox()
        self.sd_sg_polyorder.setRange(2, 10)
        self.sd_sg_polyorder.setValue(3)
        form.addRow("SG polynomial order:", self.sd_sg_polyorder)

        self.sd_pre_smoothing = QtWidgets.QComboBox()
        self.sd_pre_smoothing.addItems(["none", "gaussian", "moving_average", "savitzky_golay"])
        form.addRow("Pre-smoothing method:", self.sd_pre_smoothing)

        self.sd_pre_smoothing_window = QtWidgets.QSpinBox()
        self.sd_pre_smoothing_window.setRange(1, 201)
        self.sd_pre_smoothing_window.setValue(5)
        form.addRow("Pre-smoothing window:", self.sd_pre_smoothing_window)

        self.sd_pre_smoothing_iterations = QtWidgets.QSpinBox()
        self.sd_pre_smoothing_iterations.setRange(1, 20)
        self.sd_pre_smoothing_iterations.setValue(1)
        form.addRow("Pre-smoothing iterations:", self.sd_pre_smoothing_iterations)

        self.sd_min_apex_intensity = QtWidgets.QDoubleSpinBox()
        self.sd_min_apex_intensity.setRange(0, 1e12)
        self.sd_min_apex_intensity.setValue(0.0)
        self.sd_min_apex_intensity.setDecimals(1)
        form.addRow("Min apex intensity:", self.sd_min_apex_intensity)

        self.tabWidget.addTab(tab, "Second Derivative")

    def _build_matched_filter_tab(self):
        tab = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(tab)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)

        self.mf_expected_peak_width = QtWidgets.QSpinBox()
        self.mf_expected_peak_width.setRange(3, 500)
        self.mf_expected_peak_width.setValue(15)
        form.addRow("Expected peak width (scans):", self.mf_expected_peak_width)

        self.mf_correlation_threshold = QtWidgets.QDoubleSpinBox()
        self.mf_correlation_threshold.setRange(0.0, 1.0)
        self.mf_correlation_threshold.setValue(0.5)
        self.mf_correlation_threshold.setDecimals(3)
        self.mf_correlation_threshold.setSingleStep(0.05)
        form.addRow("Correlation threshold:", self.mf_correlation_threshold)

        self.mf_template_type = QtWidgets.QComboBox()
        self.mf_template_type.addItems(["gaussian"])
        form.addRow("Template type:", self.mf_template_type)

        self.mf_min_distance = QtWidgets.QSpinBox()
        self.mf_min_distance.setRange(1, 500)
        self.mf_min_distance.setValue(5)
        form.addRow("Min distance between peaks (scans):", self.mf_min_distance)

        self.tabWidget.addTab(tab, "Matched Filter")

    # ------------------------------------------------------------------
    # Settings management
    # ------------------------------------------------------------------

    _ALGORITHM_TAB_MAP = {
        "wavelettransform": 0,
        "gradientdescent": 1,
        "gmm": 2,
        "secondderivative": 3,
        "matchedfilter": 4,
    }

    _TAB_ALGORITHM_MAP = {v: k for k, v in _ALGORITHM_TAB_MAP.items()}

    def _load_settings(self, s: dict):
        """Populate controls from a settings dictionary."""
        algo = s.get("algorithm", "wavelettransform").lower()
        tab_idx = self._ALGORITHM_TAB_MAP.get(algo, 0)
        self.tabWidget.setCurrentIndex(tab_idx)

        # Wavelet Transform
        self.wt_min_scale.setValue(s.get("wt_min_scale", 1))
        self.wt_max_scale.setValue(s.get("wt_max_scale", 64))
        self.wt_num_scales.setValue(s.get("wt_num_scales", 32))
        self.wt_snr_threshold.setValue(s.get("wt_snr_threshold", 1.0))
        self.wt_min_ridge_length.setValue(s.get("wt_min_ridge_length", 4))
        self.wt_gap_threshold.setValue(s.get("wt_gap_threshold", 2))

        # Gradient Descent
        sm = s.get("gd_smoothing_method", "gaussian")
        idx = self.gd_smoothing_method.findText(sm)
        if idx >= 0:
            self.gd_smoothing_method.setCurrentIndex(idx)
        self.gd_smoothing_window.setValue(s.get("gd_smoothing_window", 5))
        self.gd_smoothing_iterations.setValue(s.get("gd_smoothing_iterations", 1))
        self.gd_min_intensity.setValue(s.get("gd_min_intensity", 1000.0))
        self.gd_min_flank_intensity.setValue(s.get("gd_min_flank_intensity", 10.0))
        self.gd_min_increase_ratio.setValue(s.get("gd_min_increase_ratio", 0.05))
        self.gd_consecutive_scans.setValue(s.get("gd_consecutive_scans", 3))
        self.gd_min_width.setValue(s.get("gd_min_width", 3))
        self.gd_max_width.setValue(s.get("gd_max_width", 300))

        # GMM
        dist = s.get("gmm_distribution", "gaussian")
        idx = self.gmm_distribution.findText(dist)
        if idx >= 0:
            self.gmm_distribution.setCurrentIndex(idx)
        self.gmm_baseline_fraction.setValue(s.get("gmm_baseline_fraction", 0.1))
        self.gmm_min_chunk_size.setValue(s.get("gmm_min_chunk_size", 5))
        self.gmm_max_peaks_per_chunk.setValue(s.get("gmm_max_peaks_per_chunk", 10))
        self.gmm_fit_max_iterations.setValue(s.get("gmm_fit_max_iterations", 5000))

        # Second Derivative
        self.sd_sg_window.setValue(s.get("sd_sg_window", 11))
        self.sd_sg_polyorder.setValue(s.get("sd_sg_polyorder", 3))
        ps = s.get("sd_pre_smoothing", "none")
        idx = self.sd_pre_smoothing.findText(ps)
        if idx >= 0:
            self.sd_pre_smoothing.setCurrentIndex(idx)
        self.sd_pre_smoothing_window.setValue(s.get("sd_pre_smoothing_window", 5))
        self.sd_pre_smoothing_iterations.setValue(s.get("sd_pre_smoothing_iterations", 1))
        self.sd_min_apex_intensity.setValue(s.get("sd_min_apex_intensity", 0.0))

        # Matched Filter
        self.mf_expected_peak_width.setValue(s.get("mf_expected_peak_width", 15))
        self.mf_correlation_threshold.setValue(s.get("mf_correlation_threshold", 0.5))
        tt = s.get("mf_template_type", "gaussian")
        idx = self.mf_template_type.findText(tt)
        if idx >= 0:
            self.mf_template_type.setCurrentIndex(idx)
        self.mf_min_distance.setValue(s.get("mf_min_distance", 5))

    def _collect_settings(self) -> dict:
        """Read all controls and return a flat settings dict."""
        tab_idx = self.tabWidget.currentIndex()
        algo = self._TAB_ALGORITHM_MAP.get(tab_idx, "wavelettransform")
        s = {"algorithm": algo}

        # Wavelet Transform
        s["wt_min_scale"] = self.wt_min_scale.value()
        s["wt_max_scale"] = self.wt_max_scale.value()
        s["wt_num_scales"] = self.wt_num_scales.value()
        s["wt_snr_threshold"] = self.wt_snr_threshold.value()
        s["wt_min_ridge_length"] = self.wt_min_ridge_length.value()
        s["wt_gap_threshold"] = self.wt_gap_threshold.value()

        # Gradient Descent
        s["gd_smoothing_method"] = self.gd_smoothing_method.currentText()
        s["gd_smoothing_window"] = self.gd_smoothing_window.value()
        s["gd_smoothing_iterations"] = self.gd_smoothing_iterations.value()
        s["gd_min_intensity"] = self.gd_min_intensity.value()
        s["gd_min_flank_intensity"] = self.gd_min_flank_intensity.value()
        s["gd_min_increase_ratio"] = self.gd_min_increase_ratio.value()
        s["gd_consecutive_scans"] = self.gd_consecutive_scans.value()
        s["gd_min_width"] = self.gd_min_width.value()
        s["gd_max_width"] = self.gd_max_width.value()

        # GMM
        s["gmm_distribution"] = self.gmm_distribution.currentText()
        s["gmm_baseline_fraction"] = self.gmm_baseline_fraction.value()
        s["gmm_min_chunk_size"] = self.gmm_min_chunk_size.value()
        s["gmm_max_peaks_per_chunk"] = self.gmm_max_peaks_per_chunk.value()
        s["gmm_fit_max_iterations"] = self.gmm_fit_max_iterations.value()

        # Second Derivative
        s["sd_sg_window"] = self.sd_sg_window.value()
        s["sd_sg_polyorder"] = self.sd_sg_polyorder.value()
        s["sd_pre_smoothing"] = self.sd_pre_smoothing.currentText()
        s["sd_pre_smoothing_window"] = self.sd_pre_smoothing_window.value()
        s["sd_pre_smoothing_iterations"] = self.sd_pre_smoothing_iterations.value()
        s["sd_min_apex_intensity"] = self.sd_min_apex_intensity.value()

        # Matched Filter
        s["mf_expected_peak_width"] = self.mf_expected_peak_width.value()
        s["mf_correlation_threshold"] = self.mf_correlation_threshold.value()
        s["mf_template_type"] = self.mf_template_type.currentText()
        s["mf_min_distance"] = self.mf_min_distance.value()

        return s

    def _accept(self):
        self._settings = self._collect_settings()
        self.settingsAccepted.emit(self._settings)
        self.accept()

    def _restore_defaults(self):
        self._load_settings({})

    def get_settings(self) -> dict:
        """Return the last accepted settings dict."""
        return dict(self._settings)

    @staticmethod
    def get_default_settings() -> dict:
        """Return default settings."""
        return {"algorithm": "wavelettransform"}


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    dlg = PeakPickingSettingsDialog()
    if dlg.exec() == QtWidgets.QDialog.Accepted:
        print("Settings:", dlg.get_settings())
    sys.exit(0)
