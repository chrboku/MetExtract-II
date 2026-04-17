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
        self.wt_min_scale.setValue(3)
        self.wt_min_scale.setToolTip("Minimum wavelet scale in scans. Defines the narrowest peak width the algorithm will detect. Smaller values find narrower peaks.")
        lbl = QtWidgets.QLabel("Min wavelet scale (scans):")
        lbl.setToolTip(self.wt_min_scale.toolTip())
        form.addRow(lbl, self.wt_min_scale)

        self.wt_max_scale = QtWidgets.QSpinBox()
        self.wt_max_scale.setRange(2, 1000)
        self.wt_max_scale.setValue(11)
        self.wt_max_scale.setToolTip("Maximum wavelet scale in scans. Defines the widest peak width the algorithm will detect. Larger values find broader peaks.")
        lbl = QtWidgets.QLabel("Max wavelet scale (scans):")
        lbl.setToolTip(self.wt_max_scale.toolTip())
        form.addRow(lbl, self.wt_max_scale)

        self.wt_num_scales = QtWidgets.QSpinBox()
        self.wt_num_scales.setRange(4, 256)
        self.wt_num_scales.setValue(8)
        self.wt_num_scales.setToolTip("Number of wavelet scales evaluated between min and max. More scales give finer resolution but increase computation time.")
        lbl = QtWidgets.QLabel("Number of scales:")
        lbl.setToolTip(self.wt_num_scales.toolTip())
        form.addRow(lbl, self.wt_num_scales)

        self.wt_snr_threshold = QtWidgets.QDoubleSpinBox()
        self.wt_snr_threshold.setRange(0.0, 1000.0)
        self.wt_snr_threshold.setValue(3.0)
        self.wt_snr_threshold.setDecimals(2)
        self.wt_snr_threshold.setToolTip("Signal-to-noise ratio threshold for accepting a peak. Higher values reject weaker peaks and reduce false positives.")
        lbl = QtWidgets.QLabel("SNR threshold:")
        lbl.setToolTip(self.wt_snr_threshold.toolTip())
        form.addRow(lbl, self.wt_snr_threshold)

        self.wt_min_ridge_length = QtWidgets.QSpinBox()
        self.wt_min_ridge_length.setRange(1, 100)
        self.wt_min_ridge_length.setValue(4)
        self.wt_min_ridge_length.setToolTip("Minimum number of consecutive wavelet scales a ridge must span to be considered a real peak. Filters out noise ridges.")
        lbl = QtWidgets.QLabel("Min ridge length:")
        lbl.setToolTip(self.wt_min_ridge_length.toolTip())
        form.addRow(lbl, self.wt_min_ridge_length)

        self.wt_gap_threshold = QtWidgets.QSpinBox()
        self.wt_gap_threshold.setRange(0, 50)
        self.wt_gap_threshold.setValue(2)
        self.wt_gap_threshold.setToolTip("Maximum number of consecutive scale gaps allowed in a ridge before it is terminated. Small values enforce stricter continuity.")
        lbl = QtWidgets.QLabel("Gap threshold:")
        lbl.setToolTip(self.wt_gap_threshold.toolTip())
        form.addRow(lbl, self.wt_gap_threshold)

        self.tabWidget.addTab(tab, "Wavelet Transform")

    def _build_gradient_tab(self):
        tab = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(tab)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)

        self.gd_smoothing_method = QtWidgets.QComboBox()
        self.gd_smoothing_method.addItems(["gaussian", "moving_average", "savitzky_golay"])
        self.gd_smoothing_method.setToolTip("Smoothing filter applied to the EIC before peak detection. Gaussian is recommended for most LC-HRMS data.")
        lbl = QtWidgets.QLabel("Smoothing method:")
        lbl.setToolTip(self.gd_smoothing_method.toolTip())
        form.addRow(lbl, self.gd_smoothing_method)

        self.gd_smoothing_window = QtWidgets.QSpinBox()
        self.gd_smoothing_window.setRange(1, 201)
        self.gd_smoothing_window.setValue(5)
        self.gd_smoothing_window.setToolTip("Width of the smoothing window in scans. Larger values produce stronger smoothing but may merge close peaks.")
        lbl = QtWidgets.QLabel("Smoothing window size:")
        lbl.setToolTip(self.gd_smoothing_window.toolTip())
        form.addRow(lbl, self.gd_smoothing_window)

        self.gd_smoothing_iterations = QtWidgets.QSpinBox()
        self.gd_smoothing_iterations.setRange(1, 20)
        self.gd_smoothing_iterations.setValue(1)
        self.gd_smoothing_iterations.setToolTip("Number of times the smoothing filter is applied. Multiple passes give stronger noise reduction.")
        lbl = QtWidgets.QLabel("Smoothing iterations:")
        lbl.setToolTip(self.gd_smoothing_iterations.toolTip())
        form.addRow(lbl, self.gd_smoothing_iterations)

        self.gd_min_intensity = QtWidgets.QDoubleSpinBox()
        self.gd_min_intensity.setRange(0, 1e12)
        self.gd_min_intensity.setValue(1000.0)
        self.gd_min_intensity.setDecimals(1)
        self.gd_min_intensity.setToolTip("Minimum intensity at the peak apex for a peak to be reported. Peaks below this threshold are discarded.")
        lbl = QtWidgets.QLabel("Min apex intensity:")
        lbl.setToolTip(self.gd_min_intensity.toolTip())
        form.addRow(lbl, self.gd_min_intensity)

        self.gd_min_flank_intensity = QtWidgets.QDoubleSpinBox()
        self.gd_min_flank_intensity.setRange(0, 1e12)
        self.gd_min_flank_intensity.setValue(10.0)
        self.gd_min_flank_intensity.setDecimals(1)
        self.gd_min_flank_intensity.setToolTip("Minimum intensity at peak boundaries. The descent from the apex stops when intensity drops below this value.")
        lbl = QtWidgets.QLabel("Min flank intensity:")
        lbl.setToolTip(self.gd_min_flank_intensity.toolTip())
        form.addRow(lbl, self.gd_min_flank_intensity)

        self.gd_min_increase_ratio = QtWidgets.QDoubleSpinBox()
        self.gd_min_increase_ratio.setRange(0.0, 1.0)
        self.gd_min_increase_ratio.setValue(0.05)
        self.gd_min_increase_ratio.setDecimals(3)
        self.gd_min_increase_ratio.setSingleStep(0.01)
        self.gd_min_increase_ratio.setToolTip("Minimum relative intensity increase required per scan during slope descent. Used to detect when the signal stalls.")
        lbl = QtWidgets.QLabel("Min increase ratio:")
        lbl.setToolTip(self.gd_min_increase_ratio.toolTip())
        form.addRow(lbl, self.gd_min_increase_ratio)

        self.gd_consecutive_scans = QtWidgets.QSpinBox()
        self.gd_consecutive_scans.setRange(1, 50)
        self.gd_consecutive_scans.setValue(3)
        self.gd_consecutive_scans.setToolTip("Number of consecutive scans with insufficient increase to trigger a stall. Descent stops when this many flat scans are observed.")
        lbl = QtWidgets.QLabel("Consecutive scans (stall check):")
        lbl.setToolTip(self.gd_consecutive_scans.toolTip())
        form.addRow(lbl, self.gd_consecutive_scans)

        self.gd_min_width = QtWidgets.QSpinBox()
        self.gd_min_width.setRange(1, 500)
        self.gd_min_width.setValue(3)
        self.gd_min_width.setToolTip("Minimum peak width in scans. Peaks narrower than this are discarded as likely noise spikes.")
        lbl = QtWidgets.QLabel("Min peak width (scans):")
        lbl.setToolTip(self.gd_min_width.toolTip())
        form.addRow(lbl, self.gd_min_width)

        self.gd_max_width = QtWidgets.QSpinBox()
        self.gd_max_width.setRange(2, 5000)
        self.gd_max_width.setValue(300)
        self.gd_max_width.setToolTip("Maximum peak width in scans. Peaks wider than this are discarded as likely baseline drift artifacts.")
        lbl = QtWidgets.QLabel("Max peak width (scans):")
        lbl.setToolTip(self.gd_max_width.toolTip())
        form.addRow(lbl, self.gd_max_width)

        self.tabWidget.addTab(tab, "Gradient Descent")

    def _build_gmm_tab(self):
        tab = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(tab)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)

        self.gmm_distribution = QtWidgets.QComboBox()
        self.gmm_distribution.addItems(["gaussian", "emg"])
        self.gmm_distribution.setToolTip("Distribution model used for curve fitting. 'gaussian' fits symmetric peaks; 'emg' (Exponentially Modified Gaussian) fits asymmetric/tailing peaks.")
        lbl = QtWidgets.QLabel("Distribution type:")
        lbl.setToolTip(self.gmm_distribution.toolTip())
        form.addRow(lbl, self.gmm_distribution)

        self.gmm_baseline_fraction = QtWidgets.QDoubleSpinBox()
        self.gmm_baseline_fraction.setRange(0.0, 10.0)
        self.gmm_baseline_fraction.setValue(0.1)
        self.gmm_baseline_fraction.setDecimals(3)
        self.gmm_baseline_fraction.setSingleStep(0.01)
        self.gmm_baseline_fraction.setToolTip("Fraction of the maximum intensity used to define the baseline. Regions below this threshold are used to split the EIC into chunks.")
        lbl = QtWidgets.QLabel("Baseline fraction:")
        lbl.setToolTip(self.gmm_baseline_fraction.toolTip())
        form.addRow(lbl, self.gmm_baseline_fraction)

        self.gmm_min_chunk_size = QtWidgets.QSpinBox()
        self.gmm_min_chunk_size.setRange(2, 500)
        self.gmm_min_chunk_size.setValue(5)
        self.gmm_min_chunk_size.setToolTip("Minimum number of scans in a chunk for it to be processed. Chunks shorter than this are skipped.")
        lbl = QtWidgets.QLabel("Min chunk size (scans):")
        lbl.setToolTip(self.gmm_min_chunk_size.toolTip())
        form.addRow(lbl, self.gmm_min_chunk_size)

        self.gmm_max_peaks_per_chunk = QtWidgets.QSpinBox()
        self.gmm_max_peaks_per_chunk.setRange(1, 50)
        self.gmm_max_peaks_per_chunk.setValue(10)
        self.gmm_max_peaks_per_chunk.setToolTip("Maximum number of overlapping peaks to fit in a single chunk. Higher values allow more complex deconvolutions but increase computation time.")
        lbl = QtWidgets.QLabel("Max peaks per chunk:")
        lbl.setToolTip(self.gmm_max_peaks_per_chunk.toolTip())
        form.addRow(lbl, self.gmm_max_peaks_per_chunk)

        self.gmm_fit_max_iterations = QtWidgets.QSpinBox()
        self.gmm_fit_max_iterations.setRange(100, 100000)
        self.gmm_fit_max_iterations.setValue(5000)
        self.gmm_fit_max_iterations.setSingleStep(500)
        self.gmm_fit_max_iterations.setToolTip("Maximum iterations for the curve fitting optimizer. Increase if convergence warnings appear for complex peaks.")
        lbl = QtWidgets.QLabel("Fit max iterations:")
        lbl.setToolTip(self.gmm_fit_max_iterations.toolTip())
        form.addRow(lbl, self.gmm_fit_max_iterations)

        self.tabWidget.addTab(tab, "GMM (Curve Fitting)")

    def _build_second_derivative_tab(self):
        tab = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(tab)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)

        self.sd_sg_window = QtWidgets.QSpinBox()
        self.sd_sg_window.setRange(5, 201)
        self.sd_sg_window.setValue(11)
        self.sd_sg_window.setSingleStep(2)
        self.sd_sg_window.setToolTip("Window length for the Savitzky-Golay second derivative filter. Must be odd. Larger windows smooth more aggressively.")
        lbl = QtWidgets.QLabel("Savitzky-Golay window:")
        lbl.setToolTip(self.sd_sg_window.toolTip())
        form.addRow(lbl, self.sd_sg_window)

        self.sd_sg_polyorder = QtWidgets.QSpinBox()
        self.sd_sg_polyorder.setRange(2, 10)
        self.sd_sg_polyorder.setValue(3)
        self.sd_sg_polyorder.setToolTip("Polynomial order for the Savitzky-Golay filter. Must be less than the window length. Higher orders preserve sharper peak shapes.")
        lbl = QtWidgets.QLabel("SG polynomial order:")
        lbl.setToolTip(self.sd_sg_polyorder.toolTip())
        form.addRow(lbl, self.sd_sg_polyorder)

        self.sd_pre_smoothing = QtWidgets.QComboBox()
        self.sd_pre_smoothing.addItems(["none", "gaussian", "moving_average", "savitzky_golay"])
        self.sd_pre_smoothing.setToolTip("Optional smoothing applied to the EIC before computing the second derivative. Useful for noisy data.")
        lbl = QtWidgets.QLabel("Pre-smoothing method:")
        lbl.setToolTip(self.sd_pre_smoothing.toolTip())
        form.addRow(lbl, self.sd_pre_smoothing)

        self.sd_pre_smoothing_window = QtWidgets.QSpinBox()
        self.sd_pre_smoothing_window.setRange(1, 201)
        self.sd_pre_smoothing_window.setValue(5)
        self.sd_pre_smoothing_window.setToolTip("Window size for the optional pre-smoothing filter in scans.")
        lbl = QtWidgets.QLabel("Pre-smoothing window:")
        lbl.setToolTip(self.sd_pre_smoothing_window.toolTip())
        form.addRow(lbl, self.sd_pre_smoothing_window)

        self.sd_pre_smoothing_iterations = QtWidgets.QSpinBox()
        self.sd_pre_smoothing_iterations.setRange(1, 20)
        self.sd_pre_smoothing_iterations.setValue(1)
        self.sd_pre_smoothing_iterations.setToolTip("Number of times pre-smoothing is applied. More iterations give stronger smoothing.")
        lbl = QtWidgets.QLabel("Pre-smoothing iterations:")
        lbl.setToolTip(self.sd_pre_smoothing_iterations.toolTip())
        form.addRow(lbl, self.sd_pre_smoothing_iterations)

        self.sd_min_apex_intensity = QtWidgets.QDoubleSpinBox()
        self.sd_min_apex_intensity.setRange(0, 1e12)
        self.sd_min_apex_intensity.setValue(0.0)
        self.sd_min_apex_intensity.setDecimals(1)
        self.sd_min_apex_intensity.setToolTip("Minimum intensity at the peak apex. Peaks with an apex below this value are discarded.")
        lbl = QtWidgets.QLabel("Min apex intensity:")
        lbl.setToolTip(self.sd_min_apex_intensity.toolTip())
        form.addRow(lbl, self.sd_min_apex_intensity)

        self.tabWidget.addTab(tab, "Second Derivative")

    def _build_matched_filter_tab(self):
        tab = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(tab)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)

        self.mf_expected_peak_width = QtWidgets.QSpinBox()
        self.mf_expected_peak_width.setRange(3, 500)
        self.mf_expected_peak_width.setValue(15)
        self.mf_expected_peak_width.setToolTip("Expected peak width in scans, used to generate the Gaussian template for cross-correlation.")
        lbl = QtWidgets.QLabel("Expected peak width (scans):")
        lbl.setToolTip(self.mf_expected_peak_width.toolTip())
        form.addRow(lbl, self.mf_expected_peak_width)

        self.mf_correlation_threshold = QtWidgets.QDoubleSpinBox()
        self.mf_correlation_threshold.setRange(0.0, 1.0)
        self.mf_correlation_threshold.setValue(0.5)
        self.mf_correlation_threshold.setDecimals(3)
        self.mf_correlation_threshold.setSingleStep(0.05)
        self.mf_correlation_threshold.setToolTip("Minimum normalized cross-correlation score for accepting a peak. Higher values require a closer match to the template shape.")
        lbl = QtWidgets.QLabel("Correlation threshold:")
        lbl.setToolTip(self.mf_correlation_threshold.toolTip())
        form.addRow(lbl, self.mf_correlation_threshold)

        self.mf_template_type = QtWidgets.QComboBox()
        self.mf_template_type.addItems(["gaussian"])
        self.mf_template_type.setToolTip("Shape of the template used for cross-correlation. Currently only Gaussian is supported.")
        lbl = QtWidgets.QLabel("Template type:")
        lbl.setToolTip(self.mf_template_type.toolTip())
        form.addRow(lbl, self.mf_template_type)

        self.mf_min_distance = QtWidgets.QSpinBox()
        self.mf_min_distance.setRange(1, 500)
        self.mf_min_distance.setValue(5)
        self.mf_min_distance.setToolTip("Minimum distance in scans between two adjacent peaks. Peaks closer than this are merged, keeping the stronger one.")
        lbl = QtWidgets.QLabel("Min distance between peaks (scans):")
        lbl.setToolTip(self.mf_min_distance.toolTip())
        form.addRow(lbl, self.mf_min_distance)

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
        self.wt_min_scale.setValue(s.get("wt_min_scale", 3))
        self.wt_max_scale.setValue(s.get("wt_max_scale", 11))
        self.wt_num_scales.setValue(s.get("wt_num_scales", 8))
        self.wt_snr_threshold.setValue(s.get("wt_snr_threshold", 3.0))
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
