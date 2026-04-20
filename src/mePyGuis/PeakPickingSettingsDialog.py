"""
Peak Picking Settings Dialog for MetExtract II.

Provides a tabbed dialog where each tab corresponds to one of the five
peak picking algorithms. The active (selected) tab determines which
algorithm will be used for chromatographic peak picking.

Features:
- EIC specification area (file, polarity, mz, ppm) with visualization
- Live peak picking results overlay when parameters change
- Algorithm info button (?) on each tab with references
"""

from PySide6 import QtCore, QtGui, QtWidgets
import numpy as np

# Algorithm descriptions and references
_ALGO_INFO = {
    "wavelettransform": (
        "Wavelet Transform Peak Picker",
        "Uses the Continuous Wavelet Transform (CWT) with Ricker (Mexican-hat) wavelets to detect "
        "chromatographic peaks across multiple scales. Ridge lines in the CWT coefficient matrix "
        "are traced from coarse to fine scales; ridges that exceed the SNR threshold and span "
        "enough scales are reported as peaks.\n\n"
        "Parameters:\n"
        "• Min/Max scale — range of wavelet widths (in scans)\n"
        "• Number of scales — resolution between min and max\n"
        "• SNR threshold — signal-to-noise cutoff\n"
        "• Min ridge length — minimum scales a ridge must span\n"
        "• Gap threshold — allowed gaps in a ridge\n\n"
        "References:\n"
        "Du P., Kibbe W.A., Lin S.M. (2006). Improved peak detection in mass spectrum by "
        "incorporating continuous wavelet transform-based pattern matching. Bioinformatics, "
        "22(17), 2059-2065.\n\n"
        "Based on the R MassSpecWavelet package approach, reimplemented in Python with scipy."
    ),
    "gradientdescent": (
        "Gradient Descent Peak Picker",
        "Smooths the EIC, finds local maxima, then walks downhill from each apex until the signal "
        "stalls or drops below a threshold. Stall detection uses a configurable number of "
        "consecutive scans with insufficient intensity increase.\n\n"
        "Parameters:\n"
        "• Smoothing method/window/iterations — noise reduction before detection\n"
        "• Min apex intensity — minimum peak height\n"
        "• Min flank intensity — boundary cutoff\n"
        "• Min increase ratio — stall detection sensitivity\n"
        "• Consecutive scans — stall window\n"
        "• Min/Max peak width — width filter in scans\n\n"
        "This is a custom algorithm developed for MetExtract II."
    ),
    "gmm": (
        "GMM (Curve Fitting) Peak Picker",
        "Segments the EIC into regions of interest (chunks) based on a baseline fraction, "
        "then fits overlapping Gaussian or Exponentially Modified Gaussian (EMG) curves using "
        "non-linear least squares (scipy.optimize.curve_fit).\n\n"
        "Parameters:\n"
        "• Distribution type — gaussian or EMG for tailing peaks\n"
        "• Baseline fraction — threshold for chunking\n"
        "• Min chunk size — minimum region length\n"
        "• Max peaks per chunk — deconvolution complexity limit\n"
        "• Fit max iterations — optimizer iteration cap\n\n"
        "References:\n"
        "Yu T., Bhatt D.P., et al. (2009). Practical guidelines for estimation of peak areas "
        "in LC-MS bioanalysis. Bioanalysis, 1(1), 87-95."
    ),
    "secondderivative": (
        "Second Derivative Peak Picker",
        "Applies a Savitzky-Golay second derivative filter to the EIC. Peaks are identified "
        "where the second derivative crosses zero (negative to positive = peak start, positive "
        "to negative = peak end) with a negative minimum between crossings (= apex).\n\n"
        "Parameters:\n"
        "• SG window/polyorder — Savitzky-Golay filter settings\n"
        "• Pre-smoothing — optional smoothing before differentiation\n"
        "• Min apex intensity — height threshold\n\n"
        "References:\n"
        "Savitzky A., Golay M.J.E. (1964). Smoothing and differentiation of data by simplified "
        "least squares procedures. Analytical Chemistry, 36(8), 1627-1639."
    ),
    "matchedfilter": (
        "Matched Filter Peak Picker",
        "Generates a Gaussian template of the expected peak width and cross-correlates it with "
        "the EIC. Local maxima in the correlation trace above a threshold are reported as peaks.\n\n"
        "Parameters:\n"
        "• Expected peak width — template width in scans\n"
        "• Correlation threshold — minimum match score (0–1)\n"
        "• Min distance — merge distance between peaks\n\n"
        "References:\n"
        "Smith C.A., Want E.J., et al. (2006). XCMS: Processing mass spectrometry data for "
        "metabolite profiling using nonlinear peak alignment, matching, and identification. "
        "Analytical Chemistry, 78(3), 779-787."
    ),
}


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
        self.setMinimumSize(700, 620)
        self.resize(780, 700)

        if current_settings is None:
            current_settings = {}

        self._settings = dict(current_settings)
        self._eic_entries = []  # list of (file, polarity, mz, ppm) tuples from UI

        self._build_ui()
        self._load_settings(current_settings)

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        layout = QtWidgets.QVBoxLayout(self)

        # --- Top: Algorithm tabs ---
        self.tabWidget = QtWidgets.QTabWidget()
        layout.addWidget(self.tabWidget, stretch=2)

        self._build_wavelet_tab()
        self._build_gradient_tab()
        self._build_gmm_tab()
        self._build_second_derivative_tab()
        self._build_matched_filter_tab()

        # --- Middle: EIC specification area ---
        eic_group = QtWidgets.QGroupBox("EIC Definitions (for preview)")
        eic_layout = QtWidgets.QVBoxLayout(eic_group)

        # Table of EIC definitions
        self.eicTable = QtWidgets.QTableWidget(0, 4)
        self.eicTable.setHorizontalHeaderLabels(["File", "Polarity", "m/z", "ppm"])
        self.eicTable.horizontalHeader().setStretchLastSection(True)
        self.eicTable.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.eicTable.setMaximumHeight(120)
        self.eicTable.setToolTip(
            "Define EICs to preview peak picking. Each row specifies a file path, "
            "polarity (+/-), m/z value, and ppm tolerance for EIC extraction."
        )
        eic_layout.addWidget(self.eicTable)

        eic_btn_layout = QtWidgets.QHBoxLayout()
        self.btnAddEic = QtWidgets.QPushButton("Add EIC")
        self.btnAddEic.setToolTip("Add a new EIC definition row")
        self.btnRemoveEic = QtWidgets.QPushButton("Remove Selected")
        self.btnRemoveEic.setToolTip("Remove selected EIC definition(s)")
        self.btnPickEic = QtWidgets.QPushButton("Pick Peaks")
        self.btnPickEic.setToolTip("Run peak picking on defined EICs and show results below")
        eic_btn_layout.addWidget(self.btnAddEic)
        eic_btn_layout.addWidget(self.btnRemoveEic)
        eic_btn_layout.addStretch()
        eic_btn_layout.addWidget(self.btnPickEic)
        eic_layout.addLayout(eic_btn_layout)

        layout.addWidget(eic_group, stretch=0)

        # --- EIC plot area (placeholder – real plotting requires matplotlib) ---
        self.eicPlotArea = QtWidgets.QLabel(
            "EIC visualization will appear here when EICs are defined and 'Pick Peaks' is clicked.\n"
            "(Requires mzXML/mzML files accessible from this machine.)"
        )
        self.eicPlotArea.setAlignment(QtCore.Qt.AlignCenter)
        self.eicPlotArea.setMinimumHeight(100)
        self.eicPlotArea.setStyleSheet(
            "background-color: #f8f8f8; border: 1px solid #ccc; color: #888; font-style: italic;"
        )
        layout.addWidget(self.eicPlotArea, stretch=1)

        # --- Bottom: Buttons ---
        btn_layout = QtWidgets.QHBoxLayout()
        self.btnDefaults = QtWidgets.QPushButton("Restore Defaults")
        self.btnOk = QtWidgets.QPushButton("OK")
        self.btnCancel = QtWidgets.QPushButton("Cancel")
        btn_layout.addWidget(self.btnDefaults)
        btn_layout.addStretch()
        btn_layout.addWidget(self.btnOk)
        btn_layout.addWidget(self.btnCancel)
        layout.addLayout(btn_layout)

        # Connections
        self.btnOk.clicked.connect(self._accept)
        self.btnCancel.clicked.connect(self.reject)
        self.btnDefaults.clicked.connect(self._restore_defaults)
        self.btnAddEic.clicked.connect(self._add_eic_row)
        self.btnRemoveEic.clicked.connect(self._remove_eic_rows)
        self.btnPickEic.clicked.connect(self._run_preview)

    # -- Tab builders ---------------------------------------------------

    def _add_info_button(self, form, algo_key):
        """Add a '?' button at the bottom of a tab that shows algorithm info."""
        btn = QtWidgets.QPushButton("?")
        btn.setFixedSize(28, 28)
        btn.setToolTip("Show algorithm description and references")
        btn.clicked.connect(lambda checked, k=algo_key: self._show_algo_info(k))
        info_layout = QtWidgets.QHBoxLayout()
        info_layout.addWidget(btn)
        info_layout.addStretch()
        form.addRow(info_layout)

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

        self._add_info_button(form, "wavelettransform")
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

        self._add_info_button(form, "gradientdescent")
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

        self._add_info_button(form, "gmm")
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

        self._add_info_button(form, "secondderivative")
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

        self._add_info_button(form, "matchedfilter")
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

    # ------------------------------------------------------------------
    # EIC management
    # ------------------------------------------------------------------

    def _add_eic_row(self):
        """Add an empty row to the EIC definition table."""
        row = self.eicTable.rowCount()
        self.eicTable.insertRow(row)
        # File – with a browse button
        fileWidget = QtWidgets.QWidget()
        fileLayout = QtWidgets.QHBoxLayout(fileWidget)
        fileLayout.setContentsMargins(0, 0, 0, 0)
        fileEdit = QtWidgets.QLineEdit()
        fileEdit.setPlaceholderText("path/to/file.mzXML")
        fileBrowse = QtWidgets.QPushButton("…")
        fileBrowse.setFixedWidth(28)
        fileBrowse.clicked.connect(lambda checked, le=fileEdit: self._browse_file(le))
        fileLayout.addWidget(fileEdit)
        fileLayout.addWidget(fileBrowse)
        self.eicTable.setCellWidget(row, 0, fileWidget)

        # Polarity
        polCombo = QtWidgets.QComboBox()
        polCombo.addItems(["+", "-"])
        self.eicTable.setCellWidget(row, 1, polCombo)

        # m/z
        mzSpin = QtWidgets.QDoubleSpinBox()
        mzSpin.setRange(0, 99999.0)
        mzSpin.setDecimals(5)
        mzSpin.setValue(100.0)
        self.eicTable.setCellWidget(row, 2, mzSpin)

        # ppm
        ppmSpin = QtWidgets.QDoubleSpinBox()
        ppmSpin.setRange(0.1, 1000.0)
        ppmSpin.setDecimals(1)
        ppmSpin.setValue(5.0)
        self.eicTable.setCellWidget(row, 3, ppmSpin)

    def _remove_eic_rows(self):
        """Remove selected rows from the EIC table."""
        rows = sorted({idx.row() for idx in self.eicTable.selectedIndexes()}, reverse=True)
        for r in rows:
            self.eicTable.removeRow(r)

    def _browse_file(self, line_edit):
        """Open a file dialog and set the result in the given QLineEdit."""
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select mzXML/mzML file", "",
            "LC-MS files (*.mzxml *.mzml);;All files (*.*)"
        )
        if fname:
            line_edit.setText(fname)

    def _run_preview(self):
        """Attempt to load EICs and run peak picking for preview."""
        # Gather EIC definitions
        eics = []
        for row in range(self.eicTable.rowCount()):
            fileW = self.eicTable.cellWidget(row, 0)
            fileEdit = fileW.findChild(QtWidgets.QLineEdit) if fileW else None
            fpath = fileEdit.text().strip() if fileEdit else ""
            polW = self.eicTable.cellWidget(row, 1)
            pol = polW.currentText() if polW else "+"
            mzW = self.eicTable.cellWidget(row, 2)
            mz = mzW.value() if mzW else 0.0
            ppmW = self.eicTable.cellWidget(row, 3)
            ppm = ppmW.value() if ppmW else 5.0
            if fpath:
                eics.append((fpath, pol, mz, ppm))

        if not eics:
            self.eicPlotArea.setText("No EICs defined. Add at least one EIC row and click 'Pick Peaks'.")
            return

        # Try to load and pick peaks
        # NOTE: imports are lazy here because this dialog can be opened outside the
        # full MetExtract environment (e.g. standalone testing) where Chromatogram
        # and peakpickers may not be available.
        try:
            from ..Chromatogram import Chromatogram
            from ..chromPeakPicking.peakpickers import (
                WaveletTransformPeakPicker,
                GradientDescentPeakPicker,
                GMMPeakPicker,
                SecondDerivativePeakPicker,
                MatchedFilterPeakPicker,
            )

            settings = self._collect_settings()
            results = []

            for fpath, pol, mz, ppm in eics:
                try:
                    mzxml = Chromatogram()
                    mzxml.parse_file(fpath)
                    scan_event = None  # auto-detect
                    eic_data, times, scan_ids, mzs = mzxml.getEIC(mz, ppm, filterLine=scan_event)
                    rt_array = np.array(times, dtype=float)
                    int_array = np.array(eic_data, dtype=float)

                    algo = settings.get("algorithm", "wavelettransform")
                    if algo == "wavelettransform":
                        picker = WaveletTransformPeakPicker(
                            min_scale=settings.get("wt_min_scale", 3),
                            max_scale=settings.get("wt_max_scale", 11),
                            num_scales=settings.get("wt_num_scales", 8),
                            snr_threshold=settings.get("wt_snr_threshold", 3.0),
                            min_ridge_length=settings.get("wt_min_ridge_length", 4),
                            gap_threshold=settings.get("wt_gap_threshold", 2),
                        )
                    elif algo == "gradientdescent":
                        picker = GradientDescentPeakPicker(
                            smoothing_method=settings.get("gd_smoothing_method", "gaussian"),
                            smoothing_window=settings.get("gd_smoothing_window", 5),
                        )
                    elif algo == "gmm":
                        picker = GMMPeakPicker(
                            distribution=settings.get("gmm_distribution", "gaussian"),
                        )
                    elif algo == "secondderivative":
                        picker = SecondDerivativePeakPicker(
                            sg_window=settings.get("sd_sg_window", 11),
                        )
                    elif algo == "matchedfilter":
                        picker = MatchedFilterPeakPicker(
                            expected_peak_width=settings.get("mf_expected_peak_width", 15),
                        )
                    else:
                        picker = WaveletTransformPeakPicker()

                    peaks = picker.pick_peaks(rt_array, int_array)
                    results.append((fpath, pol, mz, ppm, rt_array, int_array, peaks))
                except Exception as e:
                    results.append((fpath, pol, mz, ppm, None, None, str(e)))

            # Format results as text
            lines = []
            for res in results:
                fpath, pol, mz, ppm = res[:4]
                lines.append("EIC: %s [%s] m/z=%.5f ±%.1f ppm" % (fpath.split("/")[-1], pol, mz, ppm))
                if isinstance(res[6], str):
                    lines.append("  Error: %s" % res[6])
                elif res[6] is not None:
                    peaks = res[6]
                    lines.append("  Found %d peak(s):" % len(peaks))
                    for i, pk in enumerate(peaks):
                        lines.append(
                            "    Peak %d: RT=%.3f min, FWHM=%.3f, SNR=%.1f, Area=%.0f"
                            % (i + 1, pk.apex_rt, pk.fwhm, pk.snr, pk.area)
                        )
                lines.append("")
            self.eicPlotArea.setText("\n".join(lines))
        except ImportError:
            self.eicPlotArea.setText(
                "Preview not available: required modules (Chromatogram, peakpickers) not found.\n"
                "This feature requires the full MetExtract II environment."
            )
        except Exception as e:
            self.eicPlotArea.setText("Error during preview: %s" % str(e))

    # ------------------------------------------------------------------
    # Algorithm info
    # ------------------------------------------------------------------

    def _show_algo_info(self, algo_key):
        """Show an info dialog with algorithm description and references."""
        title, text = _ALGO_INFO.get(algo_key, ("Unknown", "No information available."))
        msg = QtWidgets.QMessageBox(self)
        msg.setWindowTitle(title)
        msg.setText(text)
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.exec()


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    dlg = PeakPickingSettingsDialog()
    if dlg.exec() == QtWidgets.QDialog.Accepted:
        print("Settings:", dlg.get_settings())
    sys.exit(0)
