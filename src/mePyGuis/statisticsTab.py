# -*- coding: utf-8 -*-
"""
Statistics Tab GUI for MetExtract II

This module provides the GUI components for the Statistics tab, including:
- Tree view for selecting analysis methods
- Data Quality Overview visualizations
- PCA, HCA, and Heat Map plots
- Interactive Volcano Plots with rectangular selection

Copyright (C) 2015 MetExtract Team
License: GNU General Public License v2 (GPLv2)
"""

from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import Qt, Signal, QRectF
from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QSplitter,
    QTreeWidget,
    QTreeWidgetItem,
    QTableWidget,
    QTableWidgetItem,
    QDialog,
    QDialogButtonBox,
    QComboBox,
    QLabel,
    QPushButton,
    QFrame,
    QScrollArea,
    QGroupBox,
    QGridLayout,
    QHeaderView,
    QAbstractItemView,
)

import matplotlib

matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from typing import List, Dict, Tuple, Optional, Any
import logging

from .statisticsModule import StatisticsData, DataQualityAnalysis, MultivariateAnalysis, UnivariateAnalysis, SelectionManager


class AddComparisonDialog(QDialog):
    """Dialog for adding a new volcano plot comparison."""

    def __init__(self, available_groups: List[str], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Add Volcano Plot Comparison")
        self.setMinimumWidth(300)

        layout = QVBoxLayout(self)

        # Group 1 selection
        layout.addWidget(QLabel("Select first group (Control/Reference):"))
        self.group1_combo = QComboBox()
        self.group1_combo.addItems(available_groups)
        layout.addWidget(self.group1_combo)

        # Group 2 selection
        layout.addWidget(QLabel("Select second group (Treatment/Condition):"))
        self.group2_combo = QComboBox()
        self.group2_combo.addItems(available_groups)
        layout.addWidget(self.group2_combo)

        # Buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def get_groups(self) -> Tuple[str, str]:
        """Return selected groups."""
        return (self.group1_combo.currentText(), self.group2_combo.currentText())


class InteractiveVolcanoCanvas(FigureCanvas):
    """Canvas for volcano plot with interactive rectangular selection."""

    selectionChanged = Signal(list, bool)  # Signal emitted when selection changes

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)

        self.volcano_data = None
        self.scatter = None
        self.highlighted_indices = []
        self.rect_selector = None

        # Set up matplotlib event handling
        self.fig.canvas.mpl_connect("key_press_event", self._on_key_press)

        # Initialize rectangle selector
        self._setup_rect_selector()

    def _setup_rect_selector(self):
        """Set up the rectangle selector for feature selection."""
        self.rect_selector = RectangleSelector(
            self.axes,
            self._on_select,
            useblit=True,
            button=[1],  # Left mouse button
            minspanx=5,
            minspany=5,
            spancoords="pixels",
            interactive=True,
        )

    def _on_select(self, eclick, erelease):
        """Handle rectangle selection."""
        if self.volcano_data is None:
            return

        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        if x1 is None or y1 is None or x2 is None or y2 is None:
            return

        # Get data bounds
        x_min, x_max = min(x1, x2), max(x1, x2)
        y_min, y_max = min(y1, y2), max(y1, y2)

        # Find points within selection rectangle
        log2_fc = self.volcano_data["log2_fc"]
        neg_log10_pval = self.volcano_data["neg_log10_pval"]

        selected_mask = (log2_fc >= x_min) & (log2_fc <= x_max) & (neg_log10_pval >= y_min) & (neg_log10_pval <= y_max)

        selected_indices = np.where(selected_mask)[0].tolist()

        # Check if Ctrl is pressed for additive selection
        modifiers = QtWidgets.QApplication.keyboardModifiers()
        additive = modifiers == Qt.ControlModifier

        self.selectionChanged.emit(selected_indices, additive)

    def _on_key_press(self, event):
        """Handle key press events."""
        if event.key == "escape":
            self.highlighted_indices = []
            self.update_highlighting()

    def set_volcano_data(self, data: Dict[str, Any]):
        """Set the volcano plot data and redraw."""
        self.volcano_data = data
        self.draw_volcano()

    def draw_volcano(self):
        """Draw the volcano plot."""
        if self.volcano_data is None or not self.volcano_data.get("success", False):
            return

        self.axes.clear()

        log2_fc = self.volcano_data["log2_fc"]
        neg_log10_pval = self.volcano_data["neg_log10_pval"]
        significant = self.volcano_data["significant"]
        fc_threshold = self.volcano_data.get("fc_threshold", 1.0)
        pvalue_threshold = self.volcano_data.get("pvalue_threshold", 0.05)

        # Color points based on significance
        colors = []
        for i, sig in enumerate(significant):
            if i in self.highlighted_indices:
                colors.append("orange")
            elif sig:
                if log2_fc[i] > 0:
                    colors.append("red")
                else:
                    colors.append("blue")
            else:
                colors.append("gray")

        self.scatter = self.axes.scatter(log2_fc, neg_log10_pval, c=colors, alpha=0.7, s=30)

        # Add threshold lines
        self.axes.axhline(y=-np.log10(pvalue_threshold), color="gray", linestyle="--", alpha=0.5)
        self.axes.axvline(x=-fc_threshold, color="gray", linestyle="--", alpha=0.5)
        self.axes.axvline(x=fc_threshold, color="gray", linestyle="--", alpha=0.5)

        self.axes.set_xlabel("log₂(Fold Change)")
        self.axes.set_ylabel("-log₁₀(p-value)")
        self.axes.set_title("Volcano Plot")

        self.fig.tight_layout()
        self.draw()

    def update_highlighting(self, indices: Optional[List[int]] = None):
        """Update highlighted points."""
        if indices is not None:
            self.highlighted_indices = indices
        self.draw_volcano()

    def clear_highlighting(self):
        """Clear all highlighting."""
        self.highlighted_indices = []
        self.draw_volcano()


class MultiVolcanoWidget(QWidget):
    """Widget for displaying multiple volcano plots simultaneously."""

    featureSelected = Signal(list)  # Signal when features are selected

    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QGridLayout(self)
        self.volcano_canvases: List[InteractiveVolcanoCanvas] = []
        self.selection_manager = SelectionManager()
        self.selection_manager.register_callback(self._on_selection_changed)

    def set_comparisons(self, comparisons: List[Tuple[str, str]], data: Dict[str, Any]):
        """
        Set up volcano plots for all comparisons.

        Args:
            comparisons: List of (group1, group2) tuples
            data: Dictionary containing feature data and group info
        """
        # Clear existing canvases
        for canvas in self.volcano_canvases:
            canvas.deleteLater()
        self.volcano_canvases.clear()

        # Clear layout
        while self.layout.count():
            item = self.layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        if not comparisons:
            return

        # Calculate grid layout
        n_plots = len(comparisons)
        cols = min(3, n_plots)
        rows = (n_plots + cols - 1) // cols

        for i, (group1, group2) in enumerate(comparisons):
            canvas = InteractiveVolcanoCanvas(self)
            canvas.selectionChanged.connect(self._handle_selection)

            row = i // cols
            col = i % cols

            # Create container with label
            container = QFrame()
            container_layout = QVBoxLayout(container)
            label = QLabel(f"{group1} vs {group2}")
            label.setAlignment(Qt.AlignCenter)
            container_layout.addWidget(label)
            container_layout.addWidget(canvas)

            self.layout.addWidget(container, row, col)
            self.volcano_canvases.append(canvas)

            # Calculate volcano data for this comparison
            if "feature_data" in data and "group_info" in data:
                feature_data = data["feature_data"]
                group_info = data["group_info"]

                if group1 in group_info and group2 in group_info:
                    volcano_data = UnivariateAnalysis.calculate_volcano_data(feature_data, group_info[group1], group_info[group2])
                    canvas.set_volcano_data(volcano_data)

    def _handle_selection(self, indices: List[int], additive: bool):
        """Handle selection from any volcano plot."""
        self.selection_manager.add_selection(indices, additive)
        self.featureSelected.emit(self.selection_manager.get_selected_indices())

    def _on_selection_changed(self, indices: List[int]):
        """Update all volcano plots when selection changes."""
        for canvas in self.volcano_canvases:
            canvas.update_highlighting(indices)


class SelectedFeaturesTable(QTableWidget):
    """Table widget for displaying selected features."""

    viewFeatureRequested = Signal(int, str)  # Signal to view feature in experiment results

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setColumnCount(5)
        self.setHorizontalHeaderLabels(["Index", "m/z", "RT", "Log2 FC", "p-value"])
        self.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._show_context_menu)
        self.feature_data = []

    def update_features(self, indices: List[int], feature_metadata: Optional[Dict[str, Any]] = None):
        """Update table with selected features."""
        self.setRowCount(len(indices))
        self.feature_data = []

        for row, idx in enumerate(indices):
            self.setItem(row, 0, QTableWidgetItem(str(idx)))

            if feature_metadata:
                mz = feature_metadata.get("mz", {}).get(idx, "N/A")
                rt = feature_metadata.get("rt", {}).get(idx, "N/A")
                fc = feature_metadata.get("log2_fc", {}).get(idx, "N/A")
                pval = feature_metadata.get("pvalue", {}).get(idx, "N/A")
            else:
                mz, rt, fc, pval = "N/A", "N/A", "N/A", "N/A"

            self.setItem(row, 1, QTableWidgetItem(str(mz) if mz != "N/A" else mz))
            self.setItem(row, 2, QTableWidgetItem(f"{rt:.2f}" if isinstance(rt, (int, float)) else str(rt)))
            self.setItem(row, 3, QTableWidgetItem(f"{fc:.3f}" if isinstance(fc, (int, float)) else str(fc)))
            self.setItem(row, 4, QTableWidgetItem(f"{pval:.2e}" if isinstance(pval, (int, float)) else str(pval)))

            self.feature_data.append({"index": idx, "mz": mz, "rt": rt})

    def _show_context_menu(self, position):
        """Show context menu for feature actions."""
        menu = QtWidgets.QMenu(self)
        view_action = menu.addAction("View in Experiment Results")
        action = menu.exec_(self.mapToGlobal(position))

        if action == view_action:
            current_row = self.currentRow()
            if current_row >= 0 and current_row < len(self.feature_data):
                feature = self.feature_data[current_row]
                self.viewFeatureRequested.emit(feature["index"], "experiment_results")


class StatisticsCanvas(FigureCanvas):
    """Generic canvas for statistical visualizations."""

    def __init__(self, parent=None, width=6, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


class StatisticsTabWidget(QWidget):
    """Main widget for the Statistics tab."""

    # Signal to switch to experiment results and show a specific feature
    showFeatureInExperiment = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.stats_data = StatisticsData()
        self.selection_manager = SelectionManager()

        self._setup_ui()
        self._connect_signals()

    def _setup_ui(self):
        """Set up the user interface."""
        main_layout = QHBoxLayout(self)

        # Left panel: Tree view for method selection
        left_panel = QFrame()
        left_layout = QVBoxLayout(left_panel)
        left_panel.setMaximumWidth(300)

        left_layout.addWidget(QLabel("<b>Analysis Methods</b>"))

        self.methods_tree = QTreeWidget()
        self.methods_tree.setHeaderHidden(True)
        self._populate_methods_tree()
        left_layout.addWidget(self.methods_tree)

        # Add comparison button
        self.add_comparison_btn = QPushButton("Add Volcano Comparison")
        self.add_comparison_btn.setEnabled(False)
        left_layout.addWidget(self.add_comparison_btn)

        # Remove comparison button
        self.remove_comparison_btn = QPushButton("Remove Selected Comparison")
        self.remove_comparison_btn.setEnabled(False)
        left_layout.addWidget(self.remove_comparison_btn)

        main_layout.addWidget(left_panel)

        # Right panel: Content area with splitter
        right_splitter = QSplitter(Qt.Vertical)

        # Visualization area
        self.viz_container = QScrollArea()
        self.viz_container.setWidgetResizable(True)
        self.viz_widget = QWidget()
        self.viz_layout = QVBoxLayout(self.viz_widget)
        self.viz_container.setWidget(self.viz_widget)
        right_splitter.addWidget(self.viz_container)

        # Selected features table
        table_container = QFrame()
        table_layout = QVBoxLayout(table_container)
        table_layout.addWidget(QLabel("<b>Selected Features</b>"))
        self.features_table = SelectedFeaturesTable()
        table_layout.addWidget(self.features_table)

        view_feature_btn = QPushButton("View Selected in Experiment Results")
        view_feature_btn.clicked.connect(self._view_selected_feature)
        table_layout.addWidget(view_feature_btn)

        right_splitter.addWidget(table_container)
        right_splitter.setSizes([600, 200])

        main_layout.addWidget(right_splitter)

        # Store current visualization widgets
        self.current_canvas = None
        self.current_toolbar = None
        self.multi_volcano_widget = None

    def _populate_methods_tree(self):
        """Populate the methods tree with analysis categories."""
        # Data Quality Overview
        quality_item = QTreeWidgetItem(self.methods_tree, ["Data Quality Overview"])
        QTreeWidgetItem(quality_item, ["Feature Detection Counts"])
        QTreeWidgetItem(quality_item, ["RSD (Relative Standard Deviation)"])
        QTreeWidgetItem(quality_item, ["Abundance Histograms"])

        # Multivariate Analysis
        multivariate_item = QTreeWidgetItem(self.methods_tree, ["Multivariate Analysis"])
        QTreeWidgetItem(multivariate_item, ["PCA (Principal Component Analysis)"])
        QTreeWidgetItem(multivariate_item, ["HCA (Hierarchical Cluster Analysis)"])
        QTreeWidgetItem(multivariate_item, ["Heat Map"])

        # Univariate Analysis
        univariate_item = QTreeWidgetItem(self.methods_tree, ["Univariate Analysis"])
        self.volcano_parent_item = QTreeWidgetItem(univariate_item, ["Volcano Plots"])
        self.all_volcano_item = QTreeWidgetItem(self.volcano_parent_item, ["All Comparisons"])

        self.methods_tree.expandAll()

    def _connect_signals(self):
        """Connect signals and slots."""
        self.methods_tree.itemClicked.connect(self._on_method_selected)
        self.add_comparison_btn.clicked.connect(self._add_volcano_comparison)
        self.remove_comparison_btn.clicked.connect(self._remove_volcano_comparison)
        self.features_table.viewFeatureRequested.connect(self._on_view_feature_requested)
        self.selection_manager.register_callback(self._on_selection_changed)

    def load_experiment_data(self, experiment_data: Dict[str, Any]):
        """Load experiment data for statistical analysis."""
        if self.stats_data.load_from_experiment(experiment_data):
            self.add_comparison_btn.setEnabled(len(self.stats_data.get_group_names()) >= 2)
            logging.info("Statistics tab: Experiment data loaded successfully")
        else:
            logging.warning("Statistics tab: Failed to load experiment data")

    def _on_method_selected(self, item: QTreeWidgetItem, column: int):
        """Handle method selection from tree."""
        method_name = item.text(0)
        parent = item.parent()

        # Enable/disable buttons based on selection
        is_volcano_comparison = parent and parent.text(0) == "Volcano Plots" and method_name != "All Comparisons"
        self.remove_comparison_btn.setEnabled(is_volcano_comparison)

        # Show appropriate visualization
        if method_name == "Feature Detection Counts":
            self._show_detection_counts()
        elif method_name == "RSD (Relative Standard Deviation)":
            self._show_rsd_plot()
        elif method_name == "Abundance Histograms":
            self._show_abundance_histogram()
        elif method_name == "PCA (Principal Component Analysis)":
            self._show_pca()
        elif method_name == "HCA (Hierarchical Cluster Analysis)":
            self._show_hca()
        elif method_name == "Heat Map":
            self._show_heatmap()
        elif method_name == "All Comparisons":
            self._show_all_volcano_plots()
        elif is_volcano_comparison:
            self._show_single_volcano_plot(item)

    def _clear_visualization(self):
        """Clear the current visualization."""
        while self.viz_layout.count():
            item = self.viz_layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()

        self.current_canvas = None
        self.current_toolbar = None
        self.multi_volcano_widget = None

    def _show_detection_counts(self):
        """Show feature detection counts visualization."""
        self._clear_visualization()

        if self.stats_data.feature_data is None:
            self._show_no_data_message()
            return

        canvas = StatisticsCanvas(self, width=8, height=6)
        toolbar = NavigationToolbar(canvas, self)

        detection_counts = DataQualityAnalysis.calculate_detection_counts(self.stats_data.feature_data, self.stats_data.group_info)

        if detection_counts:
            ax = canvas.axes
            groups = list(detection_counts.keys())
            n_groups = len(groups)

            for i, (group_name, counts) in enumerate(detection_counts.items()):
                ax.hist(counts, bins=np.arange(0, max(counts) + 2) - 0.5, alpha=0.7, label=group_name)

            ax.set_xlabel("Number of replicates with detection")
            ax.set_ylabel("Number of features")
            ax.set_title("Feature Detection Counts by Group")
            ax.legend()
            canvas.fig.tight_layout()

        self.viz_layout.addWidget(toolbar)
        self.viz_layout.addWidget(canvas)
        self.current_canvas = canvas
        self.current_toolbar = toolbar

    def _show_rsd_plot(self):
        """Show RSD (Relative Standard Deviation) plot."""
        self._clear_visualization()

        if self.stats_data.feature_data is None:
            self._show_no_data_message()
            return

        canvas = StatisticsCanvas(self, width=8, height=6)
        toolbar = NavigationToolbar(canvas, self)

        rsd_values = DataQualityAnalysis.calculate_rsd(self.stats_data.feature_data, self.stats_data.group_info)

        if rsd_values:
            ax = canvas.axes

            data_to_plot = [rsd[~np.isnan(rsd)] for rsd in rsd_values.values()]
            labels = list(rsd_values.keys())

            bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)

            # Color the boxes
            colors = plt.cm.Set3(np.linspace(0, 1, len(labels)))
            for patch, color in zip(bp["boxes"], colors):
                patch.set_facecolor(color)

            ax.set_ylabel("RSD (%)")
            ax.set_title("Relative Standard Deviation by Group")
            ax.axhline(y=30, color="red", linestyle="--", alpha=0.7, label="30% threshold")
            ax.legend()
            canvas.fig.tight_layout()

        self.viz_layout.addWidget(toolbar)
        self.viz_layout.addWidget(canvas)
        self.current_canvas = canvas
        self.current_toolbar = toolbar

    def _show_abundance_histogram(self):
        """Show abundance histogram."""
        self._clear_visualization()

        if self.stats_data.feature_data is None:
            self._show_no_data_message()
            return

        canvas = StatisticsCanvas(self, width=8, height=6)
        toolbar = NavigationToolbar(canvas, self)

        counts, edges = DataQualityAnalysis.calculate_abundance_distribution(self.stats_data.feature_data)

        if len(counts) > 0:
            ax = canvas.axes
            ax.bar(edges[:-1], counts, width=np.diff(edges), align="edge", alpha=0.7, color="steelblue")
            ax.set_xlabel("log₁₀(Intensity)")
            ax.set_ylabel("Frequency")
            ax.set_title("Feature Abundance Distribution")
            canvas.fig.tight_layout()

        self.viz_layout.addWidget(toolbar)
        self.viz_layout.addWidget(canvas)
        self.current_canvas = canvas
        self.current_toolbar = toolbar

    def _show_pca(self):
        """Show PCA visualization."""
        self._clear_visualization()

        if self.stats_data.feature_data is None:
            self._show_no_data_message()
            return

        canvas = StatisticsCanvas(self, width=8, height=6)
        toolbar = NavigationToolbar(canvas, self)

        pca_result = MultivariateAnalysis.perform_pca(self.stats_data.feature_data)

        if pca_result.get("success", False):
            ax = canvas.axes
            scores = pca_result["scores"]
            sample_names = pca_result["sample_names"]
            var_ratio = pca_result["explained_variance_ratio"]

            # Color by group
            group_colors = {}
            color_idx = 0
            colors = plt.cm.Set1(np.linspace(0, 1, len(self.stats_data.group_info)))

            for group_name, samples in self.stats_data.group_info.items():
                group_colors[group_name] = colors[color_idx]
                color_idx += 1

            for i, sample in enumerate(sample_names):
                sample_color = "gray"
                sample_group = None
                for group_name, samples in self.stats_data.group_info.items():
                    if sample in samples:
                        sample_color = group_colors[group_name]
                        sample_group = group_name
                        break

                ax.scatter(scores[i, 0], scores[i, 1], c=[sample_color], s=100, alpha=0.8)
                ax.annotate(sample, (scores[i, 0], scores[i, 1]), fontsize=8, alpha=0.7)

            ax.set_xlabel(f"PC1 ({var_ratio[0]*100:.1f}%)")
            ax.set_ylabel(f"PC2 ({var_ratio[1]*100:.1f}%)")
            ax.set_title("PCA Score Plot")

            # Add legend
            legend_handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=10, label=group) for group, color in group_colors.items()]
            ax.legend(handles=legend_handles, loc="best")
            canvas.fig.tight_layout()

        self.viz_layout.addWidget(toolbar)
        self.viz_layout.addWidget(canvas)
        self.current_canvas = canvas
        self.current_toolbar = toolbar

    def _show_hca(self):
        """Show HCA dendrogram."""
        self._clear_visualization()

        if self.stats_data.feature_data is None:
            self._show_no_data_message()
            return

        canvas = StatisticsCanvas(self, width=10, height=6)
        toolbar = NavigationToolbar(canvas, self)

        hca_result = MultivariateAnalysis.perform_hca(self.stats_data.feature_data)

        if hca_result.get("success", False):
            ax = canvas.axes
            dendrogram(hca_result["linkage_matrix"], labels=hca_result["sample_names"], ax=ax, leaf_rotation=90)
            ax.set_ylabel("Distance")
            ax.set_title("Hierarchical Cluster Analysis")
            canvas.fig.tight_layout()

        self.viz_layout.addWidget(toolbar)
        self.viz_layout.addWidget(canvas)
        self.current_canvas = canvas
        self.current_toolbar = toolbar

    def _show_heatmap(self):
        """Show heat map visualization."""
        self._clear_visualization()

        if self.stats_data.feature_data is None:
            self._show_no_data_message()
            return

        canvas = StatisticsCanvas(self, width=10, height=8)
        toolbar = NavigationToolbar(canvas, self)

        heatmap_result = MultivariateAnalysis.prepare_heatmap_data(self.stats_data.feature_data)

        if heatmap_result.get("success", False):
            ax = canvas.axes
            data = heatmap_result["data"]

            # Limit number of features for visualization
            max_features = 100
            if data.shape[0] > max_features:
                # Select top variable features
                var = np.var(data, axis=1)
                top_indices = np.argsort(var)[-max_features:]
                data = data[top_indices, :]

            im = ax.imshow(data, aspect="auto", cmap="RdBu_r", interpolation="nearest")
            canvas.fig.colorbar(im, ax=ax, label="Z-score")

            ax.set_xlabel("Samples")
            ax.set_ylabel("Features")
            ax.set_title("Feature Heat Map (Top Variable Features)")

            # Set sample labels if not too many
            if len(heatmap_result["col_names"]) <= 20:
                ax.set_xticks(range(len(heatmap_result["col_names"])))
                ax.set_xticklabels(heatmap_result["col_names"], rotation=90, fontsize=8)

            canvas.fig.tight_layout()

        self.viz_layout.addWidget(toolbar)
        self.viz_layout.addWidget(canvas)
        self.current_canvas = canvas
        self.current_toolbar = toolbar

    def _show_all_volcano_plots(self):
        """Show all volcano plots simultaneously."""
        self._clear_visualization()

        if self.stats_data.feature_data is None or not self.stats_data.volcano_comparisons:
            self._show_no_data_message("No data or no comparisons defined")
            return

        self.multi_volcano_widget = MultiVolcanoWidget(self)
        self.multi_volcano_widget.featureSelected.connect(self._on_features_selected)

        data = {"feature_data": self.stats_data.feature_data, "group_info": self.stats_data.group_info}

        self.multi_volcano_widget.set_comparisons(self.stats_data.volcano_comparisons, data)
        self.viz_layout.addWidget(self.multi_volcano_widget)

    def _show_single_volcano_plot(self, item: QTreeWidgetItem):
        """Show a single volcano plot for the selected comparison."""
        self._clear_visualization()

        if self.stats_data.feature_data is None:
            self._show_no_data_message()
            return

        # Parse comparison from item text
        text = item.text(0)
        parts = text.split(" vs ")
        if len(parts) != 2:
            return

        group1, group2 = parts[0], parts[1]

        if group1 not in self.stats_data.group_info or group2 not in self.stats_data.group_info:
            return

        canvas = InteractiveVolcanoCanvas(self, width=8, height=6)
        canvas.selectionChanged.connect(lambda indices, additive: self.selection_manager.add_selection(indices, additive))
        toolbar = NavigationToolbar(canvas, self)

        volcano_data = UnivariateAnalysis.calculate_volcano_data(self.stats_data.feature_data, self.stats_data.group_info[group1], self.stats_data.group_info[group2])

        canvas.set_volcano_data(volcano_data)

        self.viz_layout.addWidget(toolbar)
        self.viz_layout.addWidget(canvas)
        self.current_canvas = canvas
        self.current_toolbar = toolbar

    def _show_no_data_message(self, message: str = "No experiment data loaded"):
        """Show a message when no data is available."""
        label = QLabel(message)
        label.setAlignment(Qt.AlignCenter)
        label.setStyleSheet("color: gray; font-size: 14px;")
        self.viz_layout.addWidget(label)

    def _add_volcano_comparison(self):
        """Add a new volcano plot comparison."""
        groups = self.stats_data.get_group_names()
        if len(groups) < 2:
            QtWidgets.QMessageBox.warning(self, "Cannot Add Comparison", "At least 2 groups are required for comparison.")
            return

        dialog = AddComparisonDialog(groups, self)
        if dialog.exec_() == QDialog.Accepted:
            group1, group2 = dialog.get_groups()
            if group1 == group2:
                QtWidgets.QMessageBox.warning(self, "Invalid Comparison", "Cannot compare a group with itself.")
                return

            if self.stats_data.add_volcano_comparison(group1, group2):
                # Add to tree
                comparison_item = QTreeWidgetItem(self.volcano_parent_item, [f"{group1} vs {group2}"])
                self.methods_tree.expandItem(self.volcano_parent_item)

    def _remove_volcano_comparison(self):
        """Remove the selected volcano comparison."""
        current_item = self.methods_tree.currentItem()
        if current_item and current_item.parent() == self.volcano_parent_item:
            text = current_item.text(0)
            if text != "All Comparisons":
                # Find and remove from data
                parts = text.split(" vs ")
                if len(parts) == 2:
                    for i, (g1, g2) in enumerate(self.stats_data.volcano_comparisons):
                        if g1 == parts[0] and g2 == parts[1]:
                            self.stats_data.remove_volcano_comparison(i)
                            break

                # Remove from tree
                self.volcano_parent_item.removeChild(current_item)

    def _on_features_selected(self, indices: List[int]):
        """Handle feature selection from volcano plots."""
        self.selection_manager.add_selection(indices, False)

    def _on_selection_changed(self, indices: List[int]):
        """Handle selection manager updates."""
        # Update features table
        metadata = None
        if self.stats_data.feature_metadata is not None:
            metadata = self.stats_data.feature_metadata.to_dict()
        self.features_table.update_features(indices, metadata)

    def _on_view_feature_requested(self, feature_index: int, target: str):
        """Handle request to view feature in experiment results."""
        self.showFeatureInExperiment.emit(feature_index)

    def _view_selected_feature(self):
        """View the currently selected feature in experiment results."""
        selected_indices = self.selection_manager.get_selected_indices()
        if selected_indices:
            self.showFeatureInExperiment.emit(selected_indices[0])
        else:
            QtWidgets.QMessageBox.information(self, "No Selection", "Please select a feature first.")
