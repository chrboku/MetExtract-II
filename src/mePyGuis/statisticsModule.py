# -*- coding: utf-8 -*-
"""
Statistics Module for MetExtract II

This module provides statistical analysis functionality for metabolomics data,
including:
- Data Quality Overview (feature detection, RSD plots, abundance histograms)
- Multivariate Analysis (PCA, HCA, Heat Map)
- Univariate Analysis (Volcano Plots with interactive selection)

Copyright (C) 2015 MetExtract Team
License: GNU General Public License v2 (GPLv2)
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy import stats
from typing import List, Dict, Tuple, Optional, Any
import logging


class StatisticsData:
    """Container for storing and managing experiment data for statistical analysis."""

    def __init__(self):
        self.feature_data: Optional[pd.DataFrame] = None
        self.group_info: Dict[str, List[str]] = {}  # group_name -> [sample_names]
        self.feature_metadata: Optional[pd.DataFrame] = None
        self.selected_features: List[int] = []
        self.volcano_comparisons: List[Tuple[str, str]] = []  # (group1, group2) pairs

    def load_from_experiment(self, experiment_data: Dict[str, Any]) -> bool:
        """
        Load data from experiment results.

        Args:
            experiment_data: Dictionary containing experiment feature data and metadata

        Returns:
            True if data was loaded successfully, False otherwise
        """
        try:
            if "features" in experiment_data and len(experiment_data["features"]) > 0:
                self.feature_data = pd.DataFrame(experiment_data["features"])
            if "groups" in experiment_data:
                self.group_info = experiment_data["groups"]
            if "metadata" in experiment_data:
                self.feature_metadata = pd.DataFrame(experiment_data["metadata"])
            return True
        except Exception as e:
            logging.error(f"Error loading experiment data: {e}")
            return False

    def get_group_names(self) -> List[str]:
        """Return list of available group names."""
        return list(self.group_info.keys())

    def add_volcano_comparison(self, group1: str, group2: str) -> bool:
        """
        Add a new volcano plot comparison.

        Args:
            group1: First group name
            group2: Second group name

        Returns:
            True if comparison was added, False if already exists
        """
        comparison = (group1, group2)
        if comparison not in self.volcano_comparisons and (group2, group1) not in self.volcano_comparisons:
            self.volcano_comparisons.append(comparison)
            return True
        return False

    def remove_volcano_comparison(self, index: int) -> bool:
        """Remove a volcano plot comparison by index."""
        if 0 <= index < len(self.volcano_comparisons):
            self.volcano_comparisons.pop(index)
            return True
        return False


class DataQualityAnalysis:
    """Performs data quality analysis on experiment data."""

    @staticmethod
    def calculate_detection_counts(data: pd.DataFrame, group_info: Dict[str, List[str]]) -> Dict[str, np.ndarray]:
        """
        Calculate feature detection counts per replicate for each group.

        Args:
            data: Feature intensity data (features x samples)
            group_info: Mapping of group names to sample names

        Returns:
            Dictionary mapping group names to arrays of detection counts
        """
        detection_counts = {}
        for group_name, samples in group_info.items():
            if len(samples) > 0:
                sample_cols = [col for col in data.columns if col in samples]
                if sample_cols:
                    # Count non-zero values per feature
                    counts = (data[sample_cols] > 0).sum(axis=1).values
                    detection_counts[group_name] = counts
        return detection_counts

    @staticmethod
    def calculate_rsd(data: pd.DataFrame, group_info: Dict[str, List[str]]) -> Dict[str, np.ndarray]:
        """
        Calculate Relative Standard Deviation (RSD/CV%) for each group.

        Args:
            data: Feature intensity data
            group_info: Mapping of group names to sample names

        Returns:
            Dictionary mapping group names to RSD arrays (one per feature)
        """
        rsd_values = {}
        for group_name, samples in group_info.items():
            sample_cols = [col for col in data.columns if col in samples]
            if len(sample_cols) >= 2:
                group_data = data[sample_cols]
                means = group_data.mean(axis=1)
                stds = group_data.std(axis=1)
                # Calculate RSD, avoiding division by zero
                rsd = np.where(means > 0, (stds / means) * 100, np.nan)
                rsd_values[group_name] = rsd
        return rsd_values

    @staticmethod
    def calculate_abundance_distribution(data: pd.DataFrame, bins: int = 50) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate abundance histogram data.

        Args:
            data: Feature intensity data
            bins: Number of histogram bins

        Returns:
            Tuple of (histogram counts, bin edges)
        """
        # Log transform for better visualization
        all_values = data.values.flatten()
        positive_values = all_values[all_values > 0]
        if len(positive_values) > 0:
            log_values = np.log10(positive_values)
            counts, edges = np.histogram(log_values, bins=bins)
            return counts, edges
        return np.array([]), np.array([])


class MultivariateAnalysis:
    """Performs multivariate statistical analyses (PCA, HCA, Heatmap)."""

    @staticmethod
    def perform_pca(data: pd.DataFrame, n_components: int = 2, scale: bool = True) -> Dict[str, Any]:
        """
        Perform Principal Component Analysis.

        Args:
            data: Feature intensity data (samples as rows after transpose)
            n_components: Number of principal components
            scale: Whether to standardize the data

        Returns:
            Dictionary with PCA results including scores, loadings, and explained variance
        """
        try:
            # Transpose so samples are rows, features are columns
            X = data.T.values

            # Handle missing values
            X = np.nan_to_num(X, nan=0)

            if scale:
                scaler = StandardScaler()
                X = scaler.fit_transform(X)

            n_components = min(n_components, min(X.shape))
            pca = PCA(n_components=n_components)
            scores = pca.fit_transform(X)

            return {
                "scores": scores,
                "loadings": pca.components_,
                "explained_variance_ratio": pca.explained_variance_ratio_,
                "sample_names": data.columns.tolist(),
                "feature_names": data.index.tolist(),
                "success": True,
            }
        except Exception as e:
            logging.error(f"PCA analysis failed: {e}")
            return {"success": False, "error": str(e)}

    @staticmethod
    def perform_hca(data: pd.DataFrame, method: str = "ward", metric: str = "euclidean") -> Dict[str, Any]:
        """
        Perform Hierarchical Cluster Analysis.

        Args:
            data: Feature intensity data
            method: Linkage method (ward, complete, average, single)
            metric: Distance metric

        Returns:
            Dictionary with HCA results including linkage matrix and dendrogram data
        """
        try:
            X = data.T.values
            X = np.nan_to_num(X, nan=0)

            # Standardize
            scaler = StandardScaler()
            X = scaler.fit_transform(X)

            linkage_matrix = linkage(X, method=method, metric=metric)

            return {
                "linkage_matrix": linkage_matrix,
                "sample_names": data.columns.tolist(),
                "success": True,
            }
        except Exception as e:
            logging.error(f"HCA analysis failed: {e}")
            return {"success": False, "error": str(e)}

    @staticmethod
    def prepare_heatmap_data(data: pd.DataFrame, scale_rows: bool = True, cluster_rows: bool = True, cluster_cols: bool = True) -> Dict[str, Any]:
        """
        Prepare data for heatmap visualization.

        Args:
            data: Feature intensity data
            scale_rows: Whether to z-score scale each row (feature)
            cluster_rows: Whether to cluster features
            cluster_cols: Whether to cluster samples

        Returns:
            Dictionary with prepared heatmap data including clustering info
        """
        try:
            heatmap_data = data.values.copy()
            heatmap_data = np.nan_to_num(heatmap_data, nan=0)

            if scale_rows:
                # Z-score scaling per row
                row_means = np.mean(heatmap_data, axis=1, keepdims=True)
                row_stds = np.std(heatmap_data, axis=1, keepdims=True)
                row_stds[row_stds == 0] = 1  # Avoid division by zero
                heatmap_data = (heatmap_data - row_means) / row_stds

            result = {
                "data": heatmap_data,
                "row_names": data.index.tolist(),
                "col_names": data.columns.tolist(),
                "success": True,
            }

            if cluster_rows and heatmap_data.shape[0] > 1:
                row_linkage = linkage(heatmap_data, method="average")
                result["row_linkage"] = row_linkage

            if cluster_cols and heatmap_data.shape[1] > 1:
                col_linkage = linkage(heatmap_data.T, method="average")
                result["col_linkage"] = col_linkage

            return result
        except Exception as e:
            logging.error(f"Heatmap preparation failed: {e}")
            return {"success": False, "error": str(e)}


class UnivariateAnalysis:
    """Performs univariate statistical analyses (Volcano Plot)."""

    @staticmethod
    def calculate_volcano_data(data: pd.DataFrame, group1_samples: List[str], group2_samples: List[str], fc_threshold: float = 1.0, pvalue_threshold: float = 0.05) -> Dict[str, Any]:
        """
        Calculate data for volcano plot visualization.

        Args:
            data: Feature intensity data
            group1_samples: Sample names for first group
            group2_samples: Sample names for second group
            fc_threshold: Log2 fold change threshold for significance
            pvalue_threshold: P-value threshold for significance

        Returns:
            Dictionary with volcano plot data including fold changes, p-values, and significance flags
        """
        try:
            # Get data for each group
            g1_cols = [col for col in data.columns if col in group1_samples]
            g2_cols = [col for col in data.columns if col in group2_samples]

            if len(g1_cols) < 2 or len(g2_cols) < 2:
                return {"success": False, "error": "Need at least 2 samples per group for statistical testing"}

            group1_data = data[g1_cols]
            group2_data = data[g2_cols]

            # Calculate means
            g1_means = group1_data.mean(axis=1)
            g2_means = group2_data.mean(axis=1)

            # Calculate log2 fold change
            # Add small value to avoid log(0)
            epsilon = 1e-10
            log2_fc = np.log2((g2_means + epsilon) / (g1_means + epsilon))

            # Perform t-test
            pvalues = []
            for i in range(len(data)):
                g1_vals = group1_data.iloc[i].values
                g2_vals = group2_data.iloc[i].values

                # Filter out zeros for comparison
                g1_valid = g1_vals[g1_vals > 0]
                g2_valid = g2_vals[g2_vals > 0]

                if len(g1_valid) >= 2 and len(g2_valid) >= 2:
                    _, pval = stats.ttest_ind(g1_valid, g2_valid)
                    pvalues.append(pval)
                else:
                    pvalues.append(1.0)

            pvalues = np.array(pvalues)
            neg_log10_pval = -np.log10(pvalues + epsilon)

            # Determine significance
            significant = (np.abs(log2_fc) >= fc_threshold) & (pvalues <= pvalue_threshold)

            return {
                "log2_fc": log2_fc.values,
                "pvalues": pvalues,
                "neg_log10_pval": neg_log10_pval,
                "significant": significant,
                "feature_indices": np.arange(len(data)),
                "feature_names": data.index.tolist() if hasattr(data.index, "tolist") else list(range(len(data))),
                "fc_threshold": fc_threshold,
                "pvalue_threshold": pvalue_threshold,
                "success": True,
            }
        except Exception as e:
            logging.error(f"Volcano plot calculation failed: {e}")
            return {"success": False, "error": str(e)}


class SelectionManager:
    """Manages feature selection across multiple visualizations."""

    def __init__(self):
        self.selected_indices: List[int] = []
        self.selection_callbacks: List[callable] = []

    def add_selection(self, indices: List[int], additive: bool = False):
        """
        Add feature indices to selection.

        Args:
            indices: List of feature indices to select
            additive: If True, add to existing selection; if False, replace selection
        """
        if additive:
            for idx in indices:
                if idx not in self.selected_indices:
                    self.selected_indices.append(idx)
        else:
            self.selected_indices = list(indices)

        self._notify_callbacks()

    def clear_selection(self):
        """Clear all selections."""
        self.selected_indices = []
        self._notify_callbacks()

    def register_callback(self, callback: callable):
        """Register a callback to be notified of selection changes."""
        if callback not in self.selection_callbacks:
            self.selection_callbacks.append(callback)

    def unregister_callback(self, callback: callable):
        """Unregister a selection callback."""
        if callback in self.selection_callbacks:
            self.selection_callbacks.remove(callback)

    def _notify_callbacks(self):
        """Notify all registered callbacks of selection change."""
        for callback in self.selection_callbacks:
            try:
                callback(self.selected_indices)
            except Exception as e:
                logging.error(f"Selection callback error: {e}")

    def get_selected_indices(self) -> List[int]:
        """Return currently selected feature indices."""
        return self.selected_indices.copy()
