import numpy as np
from . import HCA_general


def split_group_with_hca(group_ids, similarities, min_peak_correlation, min_connection_rate):
    if len(group_ids) <= 2:
        return [group_ids]

    data = []
    for feature_a in group_ids:
        row = []
        for feature_b in group_ids:
            if feature_a in similarities.keys() and feature_b in similarities[feature_a].keys():
                row.append(similarities[feature_a][feature_b])
            elif feature_a == feature_b:
                row.append(1.0)
            else:
                row.append(0.0)
        data.append(row)

    hca = HCA_general.HCA_generic()
    tree = hca.generateTree(objs=data, ids=group_ids)

    def check_sub_cluster(tree, hca, corr_threshold, cut_off_min_ratio):
        if isinstance(tree, HCA_general.HCALeaf):
            return False
        elif isinstance(tree, HCA_general.HCAComposite):
            corrs = hca.getLinkFor(tree)
            inds = hca.getIndsFor(tree)
            if len(inds) == 0:
                return False
            inds_set = set(inds)

            above_threshold = sum([corr > corr_threshold for i, corr in enumerate(corrs) if i in inds_set])
            return not (above_threshold * 1.0 / len(inds)) >= cut_off_min_ratio
        else:
            raise RuntimeError("Unknown if-branch")

    sub_clusters = hca.splitTreeWithCallback(tree, lambda cluster, hca: check_sub_cluster(cluster, hca, min_peak_correlation, min_connection_rate), recursive=True)
    return [[leaf.getID() for leaf in sub_cluster.getLeaves()] for sub_cluster in sub_clusters]


def _normalize_relative_abundance_profile(values):
    profile = []
    for value in values:
        try:
            numeric = float(value)
            if np.isnan(numeric):
                numeric = 0.0
        except Exception:
            numeric = 0.0
        profile.append(max(0.0, numeric))

    total = sum(profile)
    if total > 0:
        profile = [value / total for value in profile]
    return profile


def _profile_correlation(profile_a, profile_b):
    if profile_a is None or profile_b is None:
        return 0.0
    if len(profile_a) < 2 or len(profile_b) < 2 or len(profile_a) != len(profile_b):
        return 0.0
    if np.std(profile_a) == 0 or np.std(profile_b) == 0:
        return 0.0

    corr = np.corrcoef(profile_a, profile_b)[0, 1]
    if np.isnan(corr):
        return 0.0
    return float(corr)


def split_group_by_relative_abundance(group_ids, abundance_vectors, min_peak_correlation, min_connection_rate):
    if len(group_ids) <= 2:
        return [group_ids]

    profiles = {}
    for feature_id in group_ids:
        if feature_id in abundance_vectors:
            profiles[feature_id] = _normalize_relative_abundance_profile(abundance_vectors[feature_id])

    if len(profiles) < len(group_ids):
        return [group_ids]
    profile_lengths = {len(profile) for profile in profiles.values()}
    if len(profile_lengths) != 1:
        return [group_ids]

    similarities = {feature_id: {} for feature_id in group_ids}
    for i in range(len(group_ids)):
        for j in range(i + 1, len(group_ids)):
            feature_a = group_ids[i]
            feature_b = group_ids[j]
            similarities[feature_a][feature_b] = _profile_correlation(profiles[feature_a], profiles[feature_b])
            similarities[feature_b][feature_a] = similarities[feature_a][feature_b]

    return split_group_with_hca(
        group_ids,
        similarities,
        min_peak_correlation=min_peak_correlation,
        min_connection_rate=min_connection_rate,
    )
