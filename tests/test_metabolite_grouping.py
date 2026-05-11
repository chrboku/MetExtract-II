from src.metaboliteGrouping import split_group_by_relative_abundance


def _normalize_groups(groups):
    return sorted([sorted(group) for group in groups])


def test_split_group_by_relative_abundance_splits_dissimilar_profiles():
    groups = split_group_by_relative_abundance(
        [1, 2, 3, 4],
        {
            1: [100, 200, 300],
            2: [50, 100, 150],
            3: [300, 200, 100],
            4: [150, 100, 50],
        },
        min_peak_correlation=0.8,
        min_connection_rate=0.6,
    )

    assert _normalize_groups(groups) == [[1, 2], [3, 4]]


def test_split_group_by_relative_abundance_keeps_group_for_missing_vectors():
    groups = split_group_by_relative_abundance(
        [10, 11, 12],
        {
            10: [100, 200, 300],
            11: [120, 220, 320],
        },
        min_peak_correlation=0.8,
        min_connection_rate=0.6,
    )

    assert groups == [[10, 11, 12]]
