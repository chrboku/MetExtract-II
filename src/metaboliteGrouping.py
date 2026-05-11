import numpy as np


def _coerce_abundance_profile(values):
    profile = []
    for value in values:
        try:
            numeric = float(value)
            if np.isnan(numeric):
                numeric = 0.0
        except (TypeError, ValueError):
            numeric = 0.0
        profile.append(max(0.0, numeric))
    return profile


def _pearson(values_a, values_b):
    if len(values_a) != len(values_b) or len(values_a) < 2:
        return None
    if np.std(values_a) == 0 or np.std(values_b) == 0:
        return None
    corr = np.corrcoef(values_a, values_b)[0, 1]
    if np.isnan(corr):
        return None
    return float(corr)


def _presence_aware_profile_similarity(profile_a, profile_b, both_missing_weight=0.25):
    if profile_a is None or profile_b is None:
        return 0.0
    if len(profile_a) < 2 or len(profile_b) < 2 or len(profile_a) != len(profile_b):
        return 0.0

    quantified_a = []
    quantified_b = []
    both_missing_count = 0

    for value_a, value_b in zip(profile_a, profile_b):
        has_a = value_a > 0.0
        has_b = value_b > 0.0
        if has_a and has_b:
            quantified_a.append(value_a)
            quantified_b.append(value_b)
        elif (not has_a) and (not has_b):
            both_missing_count += 1
        # samples where only one feature has abundance are omitted completely

    quantified_count = len(quantified_a)
    if quantified_count == 0 and both_missing_count == 0:
        return 0.0

    quantified_corr = _pearson(quantified_a, quantified_b)
    if quantified_corr is None:
        quantified_corr = 1.0 if quantified_count == 1 else 0.0

    weighted_sum = quantified_count * quantified_corr + both_missing_weight * both_missing_count
    total_weight = quantified_count + both_missing_weight * both_missing_count
    if total_weight <= 0:
        return 0.0

    return float(weighted_sum / total_weight)


def _connected_components(group_ids, adjacency):
    remaining = set(group_ids)
    components = []
    while remaining:
        root = remaining.pop()
        stack = [root]
        component = {root}
        while stack:
            node = stack.pop()
            for neighbor in adjacency.get(node, set()):
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    component.add(neighbor)
                    stack.append(neighbor)
        components.append(component)
    return components


def _split_component_by_connection_rate(component_ids, adjacency, min_connection_rate):
    if len(component_ids) <= 2:
        return [sorted(component_ids)]

    component_set = set(component_ids)
    changed = True
    while changed and len(component_set) > 2:
        changed = False
        min_required_connections = max(0.0, min_connection_rate) * (len(component_set) - 1)
        to_remove = set()
        for node in component_set:
            degree = sum(1 for neighbor in adjacency.get(node, set()) if neighbor in component_set)
            if degree < min_required_connections:
                to_remove.add(node)
        if to_remove:
            component_set -= to_remove
            changed = True

    if not component_set:
        return [[feature_id] for feature_id in sorted(component_ids)]

    result = []
    kept_components = _connected_components(component_set, adjacency)
    for kept_component in kept_components:
        if kept_component == set(component_ids):
            result.append(sorted(kept_component))
        else:
            result.extend(_split_component_by_connection_rate(sorted(kept_component), adjacency, min_connection_rate))

    removed = set(component_ids) - component_set
    if removed:
        removed_components = _connected_components(removed, adjacency)
        for removed_component in removed_components:
            if len(removed_component) <= 1:
                result.append(sorted(removed_component))
            else:
                result.extend(_split_component_by_connection_rate(sorted(removed_component), adjacency, min_connection_rate))

    return result


def split_group_by_relative_abundance(group_ids, abundance_vectors, min_peak_correlation, min_connection_rate):
    if len(group_ids) <= 2:
        return [group_ids]

    profiles = {}
    for feature_id in group_ids:
        if feature_id in abundance_vectors:
            profiles[feature_id] = _coerce_abundance_profile(abundance_vectors[feature_id])

    if len(profiles) < len(group_ids):
        return [group_ids]
    profile_lengths = {len(profile) for profile in profiles.values()}
    if len(profile_lengths) != 1:
        return [group_ids]

    similarities = {feature_id: {} for feature_id in group_ids}
    adjacency = {feature_id: set() for feature_id in group_ids}
    for i in range(len(group_ids)):
        for j in range(i + 1, len(group_ids)):
            feature_a = group_ids[i]
            feature_b = group_ids[j]
            similarities[feature_a][feature_b] = _presence_aware_profile_similarity(profiles[feature_a], profiles[feature_b])
            similarities[feature_b][feature_a] = similarities[feature_a][feature_b]
            if similarities[feature_a][feature_b] >= min_peak_correlation:
                adjacency[feature_a].add(feature_b)
                adjacency[feature_b].add(feature_a)

    components = _connected_components(group_ids, adjacency)
    refined_groups = []
    for component in components:
        refined_groups.extend(_split_component_by_connection_rate(sorted(component), adjacency, min_connection_rate))

    return sorted([sorted(group) for group in refined_groups], key=lambda group: (len(group), group), reverse=True)
