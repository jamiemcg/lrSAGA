import statistics
import networkx as nx

def get_coverage(segment):
    return segment.coverage(count_tag = "FC")

def _distribution_summary(values):
    if len(values) == 0:
        return {
            "count": 0,
            "min": None,
            "max": None,
            "mean": None,
            "median": None,
            "p10": None,
            "p50": None,
            "p90": None,
            "p99": None,
        }

    values_sorted = sorted(values)

    def _percentile(percent):
        index = int((percent / 100) * len(values_sorted))
        index = max(0, min(index, len(values_sorted) - 1))
        return values_sorted[index]

    return {
        "count": len(values_sorted),
        "min": min(values_sorted),
        "max": max(values_sorted),
        "mean": statistics.mean(values_sorted),
        "median": statistics.median(values_sorted),
        "p10": _percentile(10),
        "p50": _percentile(50),
        "p90": _percentile(90),
        "p99": _percentile(99),
    }

def _safe_get_coverage(segment):
    try:
        coverage = get_coverage(segment)
    except Exception:
        return None

    if coverage is None:
        return None

    return float(coverage)

def coverage_summary(graph):
    node_coverages = []
    weighted_pairs = []
    total_length = 0
    weighted_sum = 0.0
    missing_coverage_nodes = 0

    for segment in graph.segments:
        coverage = _safe_get_coverage(segment)

        if coverage is None:
            missing_coverage_nodes += 1
            continue

        node_coverages.append(coverage)

        segment_length = int(segment.length)
        total_length += segment_length
        weighted_sum += coverage * segment_length
        weighted_pairs.append((coverage, segment_length))

    if len(node_coverages) == 0:
        return {
            "nodes_with_fc": 0,
            "nodes_missing_fc": missing_coverage_nodes,
            "mean": None,
            "median": None,
            "stdev": None,
            "min": None,
            "q1": None,
            "q3": None,
            "max": None,
            "length_weighted_mean": None,
            "length_weighted_median": None
        }

    node_coverages = sorted(node_coverages)

    if len(node_coverages) > 1:
        stdev = statistics.stdev(node_coverages)
    else:
        stdev = 0.0

    if len(node_coverages) > 1:
        q1, _, q3 = statistics.quantiles(node_coverages, n=4, method="inclusive")
    else:
        q1 = node_coverages[0]
        q3 = node_coverages[0]

    if total_length > 0 and len(weighted_pairs) > 0:
        length_weighted_mean = weighted_sum / total_length

        weighted_pairs = sorted(weighted_pairs, key=lambda pair: pair[0])
        half_point = total_length / 2
        running_length = 0
        length_weighted_median = weighted_pairs[-1][0]

        for coverage, segment_length in weighted_pairs:
            running_length += segment_length
            if running_length >= half_point:
                length_weighted_median = coverage
                break
    else:
        length_weighted_mean = None
        length_weighted_median = None

    return {
        "nodes_with_fc": len(node_coverages),
        "nodes_missing_fc": missing_coverage_nodes,
        "mean": statistics.mean(node_coverages),
        "median": statistics.median(node_coverages),
        "stdev": stdev,
        "min": min(node_coverages),
        "q1": q1,
        "q3": q3,
        "max": max(node_coverages),
        "length_weighted_mean": length_weighted_mean,
        "length_weighted_median": length_weighted_median
    }

def segment_stats(graph):
    segments = len(graph.segments)
    edges = len(graph.edges)
    dead_ends = graph.n_dead_ends

    total_length = 0
    segment_lengths = []

    for segment in graph.segments:
        segment_lengths.append(segment.length)
        total_length += segment.length

    segment_lengths = sorted(segment_lengths, reverse = True)

    average_segment_length = int(round(statistics.mean(segment_lengths), 0))
    median_segment_length = int(round(statistics.median(segment_lengths), 0))
    longest_segment = max(segment_lengths)
    shortest_segment = min(segment_lengths)

    tmp_length = 0
    half_total_length = total_length / 2

    n50 = 0

    for segment_length in segment_lengths:
        tmp_length += segment_length

        if tmp_length >= half_total_length:
            n50 = segment_length
            break

    line = (
        f"Nodes: {segments:,}; "
        f"Edges: {edges:,}; "
        f"Dead ends: {dead_ends:,}; "
        f"Total length: {total_length:,} bp; "
        f"N50: {n50:,} bp; "
        f"Average node length: {average_segment_length:,} bp; "
        f"Median node length: {median_segment_length:,} bp; "
        f"Longest node: {longest_segment:,} bp; "
        f"Shortest node: {shortest_segment:,} bp;"
    )
    print(line)

def connectivity_summary(graph):
    connectivity_counts = {}

    for segment in graph.segments:
        neighbour_count = len(set(neighbour.name for neighbour in segment.neighbours))
        connectivity_counts[neighbour_count] = connectivity_counts.get(neighbour_count, 0) + 1

    return connectivity_counts

def connectivity_stats(graph):
    connectivity_counts = connectivity_summary(graph)

    if len(connectivity_counts) == 0:
        print("No nodes found.")
        return

    for neighbour_count in sorted(connectivity_counts):
        node_count = connectivity_counts[neighbour_count]
        print(f"{node_count:,} nodes have {neighbour_count:,} neighbours")

def length_percentiles_and_n_stats(graph):
    segment_lengths = sorted([segment.length for segment in graph.segments], reverse=True)

    if len(segment_lengths) == 0:
        return {}

    total_length = sum(segment_lengths)

    segment_lengths_sorted_asc = sorted(segment_lengths)

    percentile_points = [1, 5, 10, 25, 50, 75, 90, 95, 99]
    percentiles = {}

    for p in percentile_points:
        index = int((p / 100) * len(segment_lengths_sorted_asc))
        index = max(0, min(index, len(segment_lengths_sorted_asc) - 1))
        percentiles[f"p{p}"] = segment_lengths_sorted_asc[index]

    n_stats = {}
    for n in [10, 25, 50, 75, 90]:
        target_length = (n / 100) * total_length
        cumulative = 0

        for length in segment_lengths:
            cumulative += length
            if cumulative >= target_length:
                n_stats[f"N{n}"] = length
                break
        else:
            n_stats[f"N{n}"] = segment_lengths[-1]

    return {**percentiles, **n_stats}

def coverage_percentiles(graph):
    coverages = []

    for segment in graph.segments:
        coverage = _safe_get_coverage(segment)
        if coverage is not None:
            coverages.append(coverage)

    if len(coverages) == 0:
        return {}

    coverages_sorted = sorted(coverages)
    percentile_points = [1, 5, 10, 25, 50, 75, 90, 95, 99]
    percentiles = {}

    for p in percentile_points:
        index = int((p / 100) * len(coverages_sorted))
        index = max(0, min(index, len(coverages_sorted) - 1))
        percentiles[f"p{p}"] = coverages_sorted[index]

    return percentiles

def component_summary(graph):
    nx_graph = nx.Graph()
    node_lengths = {}

    for segment in graph.segments:
        node_name = segment.name
        node_length = int(segment.length) if segment.length is not None else 0
        node_lengths[node_name] = node_length
        nx_graph.add_node(node_name)

    for edge in graph.edges:
        nx_graph.add_edge(edge.from_name, edge.to_name)

    components = list(nx.connected_components(nx_graph))

    if len(components) == 0:
        return {
            "n_components": 0,
            "node_size_summary": _distribution_summary([]),
            "bp_size_summary": _distribution_summary([]),
            "largest_component_nodes": 0,
            "largest_component_nodes_pct": 0.0,
            "largest_component_bp": 0,
            "largest_component_bp_pct": 0.0,
            "top_components": [],
        }

    node_sizes = [len(component) for component in components]
    bp_sizes = [sum(node_lengths[node] for node in component) for component in components]

    total_nodes = sum(node_sizes)
    total_bp = sum(bp_sizes)

    largest_component_nodes = max(node_sizes)
    largest_component_bp = max(bp_sizes)

    top_components = sorted(
        zip(node_sizes, bp_sizes),
        key=lambda values: (values[1], values[0]),
        reverse=True,
    )

    return {
        "n_components": len(components),
        "node_size_summary": _distribution_summary(node_sizes),
        "bp_size_summary": _distribution_summary(bp_sizes),
        "largest_component_nodes": largest_component_nodes,
        "largest_component_nodes_pct": (100.0 * largest_component_nodes / total_nodes) if total_nodes > 0 else 0.0,
        "largest_component_bp": largest_component_bp,
        "largest_component_bp_pct": (100.0 * largest_component_bp / total_bp) if total_bp > 0 else 0.0,
        "top_components": top_components,
    }

def coverage_stats(graph):
    summary = coverage_summary(graph)

    line = (
        f"Nodes with FC: {summary['nodes_with_fc']:,}; "
        f"Nodes missing FC: {summary['nodes_missing_fc']:,}; "
        f"Mean FC: {summary['mean']:.3f}x; "
        f"Median FC: {summary['median']:.3f}x; "
        f"Stdev FC: {summary['stdev']:.3f}x; "
        f"Min FC: {summary['min']:.3f}x; "
        f"Q1 FC: {summary['q1']:.3f}x; "
        f"Q3 FC: {summary['q3']:.3f}x; "
        f"Max FC: {summary['max']:.3f}x"
    )

    if summary["length_weighted_mean"] is not None:
        line += (
            f"; Length-weighted mean FC: {summary['length_weighted_mean']:.3f}x"
            f"; Length-weighted median FC: {summary['length_weighted_median']:.3f}x"
        )

    print(line)
