#!/usr/bin/env python

"""
gfa_stats - Reports GFA assembly graph coverage, length, and connectivity statistics.
"""

import argparse
import math
import sys

import gfapy

from utils.log import print_message
import utils.stats


def _format_cov(value):
    if value is None:
        return "NA"
    return f"{value:.3f}x"


def _extract_node_coverages(graph):
    coverages = []

    for segment in graph.segments:
        try:
            coverage = utils.stats.get_coverage(segment)
        except Exception:
            continue

        if coverage is None:
            continue

        coverages.append(float(coverage))

    return coverages


def _fixed_0_to_max_histogram(coverages, max_cov):
    counts = [0] * (max_cov + 1)
    overflow = 0

    for coverage in coverages:
        if coverage < 0:
            continue

        if coverage > max_cov:
            overflow += 1
            continue

        counts[int(coverage)] += 1

    histogram = []
    for coverage_value, count in enumerate(counts):
        histogram.append((float(coverage_value), float(coverage_value), count))

    histogram.append((float(max_cov + 1), math.inf, overflow))
    return histogram


def _print_ascii_histogram(coverages, max_cov=100, chart_width=50):
    histogram = _fixed_0_to_max_histogram(coverages, max_cov)

    print()
    print("Coverage distribution (FC tag)")
    print("-----------------------------")

    if len(histogram) == 0:
        print("No FC-tagged nodes found, histogram unavailable.")
        return

    max_count = max(count for _, _, count in histogram)
    max_count = max(max_count, 1)

    for low, high, count in histogram:
        bar_len = int(round((count / max_count) * chart_width))
        bar = "#" * bar_len

        if math.isinf(high):
            label = f"> {int(low) - 1:>6d}"
        elif math.isclose(low, high):
            label = f"{low:>8.3f}"
        else:
            label = f"{low:>8.3f} - {high:>8.3f}"

        print(f"{label} | {bar} ({count:,})")


def main():
    parser = argparse.ArgumentParser(
        description="Reports GFA assembly graph coverage, length, and connectivity statistics."
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input GFA file")
    parser.add_argument(
        "--max-cov",
        type=int,
        default=100,
        help="Maximum FC value to show as individual bins from 0..max-cov (default: 100); values above this go to an overflow bin",
    )
    parser.add_argument(
        "--chart-width",
        type=int,
        default=50,
        help="Maximum character width for histogram bars (default: 50)",
    )

    parser.add_argument("--vlevel", type=int, help="Validation level for GFA parsing (0, 1, 2, or 3, default: 0)", default=0, choices=[0, 1, 2, 3])

    args = parser.parse_args()

    print_message("Importing GFA file:", args.input)

    try:
        graph = gfapy.Gfa.from_file(args.input, vlevel=args.vlevel)
    except Exception as error:
        print_message("ERROR, failed to parse GFA file:", str(error))
        sys.exit(1)

    print_message("Imported:", args.input)

    summary = utils.stats.coverage_summary(graph)

    print()
    print("Summary")
    print("-------------------------")
    print(f"Nodes in graph: {len(graph.segments):,}")
    print(f"Nodes with FC: {summary['nodes_with_fc']:,}")
    print(f"Nodes missing FC: {summary['nodes_missing_fc']:,}")
    print(f"Mean FC: {_format_cov(summary['mean'])}")
    print(f"Median FC: {_format_cov(summary['median'])}")
    print(f"Stdev FC: {_format_cov(summary['stdev'])}")
    print(f"Min FC: {_format_cov(summary['min'])}")
    print(f"Q1 FC: {_format_cov(summary['q1'])}")
    print(f"Q3 FC: {_format_cov(summary['q3'])}")
    print(f"Max FC: {_format_cov(summary['max'])}")
    print(f"Length-weighted mean FC: {_format_cov(summary['length_weighted_mean'])}")
    print(f"Length-weighted median FC: {_format_cov(summary['length_weighted_median'])}")

    cov_percentiles = utils.stats.coverage_percentiles(graph)

    print()
    print("Coverage percentiles (FC tag)")
    print("-----------------------------")

    if len(cov_percentiles) == 0:
        print("No FC-tagged nodes found.")
    else:
        percentiles_order = ["p1", "p5", "p10", "p25", "p50", "p75", "p90", "p95", "p99"]
        for key in percentiles_order:
            if key in cov_percentiles:
                print(f"  {key}: {_format_cov(cov_percentiles[key])}")

    length_stats = utils.stats.length_percentiles_and_n_stats(graph)

    print()
    print("Length percentiles and N-statistics")
    print("-----------------------------------")

    if len(length_stats) == 0:
        print("No nodes found.")
    else:
        percentiles_order = ["p1", "p5", "p10", "p25", "p50", "p75", "p90", "p95", "p99"]
        n_order = ["N10", "N25", "N50", "N75", "N90"]

        print("Percentiles:")
        for key in percentiles_order:
            if key in length_stats:
                print(f"  {key}: {length_stats[key]:,} bp")

        print("N-statistics:")
        for key in n_order:
            if key in length_stats:
                print(f"  {key}: {length_stats[key]:,} bp")

    connectivity_counts = utils.stats.connectivity_summary(graph)

    print()
    print("Connectivity summary")
    print("--------------------")

    if len(connectivity_counts) == 0:
        print("No nodes found.")
    else:
        for neighbour_count in sorted(connectivity_counts):
            node_count = connectivity_counts[neighbour_count]
            print(f"{node_count:,} nodes have {neighbour_count:,} neighbours")

    component_stats = utils.stats.component_summary(graph)

    print()
    print("Component summary")
    print("-----------------")
    print(f"Connected components: {component_stats['n_components']:,}")

    if component_stats["n_components"] == 0:
        print("No components found.")
    else:
        node_summary = component_stats["node_size_summary"]
        bp_summary = component_stats["bp_size_summary"]

        print("Component size by nodes:")
        print(
            f"  min={node_summary['min']:,}; p10={node_summary['p10']:,}; p50={node_summary['p50']:,}; "
            f"p90={node_summary['p90']:,}; p99={node_summary['p99']:,}; max={node_summary['max']:,}; "
            f"mean={node_summary['mean']:.2f}; median={node_summary['median']:.2f}"
        )

        print("Component size by bp:")
        print(
            f"  min={bp_summary['min']:,} bp; p10={bp_summary['p10']:,} bp; p50={bp_summary['p50']:,} bp; "
            f"p90={bp_summary['p90']:,} bp; p99={bp_summary['p99']:,} bp; max={bp_summary['max']:,} bp; "
            f"mean={bp_summary['mean']:.2f} bp; median={bp_summary['median']:.2f} bp"
        )

        print(
            f"Largest component: {component_stats['largest_component_nodes']:,} nodes "
            f"({component_stats['largest_component_nodes_pct']:.2f}% of nodes); "
            f"{component_stats['largest_component_bp']:,} bp "
            f"({component_stats['largest_component_bp_pct']:.2f}% of bp)"
        )

        print("Top components (by bp):")
        for index, (nodes, bp) in enumerate(component_stats["top_components"][:10], start=1):
            print(f"  {index}. {nodes:,} nodes; {bp:,} bp")

    node_coverages = _extract_node_coverages(graph)
    _print_ascii_histogram(
        node_coverages,
        max_cov=max(0, args.max_cov),
        chart_width=max(1, args.chart_width),
    )


if __name__ == "__main__":
    main()
