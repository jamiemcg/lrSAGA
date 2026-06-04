#!/usr/bin/env python

"""
lrSAGA - long read Single Amplified Genome Assembly
https://github.com/jamiemcg/lrSAGA/

Version 0.1
Jamie McGowan, 2026
"""

import argparse
import datetime
import os
import sys

import gfapy

from utils.log import print_message
import utils.stats
import utils.simplification

def main():
    parser = argparse.ArgumentParser(description = "lrSAGA - simplify and process a GFA assembly graph to remove chimeras",
                                     add_help=False)

    required = parser.add_argument_group("Required arguments")

    required.add_argument("-i", "--input", type=str, help="Input GFA file (required)", required=True)
    required.add_argument("-o", "--output", type=str, help="Output GFA file (required)", required=True)

    optional = parser.add_argument_group("Optional arguments")

    optional.add_argument("-h", "--help", action="help", help="Show this help message and exit")

    optional.add_argument("-v", "--verbose", help="Write verbose log file (appends to <output>.verbose.log)", action="store_true")

    optional.add_argument("--save-intermediate", help="Saves intermediate GFA files to <output>.intermediate.gfa", action="store_true")

    optional.add_argument("--tip-len", type=int, help="Remove tips shorter than --tip-len (requires --tip-cov)", required="--tip-cov" in sys.argv)
    optional.add_argument("--tip-cov", type=float, help="Remove tips with coverage less than --tip-cov (requires --tip-len)", required="--tip-len" in sys.argv)

    optional.add_argument("--pop-bubbles", help="Pop simple bubbles", action="store_true")
    optional.add_argument("--max-bubble-len", type=int, help="If defined, don't discard bubble longer than this")

    optional.add_argument("--inverted-chimera-len", type=int, help="Remove inverted chimeras shorter than --inverted-chimera-len (requires --inverted-chimera-cov)", required="--inverted-chimera-cov" in sys.argv)
    optional.add_argument("--inverted-chimera-cov", type=float, help="Remove inverted chimeras with coverage less than --inverted-chimera-cov (requires --inverted-chimera-len)", required="--inverted-chimera-len" in sys.argv)

    optional.add_argument("--transitive-node-len", type=int, help="Remove transitive nodes shorter than --transitive-node-len (requires --transitive-node-cov)", required="--transitive-node-cov" in sys.argv)
    optional.add_argument("--transitive-node-cov", type=float, help="Remove transitive nodes with coverage less than --transitive-node-cov (requires --transitive-node-len)", required="--transitive-node-len" in sys.argv)

    optional.add_argument("--by-property-len", type=int, help="Remove nodes with length shorter than --by-property-len (requires --by-property-cov)", required="--by-property-cov" in sys.argv)
    optional.add_argument("--by-property-cov", type=float, help="Remove nodes with coverage less than --by-property-cov (requires --by-property-len)", required="--by-property-len" in sys.argv)

    optional.add_argument("--isolated-len", type=int, help="Remove isolated nodes (no edges) with length shorter than --isolated-len (requires --isolated-cov)", required="--isolated-cov" in sys.argv)
    optional.add_argument("--isolated-cov", type=float, help="Remove isolated nodes (no edges) with coverage less than --isolated-cov (requires --isolated-len)", required="--isolated-len" in sys.argv)

    optional.add_argument("--remove-list", type=str, help="Provide a plain text file containing a list of nodes to remove")

    optional.add_argument("--vlevel", type=int, help="Validation level for GFA parsing (0, 1, 2, or 3, default: 0)", default=0, choices=[0, 1, 2, 3])

    args = parser.parse_args()

    print_message("Starting lrSAGA")
    print_message("User provided arguments:", sys.argv)
    print_message("Parsed arguments:", vars(args))
    print()

    if os.path.exists(args.output):
        print_message("ERROR, output file already exists:", args.output)
        sys.exit(1)

    if args.input == args.output:
        print_message("ERROR, you specified the same parameter for input and output files")
        sys.exit(1)

    print_message("Importing GFA file:", args.input)
    graph = gfapy.Gfa.from_file(args.input, vlevel=args.vlevel)
    print_message("Imported:", args.input)

    if args.verbose:
        print_message("Writing verbose log to:", args.output + ".verbose.log")
        verbose_log_file = open(args.output + ".verbose.log", "a")
        timestamp = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        verbose_log_file.write("#[" + timestamp + "] " +str(vars(args)) + "\n")

    print()
    utils.stats.segment_stats(graph)
    print()

    round_num = 1

    while True:
        user_specified_list = set()
        tips = set()
        inverted_chimeras = set()
        bubbles = set()
        transitive = set()
        by_property = set()
        isolated = set()

        print_message("ROUND", str(round_num))

        if args.verbose:
            verbose_log_file.write("#round " + str(round_num) + "\n")

        # Remove the user provided list of nodes. Only do this on the first round
        if round_num == 1 and args.remove_list:
            print_message("Removing nodes specified by user in file: " + args.remove_list)
            if not os.path.exists(args.remove_list):
                print("Error", args.remove_list, "does not exist")
                sys.exit()
            
            with open(args.remove_list, "r") as f:
                for line in f:
                    node_id = line.strip()

                    if node_id in graph.segment_names:
                        user_specified_list.add(node_id)
                    else:
                        print_message("Warning: '" + node_id + "' not found in graph")

            print_message(str(len(user_specified_list)) + " nodes specified by user identified")

        if args.tip_len:
            print_message("Identifying tips <= " + str(args.tip_len) + " bp and <= " + str(args.tip_cov) + "x")
            tips = utils.simplification.identify_tips(graph, length_cutoff=args.tip_len, coverage_cutoff=args.tip_cov)
            print_message(str(len(tips)) + " tips identified")

            if args.verbose:
                for node in tips:
                    verbose_log_file.write("\t".join([node.name, str(len(node.sequence)), str(utils.stats.get_coverage(node)), "tip"]) + "\n")

        if args.inverted_chimera_len:
            print_message("Identifying inverted chimeras <= " + str(args.inverted_chimera_len) + " bp and <= "  + str(args.inverted_chimera_cov) + "x")
            inverted_chimeras = utils.simplification.identify_putative_chimeras(graph, length_cutoff=args.inverted_chimera_len, coverage_cutoff=args.inverted_chimera_cov)
            print_message(str(len(inverted_chimeras)) + " inverted chimeras identified")

            if args.verbose:
                for node in inverted_chimeras:
                    verbose_log_file.write("\t".join([node.name, str(len(node.sequence)), str(utils.stats.get_coverage(node)), "inverted_chimera"]) + "\n")

        if args.pop_bubbles:
            print_message("Identifying simple bubbles")
            if args.max_bubble_len:
                # print("defined max bubble len")
                bubbles = utils.simplification.identify_simple_bubbles(graph, args.max_bubble_len)
            else:
                bubbles = utils.simplification.identify_simple_bubbles(graph)
            print_message(str(len(bubbles)) + " bubbles identified")

            if args.verbose:
                for node in bubbles:
                    verbose_log_file.write("\t".join([node.name, str(len(node.sequence)), str(utils.stats.get_coverage(node)), "bubble"]) + "\n")

        if args.by_property_len:
            print_message("Identifying nodes <= " + str(args.by_property_len) + " bp and <= " + str(args.by_property_cov) + "x")
            by_property = utils.simplification.identify_nodes_by_length_and_coverage(graph, length_cutoff=args.by_property_len, coverage_cutoff=args.by_property_cov)
            print_message(str(len(by_property)) + " such nodes identified")

            if args.verbose:
                for node in by_property:
                    verbose_log_file.write("\t".join([node.name, str(len(node.sequence)), str(utils.stats.get_coverage(node)), "by_property"]) + "\n")

        if args.isolated_len:
            print_message("Identifying isolated nodes <= " + str(args.isolated_len) + " bp and <= " + str(args.isolated_cov) + "x")
            isolated = utils.simplification.identify_nodes_by_property(graph, length_cutoff=args.isolated_len, coverage_cutoff=args.isolated_cov, num_input_edges=0, num_output_edges=0)
            print_message(str(len(isolated)), "isolated nodes identified")

            if args.verbose:
                for node in isolated:
                    verbose_log_file.write("\t".join([node.name, str(len(node.sequence)), str(utils.stats.get_coverage(node)), "isolated"]) + "\n")

        if args.transitive_node_len:
            print_message("Identifying transitive nodes <= " + str(args.transitive_node_len) + " bp and <= " + str(args.transitive_node_cov) + "x")
            transitive = utils.simplification.identify_transitive_nodes(graph, length_cutoff=args.transitive_node_len, coverage_cutoff=args.transitive_node_cov)
            print_message(str(len(transitive)) + " transitive nodes identified")

            if args.verbose:
                for node in transitive:
                    verbose_log_file.write("\t".join([node.name, str(len(node.sequence)), str(utils.stats.get_coverage(node)), "transitive"]) + "\n")

        # Join all sets of nodes to be removed into one set
        to_remove = tips | inverted_chimeras | bubbles | by_property | isolated | transitive | user_specified_list

        if len(to_remove) == 0:
            print_message("Identified 0 nodes to remove. Exit.")

            break

        print_message("Removing", len(to_remove), "nodes")
        print()
        utils.simplification.remove_from_graph(graph, to_remove)

        utils.stats.segment_stats(graph)

        print()
        print_message("Identifying linear paths to merge")

        merged_linear_nodes = utils.simplification.manual_merge_linear_paths(graph)

        if args.save_intermediate:
            base, ext = os.path.splitext(args.output)
            intermediate_filename = base + "_intermediate" + ext

            print_message("Writing intermediate graph to file:", intermediate_filename)
            graph.to_file(intermediate_filename)

        if args.verbose:
            verbose_log_file.write("#merged_linear_nodes: " + ", ".join(merged_linear_nodes) + "\n")
            verbose_log_file.flush()

        print()

        utils.stats.segment_stats(graph)

        print()

        round_num += 1

    print_message("Writing graph to file:", args.output)
    graph.to_file(args.output)

    if args.verbose:
        verbose_log_file.write("#end\n")
        verbose_log_file.close()

if __name__ == "__main__":
    main()
