import random

import networkx as nx

from utils.stats import get_coverage
from utils.log import print_message

# There is an issue with Gfapy when using merge_linear_paths
# This function merges each linear path individually
def manual_merge_linear_paths(graph):
    linear_paths = graph.linear_paths()

    removed_segments = set()

    for path in linear_paths:
        total_FC = 0
        path_segments = []

        for segment in path:
            removed_segments.add(segment.name)
            path_segments.append(segment.name)
            total_FC += graph.segment(segment.name).FC

            # Disconnect self-links
            if segment.name in graph.segment(segment.name).neighbours:
                for edge in graph.segment(segment.name).edges:
                    if segment.name == edge.from_name and segment.name == edge.to_name:
                        try:
                            edge.disconnect()
                        except:
                            pass

        new_name = "_".join(path_segments)

        graph.merge_linear_path(path, merged_name=new_name)
        graph.segment(new_name).set("FC", total_FC)

    for edge in graph.edges:
        if edge.from_name in removed_segments or edge.to_name in removed_segments:
            graph.rm(edge)

    for segment in graph.segments:
        if segment.length is None:
            # print(f"[{segment.name} length is None]")
            # if segment.virtual:
                # print(f"[{segment.name} is virtual]")
            graph.rm(segment)
            removed_segments.add(segment.name)

    print_message("Merged " + str(len(linear_paths)) + " linear paths including " + str(len(removed_segments)) + " nodes")
    return removed_segments

def remove_from_graph(graph, to_remove):
    for segment in to_remove:
        graph.rm(segment)

    # Also check for links that may not have been deleted
    for link in graph.edges:
        if link.from_segment in to_remove or link.to_segment in to_remove:
            graph.rm(link)

# Tips are segments that only have a neighbour on one side
# Also removing tips that form loops
def identify_tips(graph, length_cutoff, coverage_cutoff):
    tips = set()

    for segment in graph.segments:
        if len(set(segment.neighbours)) == 1:
            if get_coverage(segment) <= coverage_cutoff and segment.length <= length_cutoff:
                if len(segment.neighbours_L) == 0 and len(segment.neighbours_R) == 1:
                    tips.add(segment)
                elif len(segment.neighbours_R) == 0 and len(segment.neighbours_L) == 1:
                    tips.add(segment)
                elif len(segment.neighbours_L) > 0 and len(
                        segment.neighbours_R) > 0 and segment.neighbours_L == segment.neighbours_R:
                    tips.add(segment)
                elif len(set(segment.neighbours)) == 1 and len(segment.neighbours) > 1:
                    tips.add(segment)
    return tips

# NOTE this removes more than just bubbles (e.g., some chimeras) (more aggressive)
# Checks for subset of neighbour sharing instead of exact neighbour sharing
# [Consider renaming this function]
def identify_simple_bubbles(graph, max_bubble_len=None):
    bubbles = set()

    # Store all neighbours in a dict
    neighbours = {}
    for segment in graph.segments:
        neighbours[segment.name] = set(_.name for _ in segment.neighbours)

    # Visit each segment
    for segment in graph.segments:
        # Visit each neighbour
        for neighbour1 in neighbours[segment.name]:
            # Visit each neighbour of that neighbour
            for neighbour2 in neighbours[neighbour1]:
                if neighbour2 != segment.name and neighbours[neighbour2].issubset(neighbours[segment.name]):
                    if (len(segment.neighbours_L) > 0 and len(segment.neighbours_R) > 0 and
                            len(graph.segment(neighbour2).neighbours_L) > 0 and
                            len(graph.segment(neighbour2).neighbours_R) > 0):

                        # Decide which segment to remove based on coverage
                        x = segment
                        y = graph.segment(neighbour2)

                        # if either of these nodes is longer than max_bubble_len, don't discard
                        if max_bubble_len and (x.length > max_bubble_len or y.length > max_bubble_len):
                            # print_message("Skipping, max_bubble_len defined", max_bubble_len, "X:", x.name, x.length, "; Y:", y.name, y.length)
                            continue

                        # Ensure this isn't a bubble where we've already flagged one of the member segments for deletion
                        if x not in bubbles and y not in bubbles:
                            x_coverage = get_coverage(x)
                            y_coverage = get_coverage(y)

                            # Remove the segment with lower coverage. If both segments have the same coverage, remove the lexicographically smaller name
                            if x_coverage == y_coverage:
                                bubbles.add(min((x, y), key=lambda s: s.name))
                            elif x_coverage < y_coverage:
                                bubbles.add(x)
                            else:
                                bubbles.add(y)
    return bubbles

def identify_putative_chimeras(graph, length_cutoff, coverage_cutoff):
    chimeras = set()
    neighbours = {}

    # First store all neighbours in a dict
    for segment in graph.segments:
        neighbours[segment.name] = set(_.name for _ in segment.neighbours)

    for segment in neighbours:
        if len(neighbours[segment]) == 2:
            x = list(neighbours[segment])[0]
            y = list(neighbours[segment])[1]

            if x in neighbours[y] and y in neighbours[x]:
                if get_coverage(graph.segment(segment)) <= coverage_cutoff and graph.segment(
                        segment).length <= length_cutoff:
                    chimeras.add(graph.segment(segment))

    return chimeras

def identify_transitive_nodes(graph, length_cutoff, coverage_cutoff):
    segments = set()
    n_segments = len(graph.segments)

    # Precompute neighbors
    neighbours = {seg.name: set(_.name for _ in seg.neighbours) for seg in graph.segments}
    neighbours_L = {seg.name: set(_.name for _ in seg.neighbours_L) for seg in graph.segments}
    neighbours_R = {seg.name: set(_.name for _ in seg.neighbours_R) for seg in graph.segments}

    # Build NetworkX graph once
    nx_graph = nx.Graph()
    for edge in graph.edges:
        nx_graph.add_edge(edge.from_name, edge.to_name)

    for segment_name, segment_neighbors in neighbours.items():
        if get_coverage(graph.segment(segment_name)) <= coverage_cutoff:
            if len(segment_neighbors) == 2:
                if len(neighbours_L[segment_name]) == 1 and len(neighbours_R[segment_name]) == 1:
                    neighbour_L = list(neighbours_L[segment_name])[0]
                    neighbour_R = list(neighbours_R[segment_name])[0]

                    # Temporarily remove the edges/nodes from the graph
                    nx_graph.remove_edge(segment_name, neighbour_L)
                    nx_graph.remove_edge(segment_name, neighbour_R)
                    # nx_graph.remove_node(segment_name)
                    if nx.has_path(nx_graph, neighbour_L, neighbour_R):
                        # transitive node, remove
                        nx_graph.remove_node(segment_name)
                        segments.add(graph.segment(segment_name))
                    else:
                        # node segment wasn't transitive, restore node in graph
                        # nx_graph.add_node(segment_name)
                        nx_graph.add_edge(segment_name, neighbour_L)
                        nx_graph.add_edge(segment_name, neighbour_R)

    return segments

def identify_nodes_by_property(graph, length_cutoff, coverage_cutoff, num_input_edges, num_output_edges):
    segments = set()

    for segment in graph.segments:
        if segment.length <= length_cutoff:
            if get_coverage(segment) <= coverage_cutoff:
                if len(segment.neighbours_L) == num_input_edges:
                    if len(segment.neighbours_R) == num_output_edges:
                        segments.add(segment)

    return segments

def identify_nodes_by_length_and_coverage(graph, length_cutoff, coverage_cutoff):
    segments = set()

    for segment in graph.segments:
        if segment.length <= length_cutoff:
            if get_coverage(segment) <= coverage_cutoff:
                segments.add(segment)

    return segments
