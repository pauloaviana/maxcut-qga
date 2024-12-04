import networkx as nx
from collections import Counter
import metis

# Get the edges that were cut by METIS
def get_cut_edges(graph, parts):
    partition_dict = {}
    for i, part in enumerate(parts):
        node = list(G.nodes())[i]  # Get the actual node from the NodeView
        partition_dict[node] = part

    # Step 3: Identify cut edges
    cut_edges = []
    for u, v in graph.edges():
        if partition_dict[u] != partition_dict[v]:  # If u and v are in different partitions
            cut_edges.append((u, v))
    return cut_edges

# Partitions graph into N parts according to pre-defined circuit depth
def create_subgraphs(graph, max_circuit_depth):
    part_size = 2
    (edgecuts, subgraphs) = metis.part_graph(graph, nparts=part_size)
    numb_vertices = Counter(subgraphs)
    req_depth_size = max(numb_vertices.values())
    while max_circuit_depth < req_depth_size:
        part_size += 1
        (edgecuts, subgraphs) = metis.part_graph(graph, nparts=part_size)
        numb_vertices = Counter(subgraphs)
        req_depth_size = max(numb_vertices.values())
    numb_part = max(subgraphs)
    parts = []
    for k in range(numb_part+1):
        l = [i for i, p in enumerate(subgraphs) if p == k]
        parts.append(l)
    return edgecuts, parts

def graph_merge(graph, subgraphs):
    final_graph = nx.Graph()

    # Add a vertex for each subgraph
    for i in range(len(subgraphs)):
        final_graph.add_node(i)

    # Add edges between subgraphs based on edges in the original graph
    for i in range(len(subgraphs)):
        for j in range(i+1, len(subgraphs)):
            total_weight = 0
            # Sum edges between subgraph i and subgraph j in the original graph
            for u in subgraphs[i]:
                for v in subgraphs[j]:
                    if graph.has_edge(u, v):
                        total_weight += graph[u][v].get('weight', 1)  # Assuming weighted graph
            if total_weight > 0:
                final_graph.add_edge(i, j, weight=total_weight)
    return final_graph


