import networkx as nx
import itertools
import metis
from collections import Counter

class MaxCut:
    def __init__(self, graph):
        self.graph = graph

    def evaluate(self, binary_list):
        # Ensure the binary list matches the number of nodes in G
        if len(binary_list) != self.graph.number_of_nodes():
            raise ValueError("Binary list length must match the number of nodes in the graph.")
    
        # Create two sets based on the binary list
        set1 = {node for node, value in zip(self.graph.nodes(), binary_list) if value == 0}
        set2 = {node for node, value in zip(self.graph.nodes(), binary_list) if value == 1}
    
        # Calculate the weight of the edges between set1 and set2
        maxcut_fitness = sum(self.graph[u][v]['weight'] for u in set1 for v in set2 if self.graph.has_edge(u, v))
        return maxcut_fitness
    def lower_bound(self):
        # Sum all the edge weights in the graph
        total_weight = sum(data.get('weight', 1) for u, v, data in self.graph.edges(data=True))
        
        # Return half of the total weight as the lower bound
        return total_weight / 2


def max_cut_brute_force(G):
    # Get all nodes in the graph
    nodes = list(G.nodes())
    num_nodes = len(nodes)
    
    max_cut_value = 0
    best_partition = None

    # Generate all possible partitions by iterating through all combinations of node subsets
    for i in range(1, num_nodes//2 + 1):
        for subset in itertools.combinations(nodes, i):
            set1 = set(subset)
            set2 = set(nodes) - set1
            
            # Calculate the cut value for this partition
            cut_value = sum(1 for u in set1 for v in set2 if G.has_edge(u, v))
            
            # Update max_cut_value and best_partition if this cut is better
            if cut_value > max_cut_value:
                max_cut_value = cut_value
                best_partition = (set1, set2)

    return max_cut_value, best_partition

# Auxilliary functions

def to_binary(value, num_bits, symm=False):
    """
    value : int
    num_bits : int

    #Converts integer into binary inside a list, with each element being a digit of the number

    symm=True returns the symmetrical value, i.e. not(a) for each 'a' in list
    """
    binary_representation = []
    for _ in range(num_bits):
        binary_representation.append(value % 2)
        value //= 2
    if symm == True:
        binary_representation.reverse()
    return binary_representation

def create_subgraphs(graph, max_circuit_depth):
    part_size = 2
    (edgecuts, subgraphs) = metis.part_graph(graph, nparts=part_size)
    numb_vertices = Counter(subgraphs)
    req_depth_size = max(numb_vertices.values())
    
    # Adjust partition size to meet circuit depth constraint
    while max_circuit_depth < req_depth_size:
        part_size += 1
        (edgecuts, subgraphs) = metis.part_graph(graph, nparts=part_size)
        numb_vertices = Counter(subgraphs)
        req_depth_size = max(numb_vertices.values())
    
    # Create subgraphs as networkx graphs
    partitioned_subgraphs = []
    for part_id in range(part_size):
        nodes_in_part = [node for node, partition in enumerate(subgraphs) if partition == part_id]
        subgraph = graph.subgraph(nodes_in_part).copy()  # Create a subgraph and ensure it's a separate object
        partitioned_subgraphs.append(subgraph)
    
    return edgecuts, partitioned_subgraphs

def graph_merge(graph, subgraphs, bipartitions):
    """
    Merges subgraphs into a single graph, reducing each subgraph to two vertices
    based on bipartitions and preserving all edges, with edge weights passed in 
    the bipartition list.

    Parameters:
    - graph: The original NetworkX graph.
    - subgraphs: A list of lists, each containing the nodes of a subgraph.
    - bipartitions: A list of lists, each containing a weight and a binary string representing
                     the bipartition of the subgraph.

    Returns:
    - A new merged NetworkX graph.
    """
    final_graph = nx.Graph()
    contraction_list = []
    partition_list = []
    original_graph= graph
    

    # Reduce each subgraph into two vertices
    for i, (subgraph_nodes, bipartition_info) in enumerate(zip(subgraphs, bipartitions)):
        cut_weight, bipartition = bipartition_info
        cut_weight = evaluate(subgraphs[i], bipartition)
        
        # Ensure subgraph is a subgraph of graph
        subgraph_nodes = set(subgraphs[i].nodes)
        if not subgraph_nodes.issubset(graph.nodes):
            raise ValueError("subgraph must be a subgraph of graph")
    
        # Create final_graph with two initial vertices and an edge
        contraction = nx.Graph()
        contraction_vertices = [list(subgraphs[i].nodes)[0], list(subgraphs[i].nodes)[-1]]
        contraction.add_edge(contraction_vertices[0], contraction_vertices[-1], weight=cut_weight)  # Original vertices in final_graph
        #if graph.has_edge(contraction_vertices[0], contraction_vertices[-1]):  # Check if the edge exists
        #    graph.remove_edge(contraction_vertices[0], contraction_vertices[-1])  # Remove the edge
    
        # Map subgraph nodes to their position in the bipartitions string
        list0, list1 = [],[]
        part = list(bipartitions[i][1])
        subgraph_nodes_ordered = list(subgraph_nodes)
        if len(part) != len(subgraph_nodes_ordered):
            raise ValueError("bipartitions length must match the number of vertices in subgraph")
        for i, vertex in enumerate(subgraph_nodes_ordered):
            if int(part[i]) == 0:
                list0.append(vertex)
            else:
                list1.append(vertex)
        contraction_list.append(contraction)
        partition_list.append((list0,list1))
        final_graph = nx.compose(final_graph, contraction)

    for edge in graph.edges:
        indexes = [None,None]
        for idx, subgraph in enumerate(subgraphs):
        # Check if vertex0 or vertex1 belong to the current subgraph
            if edge[0] in subgraph.nodes:
                indexes[0] = idx
            if edge[1] in subgraph.nodes:
                indexes[1] = idx
        for i in range(2):
            if edge[0] in partition_list[indexes[0]][i]:
                v0_belongs_to = i
            if edge[1] in partition_list[indexes[1]][i]:
                v1_belongs_to = i
                
        vertex0 = list(contraction_list[indexes[0]].nodes)[v0_belongs_to]
        vertex1 = list(contraction_list[indexes[1]].nodes)[v1_belongs_to]

        if indexes[0] != indexes[1]:
            if final_graph.has_edge(vertex0, vertex1) and original_graph.has_edge(vertex0, vertex1):
                if 'weight' not in final_graph[vertex0][vertex1]:
                    final_graph[vertex0][vertex1]['weight'] = original_graph[vertex0][vertex1].get('weight', None)
                else:
                    final_graph[vertex0][vertex1]['weight'] += original_graph[vertex0][vertex1].get('weight', None)
            else:
                if original_graph.has_edge(vertex0, vertex1):
                    edge_data = original_graph[vertex0][vertex1]
                    ws = original_graph[vertex0][vertex1].get('weight', None)
                    final_graph.add_edge(vertex0, vertex1, weight = ws)
                else: final_graph.add_edge(vertex0, vertex1, weight = 1)

    return final_graph

def evaluate(graph, part):
    binary_representation = list(part)
    binary_list = [int(bit) for bit in binary_representation]
    
    # Create two sets based on the binary representation
    set1 = {node for node, value in zip(graph.nodes(), binary_list) if value == 0}
    set2 = {node for node, value in zip(graph.nodes(), binary_list) if value == 1}
    
    # Calculate the Max-Cut fitness by counting edges between set1 and set2
    maxcut_fitness = sum(graph[u][v].get('weight', None) for u in set1 for v in set2 if graph.has_edge(u, v))
    
    return maxcut_fitness
