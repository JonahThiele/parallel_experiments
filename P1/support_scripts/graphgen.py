import networkx as nx
import random

def generate_nauty_format(graph, filename):
    """Generate a graph file in the format your code expects"""
    n = graph.number_of_nodes()
    
    with open(filename, 'w') as f:
        # Write number of vertices
        f.write(f"{n}\n")
        
        # Write edges as vertex1,vertex2
        for u, v in graph.edges():
            f.write(f"{u},{v}\n")

# Generate different types of graphs

# Random graph (Erdős-Rényi)
n = 1500
p = 0.85  # edge probability
g = nx.erdos_renyi_graph(n, p)
generate_nauty_format(g, "random_1500_85.txt")

# # Scale-free graph (Barabási-Albert)
# n = 10000
# m = 3  # edges to attach from new node
# g = nx.barabasi_albert_graph(n, m)
# generate_nauty_format(g, "scalefree_10000.txt")

# # Small-world graph (Watts-Strogatz)
# n = 10000
# k = 6  # each node connected to k nearest neighbors
# p = 0.1  # rewiring probability
# g = nx.watts_strogatz_graph(n, k, p)
# generate_nauty_format(g, "smallworld_10000.txt")

# # Regular graph
# n = 10000
# d = 4  # degree
# g = nx.random_regular_graph(d, n)
# generate_nauty_format(g, "regular_10000.txt")

# # Grid graph
# side = 100  # 100x100 = 10000 nodes
# g = nx.grid_2d_graph(side, side)
# # Relabel nodes to integers
# g = nx.convert_node_labels_to_integers(g)
# generate_nauty_format(g, "grid_10000.txt")

print(f"Generated graphs with {n} nodes")