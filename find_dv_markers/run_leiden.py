import numpy as np
import leidenalg
import igraph as ig
import sys
import json
import os

input_file = 'adjX.txt'
print("Input file: " + input_file)

# Check if input file exists
if not os.path.exists(input_file):
    print(f"Error: Input file {input_file} does not exist.")
    sys.exit(1)

# Load adjacency matrix
try:
    adjX = np.loadtxt(input_file)
    print(f"Adjacency matrix shape: {adjX.shape}")
except Exception as e:
    print(f"Error loading adjacency matrix: {e}")
    sys.exit(1)

# Create a graph from the adjacency matrix
try:
    graph = ig.Graph.Adjacency((adjX > 0).tolist(), mode=ig.ADJ_UNDIRECTED)
    graph.es['weight'] = adjX[adjX.nonzero()]

    # Debug: Print the number of vertices and edges in the graph
    print(f"Graph has {graph.vcount()} vertices and {graph.ecount()} edges")

except Exception as e:
    print(f"Error creating graph: {e}")
    sys.exit(1)

# Perform Leiden clustering
try:
    partition = leidenalg.find_partition(graph, leidenalg.ModularityVertexPartition)
except Exception as e:
    print(f"Error in Leiden clustering: {e}")
    sys.exit(1)

# Get the clustering result
clusters = partition.membership

# Save the result to a file
output_file = 'clusters.json'
try:
    with open(output_file, 'w') as f:
        json.dump(clusters, f)
    print(f"Leiden clustering completed. Results saved to {output_file}")
except Exception as e:
    print(f"Error saving results to file: {e}")
    sys.exit(1)