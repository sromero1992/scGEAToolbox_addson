import numpy as np
import leidenalg
import igraph as ig
import sys
import json
import os

input_file = 'adjX.txt'
print("Input file: " + input_file)

# Check if input file exists and is not empty
if not os.path.exists(input_file):
    print(f"Error: Input file {input_file} does not exist.")
    sys.exit(1)

if os.stat(input_file).st_size == 0:
    print(f"Error: Input file {input_file} is empty.")
    sys.exit(1)

# Load adjacency matrix
try:
    adjX = np.loadtxt(input_file, delimiter='\t')
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

# Define resolution parameter
# Low resolution 0.1 to 1.0
# Moderate resolution from 1.0 to 2.0
# High resolution 2.0 to 4.0+
resolution = 1.0  # Adjust this value as needed
print(f"Resolution parameter: {resolution}")

# Perform Leiden clustering with the resolution parameter
try:
    partition = leidenalg.find_partition(
        graph, 
        leidenalg.RBConfigurationVertexPartition, 
        resolution_parameter = resolution
    )
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