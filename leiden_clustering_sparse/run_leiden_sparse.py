import numpy as np
import leidenalg
import igraph as ig
import scipy.io
import scipy.sparse
import sys
import json
import os

def load_sparse_matrix(file_path):
    data = np.loadtxt(file_path, delimiter='\t')
    rows = data[:, 0].astype(int) - 1  # Adjust for MATLAB's 1-based indexing
    cols = data[:, 1].astype(int) - 1  # Adjust for MATLAB's 1-based indexing
    values = data[:, 2]
    n = max(rows.max(), cols.max()) + 1
    adjX = scipy.sparse.coo_matrix((values, (rows, cols)), shape=(n, n))
    return adjX

input_file = 'adjX.txt'
print("Input file: " + input_file)

# Check if input file exists
if not os.path.exists(input_file):
    print(f"Error: Input file {input_file} does not exist.")
    sys.exit(1)

# Load adjacency matrix
try:
    adjX = load_sparse_matrix(input_file)
    adjX = adjX.tocsr()  # Convert to CSR format
    print(f"Adjacency matrix shape: {adjX.shape}")
except Exception as e:
    print(f"Error loading adjacency matrix: {e}")
    sys.exit(1)

# Create a graph from the adjacency matrix
try:
    sources, targets = adjX.nonzero()
    weights = adjX.data
    edges = list(zip(sources, targets))
    
    graph = ig.Graph(edges=edges, directed=False)
    graph.es['weight'] = weights

    # Debug: Print the number of vertices and edges in the graph
    print(f"Graph has {graph.vcount()} vertices and {graph.ecount()} edges")
except Exception as e:
    print(f"Error creating graph: {e}")
    sys.exit(1)

# Define resolution parameter
# Low resolution 0.1 to 1.0
# Moderate resolution from 1.0 to 2.0
# High resolution 2.0 to 4.0+
resolution = 2.0 # Adjust this value as needed
print(f"Resolution parameter: {resolution}")

# Perform Leiden clustering with the resolution parameter
try:
    partition = leidenalg.find_partition(
        graph, 
        leidenalg.RBConfigurationVertexPartition, 
        resolution_parameter=resolution
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