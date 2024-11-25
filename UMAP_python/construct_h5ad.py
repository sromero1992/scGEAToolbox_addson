import numpy as np
import scipy.sparse
import sys
import os
import anndata as ad
import pandas as pd

def load_sparse_matrix(file_path):
    """
    Load a sparse matrix from a tab-delimited file, assumed to contain the following columns:
    [row_index, col_index, value] for the sparse matrix.
    """
    try:
        # Load the sparse matrix data from the file
        data = np.loadtxt(file_path, delimiter='\t')
        # Adjust for MATLAB's 1-based indexing by subtracting 1 from the row/column indices
        rows = data[:, 0].astype(int) - 1  # Adjust for MATLAB's 1-based indexing
        cols = data[:, 1].astype(int) - 1  # Adjust for MATLAB's 1-based indexing
        values = data[:, 2]
        
        # Define shape as (n_genes, n_cells)
        n_genes = rows.max() + 1
        n_cells = cols.max() + 1
        
        # Create the sparse matrix in COO format
        X = scipy.sparse.coo_matrix((values, (rows, cols)), shape=(n_genes, n_cells))
        return X
    except Exception as e:
        print(f"Error loading sparse matrix: {e}")
        sys.exit(1)

# Load the sparse matrix
inputXf = 'sce_X.csv'
# Check if input file exists
if not os.path.exists(inputXf):
    print(f"Error: Input file {inputXf} does not exist.")
    sys.exit(1)

try:
    X = load_sparse_matrix(inputXf)
    X = X.tocsr()  # Convert to CSR format for efficient operations
    print(f"Count matrix shape: {X.shape}")
    # Transpose the matrix (genes as rows, cells as columns)
    X = X.T
    print(f"Transposed matrix shape: {X.shape}")
except Exception as e:
    print(f"Error loading adjacency matrix: {e}")
    sys.exit(1)

# Define file paths for other data
gene_file = "sce_genes.csv"
embedding_file = "sce_embeddings.csv"
cluster_file = "sce_clusters.csv"
celltype_file = "sce_celltypes.csv"
batch_file = "sce_batch.csv"
cellid_file = "sce_cell_ids.csv"

# Function to load CSV files with error handling
def load_csv(file_path, has_header=False):
    if os.path.exists(file_path):
        if has_header:
            return pd.read_csv(file_path)
        else:
            return pd.read_csv(file_path, header=None)
    else:
        print(f"Error: {file_path} does not exist.")
        sys.exit(1)


# Load the other data
genes = load_csv(gene_file)        # Gene names
cell_ids = load_csv(cellid_file)      # Cell IDs
cell_types = load_csv(celltype_file)   # Cell type annotations
batch_ids = load_csv(batch_file)        # Batch IDs
clusters = load_csv(cluster_file)      # Cluster IDs
embeddings = load_csv(embedding_file)  # Cell embeddings

# Convert genes and cell_ids to string if they are not
genes = genes.astype(str).values   # Ensure genes are strings
cell_ids = cell_ids.astype(str).values  # Ensure cell_ids are strings
cell_types = cell_types.astype(str).values  # Cell type annotations
batch_ids = batch_ids.astype(str).values
clusters = clusters.values
embeddings = embeddings.values

# Flatten the arrays to 1D
cell_types = cell_types.flatten()
batch_ids = batch_ids.flatten()
cell_ids = cell_ids.flatten()
clusters = clusters.flatten()

# Create the DataFrame
cell_data = pd.DataFrame({
    "cell_type": cell_types,
    "batch_id": batch_ids,
    "cell_id": cell_ids,
    "cluster": clusters
})

# Print the first few rows to verify
#print(cell_data.head())

# Create the AnnData object
adata = ad.AnnData(
    X=X,
    obs= cell_data,
    var=pd.DataFrame(index=genes),  # Gene names (ensure these are strings)
    obsm={"X_umap": embeddings}  # Embeddings
)

# Save as an h5ad file
output_file = "sce_data.h5ad"
try:
    adata.write(output_file)
    print(f"AnnData saved to {output_file}")
except Exception as e:
    print(f"Error saving AnnData: {e}")
    sys.exit(1)

# Cleanup: Delete the temporary CSV files
files_to_delete = [
    'sce_X.csv', 'sce_genes.csv', 'sce_embeddings.csv', 'sce_clusters.csv',
    'sce_celltypes.csv', 'sce_batch.csv', 'sce_cell_ids.csv'
]
for file in files_to_delete:
    if os.path.exists(file):
        os.remove(file)
