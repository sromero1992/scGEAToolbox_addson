import numpy as np
import scipy.sparse
import sys
import os
import anndata as ad
import pandas as pd
import glob  # <--- IMPORT GLOB TO FIND FILES

def load_sparse_matrix(file_path):
    """
    Load a sparse matrix from a tab-delimited file using pandas for robust I/O.
    The file is assumed to contain [row_index, col_index, value].
    """
    try:
        # Use pandas.read_csv for faster and more reliable file handling
        data = pd.read_csv(file_path, delimiter='\t', header=None).values

        if data.ndim == 1:
            data = data.reshape(1, -1)
        
        rows = data[:, 0].astype(int) - 1  # Adjust for MATLAB's 1-based indexing
        cols = data[:, 1].astype(int) - 1  # Adjust for MATLAB's 1-based indexing
        values = data[:, 2]
        
        n_genes = rows.max() + 1
        n_cells = cols.max() + 1
        
        X = scipy.sparse.coo_matrix((values, (rows, cols)), shape=(n_genes, n_cells))
        return X.tocsr()
    except Exception as e:
        print(f"Error loading sparse matrix from {file_path}: {e}")
        sys.exit(1)

# Function to load CSV files with error handling
def load_csv(file_path, header=None, dtype=None):
    if os.path.exists(file_path):
        return pd.read_csv(file_path, header=header, dtype=dtype)
    else:
        print(f"Error: {file_path} does not exist.")
        sys.exit(1)

# --- 1. Load Main Expression Matrix ---
inputXf = 'sce_X.csv'
if not os.path.exists(inputXf):
    print(f"Error: Input file {inputXf} does not exist.")
    sys.exit(1)

X = load_sparse_matrix(inputXf)
print(f"Count matrix shape (genes x cells): {X.shape}")
X = X.T  # Transpose to get (cells x genes) for AnnData
print(f"Transposed matrix shape (cells x genes): {X.shape}")

# --- 2. Load Gene Information (var) ---
gene_file = "sce_genes.csv"
genes_df = load_csv(gene_file, dtype=str)
# Use the first column as the index for the var DataFrame
var_df = pd.DataFrame(index=genes_df.iloc[:, 0].values)

# --- 3. Load Core Cell Metadata (obs) ---
cell_ids = load_csv("sce_cell_ids.csv", dtype=str).values.flatten()
cell_types = load_csv("sce_celltypes.csv", dtype=str).values.flatten()
batch_ids = load_csv("sce_batch.csv", dtype=str).values.flatten()
clusters = load_csv("sce_clusters.csv").values.flatten()

# --- START OF CORRECTED SECTION ---
# Create the initial DataFrame. Pass raw numpy arrays and let the DataFrame
# constructor align them with the 'index' argument.
obs_df = pd.DataFrame({
    "cell_type": cell_types,
    "batch_id": batch_ids,
    "cluster": clusters
}, index=cell_ids)

# It's good practice to convert string-like columns to 'category' for memory efficiency
obs_df['cell_type'] = obs_df['cell_type'].astype('category')
obs_df['batch_id'] = obs_df['batch_id'].astype('category')
obs_df['cluster'] = obs_df['cluster'].astype('category')

# --- END OF CORRECTED SECTION ---

print("Loaded core metadata. obs shape:", obs_df.shape)

# --- 4. DISCOVER AND LOAD ADDITIONAL METADATA (NEW SECTION) ---
# Find all files matching the pattern *_sce_metadata.csv
metadata_files = glob.glob('*_sce_metadata.csv')
print(f"Found {len(metadata_files)} additional metadata file(s): {metadata_files}")

for file_path in metadata_files:
    try:
        # Extract the attribute name from the filename
        # e.g., 'SubCellType_sce_metadata.csv' -> 'SubCellType'
        attribute_name = os.path.basename(file_path).removesuffix('_sce_metadata.csv')
        
        print(f"Loading attribute '{attribute_name}' from {file_path}...")
        
        # Load the data, ensuring it's treated as string/categorical
        additional_data = load_csv(file_path, dtype=str).values.flatten()
        
        # Add it as a new column to the obs dataframe
        if len(additional_data) == len(obs_df):
            obs_df[attribute_name] = pd.Series(additional_data, index=obs_df.index, dtype="category")
        else:
            print(f"Warning: Skipping '{attribute_name}' due to length mismatch. "
                  f"Expected {len(obs_df)}, got {len(additional_data)}.")

    except Exception as e:
        print(f"Could not process metadata file {file_path}: {e}")

print("Final cell metadata (obs) columns:", obs_df.columns.tolist())
print(obs_df.head())

# --- 5. Load Embeddings (obsm) ---
embeddings = load_csv("sce_embeddings.csv").values

# --- 6. Create the AnnData Object ---
adata = ad.AnnData(
    X=X,
    obs=obs_df,
    var=var_df,
    obsm={"X_umap": embeddings}
)
print("\nAnnData object created successfully:")
print(adata)

# --- 7. Save and Clean Up ---
output_file = "sce_data.h5ad"
try:
    adata.write(output_file, compression="gzip")
    print(f"\nAnnData object saved to {output_file}")
except Exception as e:
    print(f"Error saving AnnData: {e}")
    sys.exit(1)

# Cleanup: Delete the temporary CSV files
files_to_delete = [
    'sce_X.csv', 'sce_genes.csv', 'sce_embeddings.csv', 'sce_clusters.csv',
    'sce_celltypes.csv', 'sce_batch.csv', 'sce_cell_ids.csv'
]
# Add the dynamically found metadata files to the cleanup list
files_to_delete.extend(metadata_files)

print("\nCleaning up temporary files...")
for file in files_to_delete:
    try:
        if os.path.exists(file):
            os.remove(file)
            # print(f" - Deleted {file}") # Uncomment for verbose cleanup
    except Exception as e:
        print(f"Could not delete {file}: {e}")

print("Done.")
