import scanpy as sc
import pandas as pd
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

def process_h5_data( path_file, output_csv="leiden_umap.csv",  npca=50,
                     ndim=2, use_hvg=False, my_res=1, mt_pct=5, min_genes0=500,
                     min_cells0=15, max_counts=100000, algo_cluster="leidenalg"):
    """
    Process single-cell data from a 10x HDF5 file and perform Leiden clustering with UMAP embedding.

    Args:
        path_file (str): Path to the 10x HDF5 file.
        output_csv (str): Path to save the Leiden and UMAP data as CSV.
        npca (int): Number of principal components to compute.
        ndim (int): Number of dimensions for UMAP (2 or 3).
        use_hvg (bool): Whether to use highly variable genes for PCA.
        my_res (float): Resolution for Leiden clustering.
        mt_pct (float): Maximum mitochondrial content percentage.
        min_genes0 (int): Minimum number of genes per cell.
        min_cells0 (int): Minimum number of cells per gene.
        max_counts (float): Maximum total counts per cell.
        algo_cluster (str): Clustering algorithm flavor ("leidenalg" or "igraph").

    Returns:
        None: Saves the clustering and UMAP data to the specified CSV file.
    """
    # Step 1: Load the data   
    adata = sc.read_10x_h5(path_file)

    adata.var_names_make_unique()  # Ensure unique variable names

    # Step 2: Annotate mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-", "Mt-"))

    # Step 3: Compute QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Step 4: Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=min_genes0)
    sc.pp.filter_genes(adata, min_cells=min_cells0)
    adata = adata[adata.obs.pct_counts_mt < mt_pct, :].copy()
    adata = adata[adata.obs.total_counts < max_counts, :].copy()

    # Step 5: Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    # Step 6: Identify highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    # Step 7: PCA (handles large datasets with `chunked=True` if needed)
    if adata.n_obs > 20000:
        sc.tl.pca(adata, n_comps=npca, use_highly_variable=use_hvg, chunked=True, chunk_size=10000, svd_solver="arpack")
    else:
        sc.tl.pca(adata, n_comps=npca, use_highly_variable=use_hvg, chunked=False, svd_solver="arpack")

    # Step 8: Neighbors and UMAP
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=npca)
    sc.tl.umap(adata, n_components=ndim)

    # Step 9: Leiden Clustering
    sc.tl.leiden(adata, flavor=algo_cluster, n_iterations=2, resolution=my_res)

    # Step 10: Extract clustering and UMAP data
    leiden = adata.obs["leiden"]
    umap = adata.obsm["X_umap"]

    # Step 11: Save to CSV
    if ndim == 2:
        data = pd.DataFrame({
            "leiden": leiden,
            "UMAP_1": umap[:, 0],
            "UMAP_2": umap[:, 1]
        })
    elif ndim == 3:
        data = pd.DataFrame({
            "leiden": leiden,
            "UMAP_1": umap[:, 0],
            "UMAP_2": umap[:, 1],
            "UMAP_3": umap[:, 2]
        })

    data.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")

    
# process_h5_data(
#     path_file="C:\\Users\\selim\\Documents\\vs_code_working_dir\\scGeatoobox_addson\\UMAP_python\\sample_filtered_feature_bc_matrix.h5",
#     output_csv="leiden_umap.csv",
#     npca=50,
#     ndim=2,
#     use_hvg=False,
#     my_res=1,
#     mt_pct=5,
#     min_genes0=500,
#     min_cells0=15,
#     max_counts=100000,
#     algo_cluster="leidenalg"
# )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process 10x HDF5 single-cell data and perform clustering with UMAP.")
    parser.add_argument("path_file", type=str, help="Path to the 10x HDF5 file.")
    parser.add_argument("--output_csv", type=str, default="leiden_umap.csv", help="Output CSV file for results.")
    parser.add_argument("--npca", type=int, default=50, help="Number of principal components for PCA.")
    parser.add_argument("--ndim", type=int, default=2, help="Number of UMAP dimensions (2 or 3).")
    parser.add_argument("--use_hvg", action="store_true", help="Use highly variable genes for PCA.")
    parser.add_argument("--my_res", type=float, default=1.0, help="Resolution for Leiden clustering.")
    parser.add_argument("--mt_pct", type=float, default=5, help="Maximum mitochondrial content percentage.")
    parser.add_argument("--min_genes0", type=int, default=500, help="Minimum number of genes per cell.")
    parser.add_argument("--min_cells0", type=int, default=15, help="Minimum number of cells per gene.")
    parser.add_argument("--max_counts", type=float, default=100000, help="Maximum total counts per cell.")
    parser.add_argument("--algo_cluster", type=str, default="leidenalg", help="Clustering algorithm ('leidenalg' or 'igraph').")
    
    args = parser.parse_args()
    
    process_h5_data(path_file,
        output_csv=args.output_csv,
        npca=args.npca,
        ndim=args.ndim,
        use_hvg=args.use_hvg,
        my_res=args.my_res,
        mt_pct=args.mt_pct,
        min_genes0=args.min_genes0,
        min_cells0=args.min_cells0,
        max_counts=args.max_counts,
        algo_cluster=args.algo_cluster
    )