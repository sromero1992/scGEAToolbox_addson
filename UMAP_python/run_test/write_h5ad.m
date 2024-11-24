function write_h5ad(sce)    

    % Save the adjacency matrix to a text file
    [i, j, val] = find(sce.X);
    writematrix([i, j, val], 'sce_X.csv', 'Delimiter', 'tab');% Gene expression matrix
    writematrix(sce.g, 'sce_genes.csv');             % Gene names
    writematrix(sce.s, 'sce_embeddings.csv');        % Embedding (UMAP, etc.)
    writematrix(sce.c, 'sce_clusters.csv');          % Cluster IDs
    writematrix(sce.c_cell_type_tx, 'sce_celltypes.csv'); % Cell type annotations
    writematrix(sce.c_batch_id, 'sce_batch.csv');    % Batch IDs
    writematrix(sce.c_cell_id, 'sce_cell_ids.csv');  % Cell IDs

end