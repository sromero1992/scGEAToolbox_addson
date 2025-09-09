function sce = assing_clus_cell(sce, iclus, cell_type)
    % Re annotation scratch
    %iclus = 16;
    %cell_type = 'Goblet unknown';
    idx = find(sce.c_cluster_id == iclus); 
    sce.c_cell_type_tx(idx) = cell_type; 
end