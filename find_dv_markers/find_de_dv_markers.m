function find_de_dv_markers(sce, fname, cells_per_clus)
    % INPUTS:
    % sce -----------------> SCE object
    % fname ---------------> file name to save processed SCEs
    % cells_per_cluster ---> cells per cluster (clusters =
    %                                           numCells/cells_per_cluster)
    %
    % OUTPUT:
    % Files...
    % Example: find_dv_markers(sce,'example_mouse',500,true)
    if nargin <2 || isempty(fname)
        fname = 'results';
    end
    if nargin < 3 || isempty(cells_per_clus)
        cells_per_clus = 500;
    end
    
    %sce = sce.qcfilter;
    sce = sce.embedcells('tsne3d', true, false, 3);
    if cells_per_clus ~= 0
        nclusters = round(sce.numcells/cells_per_clus );
        clust_type = "kmeans";
        fprintf("Working on %d clusters %s \n",nclusters, clust_type);
        sce = sce.clustercells(nclusters, clust_type, true);
    else
        nclusters = length(unique(sce.c_cluster_id));
        fprintf("Working on %d clusters and NO recluster\n", nclusters);
    end
    
    clusters = unique(sce.c_cluster_id);
    nclus = size(clusters,1);
    
    ngenes = length(sce.g);
    % Get gene indices
    dv_scored_g_tmp = strings(nclus, ngenes);
    tic
    for iclus = 1:nclus
        prog = floor(iclus/nclus*100);
        textprogressbar(prog);
    
        % Get ith-cluster cells
        cell_idx = find(clusters(iclus) == sce.c_cluster_id);
        [X, g] = sc_qcfilter0(sce.X(:,cell_idx), sce.g);
        [X, g] = sc_selectg(X, g, 1, 0.05);
    
        % Get NONE ith-cluster cells
        cell_jdx = find(clusters(iclus)~=sce.c_cluster_id);
        [X2, g2] = sc_qcfilter0(sce.X(:,cell_jdx), sce.g);
        [X2, g2 ] = sc_selectg(X2, g2, 1, 0.05);
        
        % Work only with common genes
        [gl, irows, jrows] = intersect(g, g2, 'stable');
        X = X(irows,:);
        X2 = X2(jrows,:);
    
        % Assign a boolean to normalize
        [X] = sc_norm(X,'type','libsize');
        [X2]= sc_norm(X2,'type','libsize');
    
        % Spline and gene statistics for X
        [T1, X, g, xyz1] = sc_splinefit_new(X, gl, true, false);
        [T1, idx1] = sortrows(T1,'genes','ascend');
        X = X(idx1, :);
        g = g(idx1);
    
        [T2, X2, g2, xyz2] = sc_splinefit_new(X2, gl, true, false);
        [T2, idx2] = sortrows(T2,'genes','ascend');
        X2 = X2(idx2, :);
        g2 = g2(idx2);
    
        px1 = T1.lgu; py1 = T1.lgcv; pz1 = T1.dropr;
        px2 = T2.lgu; py2 = T2.lgcv; pz2 = T2.dropr;
    
        v1=([px1 py1 pz1] - xyz1(T1.nearidx,:));
        v2=([px2 py2 pz2] - xyz2(T2.nearidx,:));
    
        diff = v1 - v2;
        DiffDist = vecnorm(diff, 2, 2);
        %DiffSign = sign(vecnorm(v2,2,2)-vecnorm(v1,2,2));
        %DiffSign = sum(sign(diff(:,1)),2);
        DiffSign = T1.lgu - T2.lgu;
    
        T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%d', iclus));
        T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%d', iclus+nclus));
        
        T = [T1 T2 table(DiffDist) table(DiffSign)];
    
        % Filtering/trimming ends
        idxx = T.(8)==1 | T.(16)==1 | T.(8) == max(T.(8)) | T.(16) == max(T.(16));
        T.DiffDist(idxx) = 0;
        % Sort by DV distance
        T = sortrows(T,"DiffDist","descend");
    
        % Get only "Up-regulated" genes in favor of ith-cluster by sign direction
        % test_val = [1.0 0.5 0.2 0.1 0];
        % for ival = 1:5
        %     idx = table2array(T(:,18)) > test_val(ival);
        %     if sum(idx) < 10
        %         continue;
        %     else
        %         %fprintf("success???")
        %         T = T(idx, : );
        %         break;
        %     end
        % end
        idx = table2array(T(:,18)) > 0;
        T = T(idx, : );
    
        % After scored and ordered
        ng = size(T,1);
        tmp = strings(1,ngenes);
        tmp(1:ng) = table2array( T(:,1) );
        dv_scored_g_tmp(iclus, :) = tmp;
        
    end
    time_end = toc;
    fprintf("DV find markers for clusters time: %f sec \n",time_end);
    dv_scored_g = dv_scored_g_tmp;
    writematrix(dv_scored_g(:,1:100), fname+"cluster_markers.csv");
    
    % New DE loop with consolidated Excel output
    tic
    fprintf("Starting DE analysis for clusters...\n");
    fname_excel = sprintf('%s_DE_all_clusters.xlsx', fname);
    
    for iclus = 1:nclus
        fprintf("Processing cluster %d/%d\n", iclus, nclus);
        % Get ith-cluster cells
        cell_idx = find(clusters(iclus) == sce.c_cluster_id);
        X = sce.X(:,cell_idx);
        g = sce.g;
        % Get NONE ith-cluster cells
        cell_jdx = find(clusters(iclus)~=sce.c_cluster_id);
        X2 = sce.X(:,cell_jdx);
        g2 = sce.g;
        % Filter for common genes
        [gl, irows, jrows] = intersect(g, g2, 'stable');
        X = X(irows,:);
        X2 = X2(jrows,:);
        % Call the DE function
        [Tde, Tup, Tdn] = sc_deg(X, X2, gl, 1, false);
        
        % Save DE results to different sheets in the same Excel file
        sheetName_all = sprintf('Cluster_%d_All_Genes', clusters(iclus));
        writetable(Tde, fname_excel, 'Sheet', sheetName_all);
        
        sheetName_up = sprintf('Cluster_%d_Up', clusters(iclus));
        writetable(Tup, fname_excel, 'Sheet', sheetName_up);
        
        sheetName_dn = sprintf('Cluster_%d_Down', clusters(iclus));
        writetable(Tdn, fname_excel, 'Sheet', sheetName_dn);
    end
    time_end_de = toc;
    fprintf("DE analysis for clusters time: %f sec \n", time_end_de);
end