function sce = find_dv_markers0(sce, fname)
    sce = sce.qcfilter;
    %sce = sce.embedcells('umap3d', true, false, 3);
    sce = leiden_annotation(sce, 'knn', 'mouse');
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
    
        % Assing a boolean to normalize
        [X] = sc_norm(X,'type','libsize');
        [X2]= sc_norm(X2,'type','libsize');
    
        % Spline and gene statistics for X
        [T1, X, g, xyz1] = sc_splinefit_new(X, gl, true, false);
        [T1, idx1] = sortrows(T1,'genes','ascend');
        X = X(idx1, :);
        g = g(idx1);

        % Spline and gene statistics for X2
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
        DiffSign = sign(vecnorm(v1,2,2)-vecnorm(v2,2,2));

        T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%d', iclus));
        T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%d', iclus+nclus));
    
        T = [T1 T2 table(DiffDist) table(DiffSign)];
    
        % Filtering/trimming ends
        idxx = T.(8)==1 | T.(16)==1 | T.(8) == max(T.(8)) | T.(16) == max(T.(16));
        T.DiffDist(idxx) = 0;
        % Sort by DV distance
        T = sortrows(T,"DiffDist","descend");
    
        % % Get only "Up-regulated" genes in favor of ith-cluster by sign direction
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
    writematrix(dv_scored_g(:,1:30)', fname+"cluster_markers.csv");
    
    % Save sce
    save(fname + ".mat",'sce','-v7.3');
end