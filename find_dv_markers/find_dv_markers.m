function find_dv_markers(sce, fname, cells_per_clus, reduce_clusters)
    % INPUTS:
    % sce -----------------> SCE object
    % fname ---------------> file name to save processed SCEs
    % cells_per_cluster ---> cells per cluster (clusters =
    %                                           numCells/cells_per_cluster)
    % reduce_clusters -----> reduce initial clusters according common DV
    %                        genes across clusters
    % 
    % OUTPUT:
    % Files...
    % Example: find_dv_markers(sce,'example_mouse',500,true) 

    if nargin <2 || isempty(fname)
        fname = 'dv_proccessed';
    end
    if nargin < 3 || isempty(cells_per_clus)
        cells_per_clus = 1000;
    end
    if nargin < 4 || isempty(reduce_clusters)
        reduce_clusters = false;
    end
    
    sce = sce.qcfilter;
    sce = sce.embedcells('umap3d', true, false, 3);
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
    
        % Assing a boolean to normalize
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

    % Save sce
    save(fname +string(cells_per_clus) + ".mat",'sce','-v7.3');

    if reduce_clusters
        % Max number of genes to take from DV genes (if all available)
        max_genes0 = 500;
        % Threshold for gene intersection
        threshold = 0.1; % 10%
        alloc = ceil( nclus * ( nclus + 1 )/ 2);
        equi_info = zeros(alloc, 3);
        eq_genes = strings(alloc, max_genes0);
        max_intergenes = 0;
        kclus = 0;
        % Scoring intersection of DVs across clusters
        for iclus = 1:nclus
            for jclus = iclus+1:nclus
                intergenes = intersect(dv_scored_g(iclus,1:max_genes0), ...
                                       dv_scored_g(jclus,1:max_genes0),'stable');
                % Remove void intersection""
                intergenes(intergenes=="")=[];
                n1 = length(intergenes);
                max_intergenes = max(n1, max_intergenes);
                rate = n1/max_genes0;
                if rate >= threshold
                    kclus = kclus + 1;
                    equi_info(kclus,:) =[iclus jclus rate];
                    eq_genes(kclus,1:n1) = intergenes;
                    %fprintf("iclus %d jclus %d with rate %f \n",iclus,jclus,rate);
                end
            end
        end
        
        % equi info contains icluster and jcluster intersected quantity/score
        equi_info = equi_info(1:kclus,:);
        % Set of genes that are intersected across icluster and jcluster
        eq_genes = eq_genes(1:kclus, 1:max_intergenes);
        
        % Get cluster ids from sce
        clusters_reduction = sce.c_cluster_id;
        
        % Second degree comparision intersection of intersections to score
        iclus = 1;
        kclus = 0;
        max_clus = length( equi_info(:,1) );
        eq_genes_final = strings(alloc, max_genes0);
        cluster_red_id = zeros(alloc,1);
        finish = true;
        while finish
        
            % Get clusters that has to do with iclus (the actual pair) 
            jclus = find( equi_info(iclus,1) == equi_info(:,1) )';
            loc_info = equi_info(jclus,:);
            loc_genes = eq_genes(jclus,:);
        
            % Sort clusters by score in iclus
            [~, idx] = sort( loc_info(:,3), 'descend');
            loc_info = loc_info(idx,:);
            loc_genes = loc_genes(idx,:);
        
            equi_info(jclus,:) = loc_info;
            eq_genes(jclus,:) = loc_genes;
        
            ninit =  length( eq_genes( iclus, eq_genes(iclus, :) ~= ""));
        
            for mclus = jclus
        
                %ninit =  length( eq_genes( iclus, eq_genes(iclus, :) ~= ""));
        
                intergenes = intersect( eq_genes(iclus,:),...
                                        eq_genes(mclus,:));
        
                intergenes(intergenes == "")=[];
        
                ninter = length(intergenes);
                % rate = ninter/ninit 
                %fprintf("ninit %d and ninter %d \n",ninit, ninter);
                if ninter/ninit > threshold 
                    eq_genes(iclus,:) = "";
                    eq_genes(mclus,:) = "";
                    eq_genes(iclus,1:ninter) = intergenes;
                    eq_genes(mclus,1:ninter) = intergenes;
                    clusters_reduction(clusters_reduction == equi_info(mclus,2)) = ...
                                        equi_info(iclus,1);
                    equi_info(mclus, 2) = equi_info(iclus,1);
                    equi_info(mclus, 3) = ninter;
                    nlastg = ninter;
                else
                    equi_info(mclus, 1) = equi_info(mclus,2);
                end
            end
            kclus = kclus + 1;
        
            cluster_red_id(kclus) = equi_info(iclus,1);
            eq_genes_final(kclus, 1:nlastg) = eq_genes(iclus,1:nlastg);
            iclus = max(jclus) + 1;
            finish = iclus < max_clus;
        end
        
        eq_genes_final = eq_genes_final(1:kclus,:);
        cluster_red_id = cluster_red_id(1:kclus);
        T = table(cluster_red_id, eq_genes_final);
        writetable(T, fname + "cluster_markers_dv_reduced.csv");    
        
        sce.c_cluster_id = clusters_reduction;
        save(fname +string(cells_per_clus) + "_dv_reduced.mat",'sce','-v7.3');
    end
end