%clear
data = load('GSM3308547_GSM330854.mat');
%data = load('Colonocytes_06_03_2024.mat');
%data = load('Colonocytes_06_05_2024_shreyan.mat');

sce = data.sce;
clear data;

cells_per_clus = 250;

sce = sce.qcfilter;
sce = sce.embedcells('umap3d', true, false, 3);
if cells_per_clus ~= 0
    % nclusters = round(sce.numcells/cells_per_clus );
    % clust_type = "kmeans";
    % fprintf("Working on %d clusters %s \n",nclusters, clust_type);
    % sce = sce.clustercells(nclusters, clust_type, true);
    sce = leiden_annotation(sce, 'knn', 'mouse');
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
    test_val = [1.0 0.5 0.2 0.1 0];
    for ival = 1:5
        idx = table2array(T(:,18)) > test_val(ival);
        if sum(idx) < 10
            continue;
        else
            %fprintf("success???")
            T = T(idx, : );
            break;
        end
    end
    %idx = table2array(T(:,18)) > 0;
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
fname = "test_x";
writematrix(dv_scored_g(:,1:100), fname+"cluster_markers.csv");

% Save sce
save(fname +string(cells_per_clus) + ".mat",'sce','-v7.3');


if reduce_clusters
    fprintf("Reducing clusters to common DV (markers) gene sets \n");
    
    % Parameters
    max_genes0 = 500; % Max number of genes to consider
    min_genes = 10; % Minimum number of common genes required
    threshold = min_genes / max_genes0; % Threshold for gene intersection
    alloc = ceil(nclus * (nclus + 1) / 2); % Allocation size for arrays
    
    % Initialization
    equi_info = zeros(alloc, 3);
    equi_genes = strings(alloc, max_genes0);
    max_intergenes = 0;
    kclus = 0;
    
    % First Degree: Scoring intersections of DVs across clusters
    for iclus = 1:nclus
        for jclus = iclus + 1:nclus
            intergenes = intersect(dv_scored_g(iclus, 1:max_genes0), dv_scored_g(jclus, 1:max_genes0), 'stable');
            intergenes(intergenes == "") = [];
            n1 = length(intergenes);
            max_intergenes = max(n1, max_intergenes);
            rate = n1 / max_genes0;
            
            if rate >= threshold
                kclus = kclus + 1;
                equi_info(kclus, :) = [iclus, jclus, rate];
                equi_genes(kclus, 1:n1) = intergenes;
            end
        end
    end
    
    % Trim the pre-allocated arrays to the actual size
    equi_info = equi_info(1:kclus, :);
    equi_genes = equi_genes(1:kclus, 1:max_intergenes);
    
    % Initialize cluster reduction array
    clusters_reduction = sce.c_cluster_id;
    
    % Second Degree: Intersection of intersections
    ieq = 1;
    keq_redux = 0;
    max_eq = size(equi_info, 1);
    eq_genes_final = strings(alloc, max_genes0);
    cluster_red_id = zeros(alloc, 1);
    isgoing = true;
    min_genes = 2;
    
    % Better create 1 to want instead of vectorized, then after working,
    % make it vectorized
    while isgoing
        
        jeq = find(equi_info(ieq, 1) == equi_info(:, 1))';
        
        loc_info = equi_info(jeq, :);
        loc_genes = equi_genes(jeq, :);
        
        % Sort cluster scores in ieq representation
        [~, idx] = sort(loc_info(:, 3), 'descend');
        loc_info = loc_info(idx, :);
        loc_genes = loc_genes(idx, :);
        
        equi_info(jeq, :) = loc_info;
        equi_genes(jeq, :) = loc_genes;
        
        ninit = sum(equi_genes(ieq, :) ~= "");
        threshold = min_genes / ninit;
        
        % Obtain intersections of AnB vs CnD as (AnB)n(CnD)
        for meq = jeq
            intergenes = intersect(equi_genes(ieq, :), equi_genes(meq, :));
            intergenes(intergenes == "") = [];
            
            interg_bool = true;
            if meq == min(jeq)
                intergenes0 = intergenes;
            else
                intergenes0 = intersect(intergenes0, equi_genes(meq, :));
                if isempty(intergenes0) || length(intergenes0) < min_genes
                    intergenes0 = equi_genes_tmp;
                    interg_bool = false;
                end
            end
            
            ninter = length(intergenes);
            rate = ninter / ninit;
            
            if rate > threshold && interg_bool
                idx = equi_info(:,1) ==  equi_info(meq, 2);
                equi_info(meq, 2) = equi_info(ieq, 1);
                equi_info(idx,1) =

                equi_info(meq, 3) = -rate;
                equi_genes_tmp = intergenes0;
                jdx = clusters_reduction == equi_info(meq, 2);
                clusters_reduction(jdx) = equi_info(ieq, 1);
            else
                equi_info(meq, 1) = equi_info(meq, 2);
            end
        end
        keq_redux = keq_redux + 1;
        % Store cluster belonging to reduction
        cluster_red_id(keq_redux) = equi_info(ieq, 1);
        % Store genes that intersect (AnB)n(CnD)
        eq_genes_final(keq_redux, 1:length(equi_genes_tmp)) = equi_genes_tmp;
        ieq = max(jeq) + 1;
        isgoing = ieq < max_eq;
    end
    
    % Trim final arrays
    eq_genes_final = eq_genes_final(1:keq_redux, :);
    cluster_red_id = cluster_red_id(1:keq_redux);
    
    % % Third Degree: Additional merging based on common genes
    % min_genes = 3;
    % eq_genes_final2 = strings(kclus, max_genes0);
    % cluster_red_id2 = ones(kclus, 1);
    % 
    % for iclus = 1:kclus
    %     if cluster_red_id2(iclus) < 0
    %         continue;
    %     end
    %     eq_genes_final2(iclus, :) = eq_genes_final(iclus, :);
    % 
    %     for jclus = iclus + 1:kclus
    %         if cluster_red_id2(jclus) < 0
    %             continue;
    %         end
    % 
    %         intergenes = intersect(eq_genes_final2(iclus, :), eq_genes_final(jclus, :));
    %         intergenes(intergenes == "") = [];
    % 
    %         if length(intergenes) >= min_genes
    %             jdx = clusters_reduction == cluster_red_id(jclus);
    %             clusters_reduction(jdx) = cluster_red_id(iclus);
    %             nint = length(intergenes);
    %             eq_genes_final2(iclus, 1:nint) = intergenes;
    %             cluster_red_id2(jclus) = -1;
    %         end
    %     end
    % end
    
    idx = cluster_red_id2 > 0;
    eq_genes_final = eq_genes_final2(idx, :);
    cluster_red_id = cluster_red_id(idx);
    
    T = table(cluster_red_id, eq_genes_final);
    writetable(T, fname + "cluster_markers_dv_reduced.csv");
    
    % Second round reduction
    save(fname + string(cells_per_clus) + "_dv_reduced.mat", 'sce', '-v7.3');
end
sce.c_cluster_id = clusters_reduction;