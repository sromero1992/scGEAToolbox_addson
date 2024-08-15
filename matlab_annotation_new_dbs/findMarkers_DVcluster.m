function dv_scored_g = findMarkers_DVcluster(sce0, iclus, topg)

    
    %dv_scored_g = strings(ngenes,1);

    % Get unique cluster id identifiers
    clusters = unique(sce.c_cluster_id);

    % Get ith-cluster cells
    cell_idx = find( clusters(iclus) == sce0.c_cluster_id);
    [X, g] = sc_qcfilter(sce0.X(:,cell_idx), sce0.g);
    [X, g] = sc_selectg(X, g, 1, 0.05);

    % Get NONE ith-cluster cells
    cell_jdx = find( clusters(iclus) ~= sce0.c_cluster_id);
    [X2, g2 ] = sc_qcfilter(sce0.X(:,cell_jdx), sce0.g);
    [X2, g2 ] = sc_selectg(X2, g2, 1, 0.05);
    
    % Work only with common genes
    [gl, irows, jrows] = intersect(g, g2, 'stable');
    X = X(irows,:);
    X2 = X2(jrows,:);

    % Assing a boolean to normalize
    [X] = sc_norm(X,'type','libsize');
    [X2]= sc_norm(X2,'type','libsize');

    % Spline and gene statistics for X
    [T1, X, g, xyz1] = sc_splinefit(X, gl, true, false);
    [T1, idx1] = sortrows(T1,'genes','ascend');
    % X = X(idx1, :);
    % g = g(idx1);

    [T2, X2, g2, xyz2] = sc_splinefit(X2, gl, true, false);
    [T2, idx2] = sortrows(T2,'genes','ascend');
    % X2 = X2(idx2, :);
    % g2 = g2(idx2);

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
    T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%d', iclus*2));
    
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

    % After scored and ordered
    ng = size(T,1);
    ming = min(ng,topg);
    dv_scored_g = table2array( T(1:ming,1) );

end