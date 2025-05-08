function sce_tmp = annotation_function2(sce_tmp, res)
    if nargin < 2; res = 2.0; end

    % Read and process the database file only once
    dbs_wd = which('annotation_function2');
    dbs_wd = erase(dbs_wd, 'annotation_function2.m');
    dbs_wd = strcat(dbs_wd, 'markers.txt');
    if ispc
        dbs_wd = strrep(dbs_wd, "\", "\\");
    end
    T = readtable(dbs_wd, 'ReadVariableNames', false, 'Delimiter', 'tab');
    T = convertvars(T, @iscellstr, "string");
    T.Var2 = lower(string(T.Var2));

    % Remove db info
    T = T(:, 1:end);
    sce_tmp.c_cell_type_tx = repmat("", length(sce_tmp.c_cluster_id), 1);

    % Perform level 1 annotation
    fprintf("Starting annotation of custom markers \n");
    all_cell_types_db = T{:, 1};
    fprintf("Preparing unique genes in database \n");
    Tsub = T(:, 2:end);
    genes = extractUniqueGenes(Tsub);
    scores = calculateGeneScores(Tsub, genes);
    fprintf("Intersected genes within DB and sce \n");
    gtmp = upper(sce_tmp.g);
    genes = upper(genes);
    %[genes, idx, idx2] = intersect(gtmp, genes, 'stable');
    %X = pkg.norm_libsize(full(sce_tmp.X(idx, :)), 1e4);
    X = normalize_library_size_cell_chunks(sce_tmp.X);
    %[~, X, gtmp] = sc_splinefit(X, gtmp, true, false);
    % Obtain gener markers positions
    %[genes, idx, ~] = intersect(gtmp, genes, 'stable');
    % Merge gene marker positions to 2000 hvgs
    %range = union(1:2000, idx);
    %gtmp = gtmp(range);
    %X = X(range, :);
    X = log1p(X);
    % Will not be sparse after this
    %X = pearson_residuals_chunk(X);
    %X = sc_transform(X);
    [genes, idx, idx2] = intersect(gtmp, genes, 'stable');
    X = X(idx, :);
    scores = scores(idx2, :)';
    clear idx idx2;
    if res ~= 0 
        sce_tmp = leiden_clustering_ann(sce_tmp, res);
    end
    ncell = size(X, 2);
    clusters = unique(sce_tmp.c_cluster_id);
    nclus = size(clusters, 1);
    fprintf("Annotating %d cells in %d clusters... \n", ncell, nclus);
    for iclus = 1:nclus
        fprintf("Working on cluster %d \n", iclus);
        cell_idx = find(clusters(iclus) == sce_tmp.c_cluster_id);
        Xtmp = full(X(:, cell_idx));
        Xtmp = scores * Xtmp;
        sum_score = sum(Xtmp, 2,'omitnan');
        [~, idx] = maxk(sum_score, 1);
        fprintf("Identified cell type: %s \n", all_cell_types_db(idx));
        sce_tmp.c_cell_type_tx(cell_idx) = all_cell_types_db(idx);
    end

    unique_annotations = unique(sce_tmp.c_cell_type_tx);
    disp(unique_annotations);
    fprintf("Annotation finished!!!\n");
end
