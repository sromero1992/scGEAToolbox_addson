function sce_tmp = annotation_levels_gca(sce_tmp, nlevels, res)
    if nargin < 2; nlevels = 2; end
    if nargin < 3; res = 1.0; end

    % Read and process the database file only once
    dbs_wd = which('annotation_levels_gca');
    dbs_wd = erase(dbs_wd, 'annotation_levels_gca.m');
    dbs_wd = strcat(dbs_wd, 'gca_markers.txt');
    if ispc
        dbs_wd = strrep(dbs_wd, "\", "\\");
    end
    T = readtable(dbs_wd, 'ReadVariableNames', false, 'Delimiter', 'tab');
    T = convertvars(T, @iscellstr, "string");
    T.Var2 = lower(string(T.Var2));

    % Remove db info
    T = T(:, 2:end);
    sce_tmp.c_cell_type_tx = repmat("", length(sce_tmp.c_cluster_id), 1);

    % Perform level 1 annotation
    fprintf("Starting level 1 annotation of gut cell atlas markers \n");
    all_cell_types_db = T{:, 1};
    split_cell_types = cellfun(@(x) split(x, "||"), all_cell_types_db, 'UniformOutput', false);
    fprintf("Preparing unique genes in database \n");
    cell_array2 = cellfun(@(x) x{1}, split_cell_types, 'UniformOutput', false);
    gene_start = 0;
    gene_end = 10;
    current_cell_types_db = string(unique(cell_array2));
    idx_lvl1 = false(length(all_cell_types_db), 1);
    for i = 1:length(current_cell_types_db)
        temp_idx = contains(all_cell_types_db, current_cell_types_db(i));
        if any(temp_idx)
            first_match_index = find(temp_idx, 1);
            idx_lvl1(first_match_index) = true;
        end
    end
    Tsub = T(idx_lvl1, gene_start + 2:gene_end + 1);
    genes = extractUniqueGenes(Tsub);
    scores = calculateGeneScores(Tsub, genes);
    fprintf("Intersected genes within DB and sce \n");
    gtmp = upper(sce_tmp.g);
    genes = upper(genes);
    [genes, idx, idx2] = intersect(gtmp, genes, 'stable');
    X = pkg.norm_libsize(full(sce_tmp.X(idx, :)), 1e4);
    X = log1p(X);
    % Xp = sc_transform(X);
    % if max(max(Xp)) > 0
    %     X = Xp;
    % else
    %     fprintf("Pearson residual transformation failed, rolling back...\n");
    % end
    % clear Xp
    X = sparse(X);
    scores = scores(idx2, :)';
    clear idx idx2;
    sce_tmp = leiden_clustering_ann(sce_tmp, res);
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
        fprintf("Identified cell type: %s \n", current_cell_types_db(idx));
        sce_tmp.c_cell_type_tx(cell_idx) = current_cell_types_db(idx);
    end

    % Perform level 2 annotation
    if nlevels > 1
        fprintf("Starting level 2 annotation of gut cell atlas markers \n");
        unique_level1_annotations = unique(sce_tmp.c_cell_type_tx);
        for i = 1:length(unique_level1_annotations)
            level1_annotation = unique_level1_annotations(i);
            idx_level1_cells = sce_tmp.c_cell_type_tx == level1_annotation;
            sce_tmp_subset = sce_tmp.selectcells(idx_level1_cells);
            sce_tmp_subset = leiden_clustering_ann(sce_tmp_subset, res);

            idx_levl1_tab = contains(all_cell_types_db, level1_annotation);
            gene_start = 10;
            sub_cell_types_db = all_cell_types_db(idx_levl1_tab);
            Tsub = T(idx_levl1_tab, gene_start + 2:end);

            genes = extractUniqueGenes(Tsub);
            scores = calculateGeneScores(Tsub, genes);
            fprintf("Intersected genes within DB and sce \n");
            gtmp = upper(sce_tmp_subset.g);
            genes = upper(genes);
            [~, idx, idx2] = intersect(gtmp, genes, 'stable');
            X = pkg.norm_libsize(full(sce_tmp_subset.X(idx, :)), 1e4);
            X = log1p(X);
            Xp = sc_transform(X);
            if max(max(Xp)) > 0
                X = Xp;
            end
            clear Xp
            X = sparse(X);
            scores = scores(idx2, :)';
            clear idx idx2;
            ncell = size(X, 2);
            clusters = unique(sce_tmp_subset.c_cluster_id);
            nclus = size(clusters, 1);
            fprintf("Annotating %d cells in %d clusters... \n", ncell, nclus);
            for iclus = 1:nclus
                fprintf("Working on cluster %d \n", iclus);
                cell_idx = find(clusters(iclus) == sce_tmp_subset.c_cluster_id);
                Xtmp = full( X(:, cell_idx) );
                Xtmp = scores * Xtmp;
                sum_score = sum(Xtmp, 2,'omitnan');
                [~, idx] = maxk(sum_score, 1);
                fprintf("Identified cell type: %s \n", sub_cell_types_db(idx));
                sce_tmp_subset.c_cell_type_tx(cell_idx) = sub_cell_types_db(idx);
            end
            sce_tmp.c_cell_type_tx(idx_level1_cells) = sce_tmp_subset.c_cell_type_tx;
        end
    end

    cell_type_tx = sce_tmp.c_cell_type_tx;
    split_cell_types = cellfun(@(x) split(x, "||"), cellstr(cell_type_tx), 'UniformOutput', false);

    % Find the maximum number of levels (columns)
    max_levels = max(cellfun(@length, split_cell_types));
    % Initialize the string matrix
    string_matrix = repmat("", length(cell_type_tx), max_levels);
    % Populate the string matrix
    for i = 1:length(split_cell_types)
        for j = 1:length(split_cell_types{i})
            string_matrix(i, j) = split_cell_types{i}{j};
        end
    end

    % Saving all the levels of annotation
    sce_tmp.list_cell_attributes{end+1} = 'cell_type_level1'; 
    sce_tmp.list_cell_attributes{end+1} = string_matrix(:,1); 
    if nlevels> 1
        sce_tmp.list_cell_attributes{end+1} = 'cell_type_level2'; 
        sce_tmp.list_cell_attributes{end+1} = string_matrix(:,2); 
        sce_tmp.list_cell_attributes{end+1} = 'cell_type_level3'; 
        sce_tmp.list_cell_attributes{end+1} = string_matrix(:,3); 
    end
    % Reconstruct sce_tmp.c_cell_type_tx from string_matrix
    num_cells = size(string_matrix, 1);
    sce_tmp.c_cell_type_tx = repmat("", num_cells, 1); % Initialize c_cell_type_tx

    for i = 1:num_cells
        % Find the last non-empty level
        last_level = find(string_matrix(i, :) ~= "", 1, 'last');
        %last_level = 2;
        % Reconstruct the string
        if ~isempty(last_level)
            sce_tmp.c_cell_type_tx(i) = string_matrix(i, last_level);
        end
    end    
    unique_annotations = unique(sce_tmp.c_cell_type_tx);
    disp(unique_annotations);
    fprintf("Annotation finished!!!\n");
end
