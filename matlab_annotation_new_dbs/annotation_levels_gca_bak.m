function sce_tmp = annotation_levels_gca_bak(sce_tmp, nlevels,  target_cells)
    

    if nargin < 2; nlevel = 1; end
    if nargin < 3; target_cells = "NONE"; end

    fprintf("Starting level %d  annotation of gut cell atlas markers \n", nlevel);

    % Marker database file
    dbs_wd = which('annotation_levels_gca');
    dbs_wd = erase(dbs_wd, 'annotation_levels_gca.m');
    dbs_wd = strcat(dbs_wd, 'gca_markers.txt');
    if ispc
        dbs_wd = strrep(dbs_wd, "\", "\\");
    end
    
    T = readtable(dbs_wd,'ReadVariableNames',false,'Delimiter','tab');
    T = convertvars(T,@iscellstr,"string");
    T.Var2 =  lower( string(T.Var2));

    %Look for specific cell types
    idx = contains(T.Var2, target_cells, "IgnoreCase", true);
    if any(idx)
        T = T( idx,: );
    end

    % Remove database information
    T = T(:,2:end);

    ncols = size(T,2);
    for i = 1:ncols 
        icolstr = table2array(T(:,i));
        if all(icolstr == "")
            % Remove database info and voids
            T = T(:,1:i);
            break;
        end
    end

    cell_types = T{:,1};
    split_cell_types = cellfun(@(x) split(x, "||"), cell_types, 'UniformOutput', false);
    fprintf("Preparing unique genes in database \n");
    if nlevel >= 1 
        % Get the first element (which is a cell array)
        first_cell_array = cellfun(@(x) x{1}, split_cell_types, 'UniformOutput', false);

        % Find the union (unique elements) within the first cell array
        current_cell_types = string(unique(first_cell_array));

        % Initialize a logical array to store the indices
        idx = false(length(cell_types), 1);
        % Loop through each element in current_cell_types
        for i = 1:length(current_cell_types)
            % Find the indices where the current cell type is found
            temp_idx = contains(cell_types, current_cell_types(i));
            % If there are any matches, store the first one
            if any(temp_idx)
                first_match_index = find(temp_idx, 1); % Find the first index
                idx(first_match_index) = true; % Set the logical array to true at that index
            end
        end

        Tsub = T(idx,2:11);
        genes = extractUniqueGenes(Tsub);
        scores = calculateGeneScores(Tsub, genes);
        
        % Finished database preparation... we may store this to do not redo it
        fprintf("Intersected genes within DB and sce \n");
        % Dataset look for common genes within DB 
        gtmp = upper(sce_tmp.g);
        genes = upper(genes);
        [genes, idx, idx2 ]= intersect(gtmp, genes, 'stable');

        %gtmp = gtmp(idx);
        % This normalization in the sub set is on purpose
        % otherwise is memory expensive
        X =  pkg.norm_libsize(full(sce_tmp.X(idx,:)), 1e4);
        X = log1p(X);
        X = sc_transform(X);
        X = sparse(X);
    
        % Scores in db x gene now and re-mapping information
        scores = scores(idx2,:)';
        clear idx idx2;

    end
   
    sce_tmp = leiden_clustering_ann(sce_tmp, 2.0);
    ncell = size(X,2);
    clusters = unique(sce_tmp.c_cluster_id);
    nclus = size(clusters,1);
    sce_tmp.c_cell_type_tx = "";
    fprintf("Annotating %d cells in %d clusters... \n",ncell, nclus);
    for iclus = 1:nclus
        fprintf("Working on cluster %d \n", iclus)
        cell_idx = find( clusters(iclus) == sce_tmp.c_cluster_id);
        Xtmp = X(:, cell_idx);
        Xtmp = scores*Xtmp;
        sum_score = sum(Xtmp,2);
        [~, idx] = maxk(sum_score,1); %only get the maximum
        fprintf("Identified cell type: %s \n", current_cell_types(idx));  
        sce_tmp.c_cell_type_tx(cell_idx) = current_cell_types(idx);
    end
    unique(sce_tmp.c_cell_type_tx)
    fprintf("Annotation finished!!!\n");
end



