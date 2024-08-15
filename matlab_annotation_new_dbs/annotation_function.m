function sce_final = annotation_function(sce_tmp,  type_marker, organism,  ...
                                    tissue_type, target_cells)
    
    if nargin < 2
        type_marker = "NONE";
        organism = "NONE";
        tissue_type = "NONE";
        target_cells = "NONE";
    elseif nargin < 3 
        organism = "NONE";
        tissue_type = "NONE";
        target_cells = "NONE";
    elseif nargin < 4
        tissue_type = "NONE";
        target_cells = "NONE";
    elseif nargin < 5
        target_cells = "NONE";  
    end

    % Type markers: "cellmarker2024" "cellmarker_aug2021" "panglao2020" "panglao_aug2021" "NONE"
    if isempty(type_marker)
        %type_marker = "cellmarker2024";
        type_marker = "NONE";
    end
    % organism : "mouse" "human" "NONE"
    if isempty(organism)
        organism = "NONE";
    end
    % Tissue type : ... not guaranteed, "NONE"
    if isempty(tissue_type)
        tissue_type = "NONE";
    end
    %tissue_type = ["colon", "intestin"];

    fprintf("Settings: type_marker %s , organism %s , tissue_type %s \n",...
              type_marker, organism, tissue_type);

    % Standarize cluster matrix to score cells?
    std_bool = true;
    % Marker database file
    markers = "C:\Users\ssromerogon\Documents\vscode_working_dir\new_cell_db\super_markers.txt";
    
    T = readtable(markers,'ReadVariableNames',false,'Delimiter','tab');
    T = convertvars(T,@iscellstr,"string");
    T.Var2 =  lower( string(T.Var2));
       
    % Look for specific database 
    idx = contains(T.Var1, type_marker);
    if any(idx)
        T = T(idx,: );
    end
    
    % Look for specific tissue_type
    idx = contains(T.Var2, tissue_type);
    if any(idx)
        T = T(idx,: );
    end
    
    % Look for specific organism
    idx = contains(T.Var2, organism);
    if any(idx)
        T = T(idx,: );
    end

    % target_cells =[ "stem" "epithelial" "enterocytes" "tuft" ...
    %                     "goblet" "neurons" "enteroendocrine" "beta" ...
    %                     "enterochromaffin"  "fibroblast" "endothelial"...
    %                     "smooth muscle" "schwann"];

    % Lymphoid
    %target_cells = [ "nk cells" "t cells" "b cells" "plasma cells"];
    % Myeloid
    %target_cells = [ "megakaryocyte" "erythrocyte" "mast" "basophil" "neutrophil" "eosinophil" "monocyte" "macrophages" "dendritic"];
    % Fibroblasts
    %target_cells = [ "fibroblasts" "myofibroblasts" "caf"];

    %Look for specific cell types
    idx = contains(T.Var2, target_cells);
    if any(idx)
        T = T( idx,: );
    end

    % Remove database information, may be needed somewhere else?
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

    fprintf("Preparing unique genes in database \n");
    % Pre-weight stage
    ndb = size(T,1);
    ngc = size(T,2);
    % Find all unique genes first
    genes = strings(ngc,1);
    iends = zeros(ndb,1);
    Ttmp = T(:,2:end);
    fprintf("DB Size : %d  %d \n",size(T))
    for idb = 1:1
        iends(idb) =  sum(table2array(Ttmp(idb,:)) ~= "");
        for jdb = idb + 1:ndb
            iends(jdb) =  sum(table2array(Ttmp(jdb,:)) ~= "");
            tmp1 = table2array(Ttmp(idb,1:iends(idb)));
            tmp2 = table2array(Ttmp(jdb,1:iends(jdb)));
            tmp = union( tmp1, tmp2 );
            genes = union( tmp, genes);
        end
        break;
    end
    clear tmp1 tmp2 tmp;
    
    % Remove voids
    genes(find(genes=="")) =[];
        
    fprintf("Pre-scoring unique genes of database... \n");
    % ntot_genes are all the possible genes in the dataset
    ntot_genes = size(genes,1);
    scores = zeros( ntot_genes, ndb);
    for idb = 1:ndb
        bool = ismember(genes, table2array( Ttmp( idb, 1:iends(idb) )) );
        scores(:,idb) = 1.0*bool;
    end
    for ig = 1:ntot_genes
        val = sum(scores(ig,:));
        scores(ig,:) = scores(ig,:)./val;
    end
    
    clear Ttmp bool;
    % Finished database preparation... we may store this to do not redo it
    
    
    fprintf("Intersected genes within DB and sce \n");
    % Dataset look for common genes within DB 
    
    gtmp = upper(sce_tmp.g);
    genes = upper(genes);
    [genes, idx, idx2 ]= intersect(gtmp, genes, 'stable');
    %gtmp = gtmp(idx);
    X = full( sce_tmp.X(idx,:) );
    % Re-mapping information
    %genes = genes(idx2);
    scores = scores(idx2,:);
    
    % Scores in db x gene now
    scores = scores';
    clear idx idx2;
    
    % Get cluster ids
    sce_tmp = sce_tmp.embedcells('umap3d', true, false, 3);

    cells_per_clus = 1000;
    nclusters = round(sce_tmp.numcells/cells_per_clus);
    clust_type = "kmeans";    
    fprintf("Working on %d clusters %s \n",nclusters, clust_type);
    sce_tmp = sce_tmp.clustercells(nclusters, clust_type, true);
    
    ncell = size(X,2);
    cell_types = T.Var2;
    cell_types = erase(cell_types,"human");
    cell_types = erase(cell_types,"mouse");
    clusters = unique(sce_tmp.c_cluster_id);
    nclus = size(clusters,1);
    
    fprintf("Annotating %d cells in %d clusters... \n",ncell, nclus);
    score_record = zeros(nclus,5);
    db_record = strings(nclus,5);
    markers_dv_save = strings(nclus,20);
    for iclus = 1:nclus
        fprintf("Working on cluster %d \n", iclus)
        cell_idx = find( clusters(iclus) == sce_tmp.c_cluster_id);
        if false % Find makers?
            % Findmarkers section
            markers_dv = findMarkers_DVcluster(sce_tmp, iclus, 5000);
            markers_dv_save(iclus,1:20) = markers_dv(1:20);
            fprintf("DV markers size : %d \n", length(markers_dv) )     
            [genes, idx, ~] = intersect(genes, markers_dv, 'stable');
            fprintf("overlaping DV markers with DB: %d \n", length(genes));
            Xtmp = X(idx, cell_idx);    
            scores = scores(:,idx);
        else    
            Xtmp = X(:, cell_idx);
        end
        % End of find markers gene swapping
        if std_bool
            Xtmp = zscore(Xtmp);
        end
        Xtmp = scores*Xtmp;
        sum_score = sum(Xtmp,2);
        [maxval, idx] = maxk(sum_score,5);
        score_record(iclus,:) = log1p(maxval');
        db_record(iclus,:) = cell_types(idx);
        fprintf("Identified cell type: %s \n", cell_types(idx(1)));  
        sce_tmp.c_cell_type_tx(cell_idx) = cell_types(idx(1));
    end
    T = table(clusters,db_record,score_record);
    writetable(T,"db_score_record.csv")
    writematrix(markers_dv_save,'markers_dv_clusters.csv')

    unique(sce_tmp.c_cell_type_tx)
    sce_final = sce_tmp;
    fprintf("Annotation finished!!!\n");
end



