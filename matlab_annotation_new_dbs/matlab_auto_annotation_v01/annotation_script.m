std_bool = true;
% 1: panglao_2021, 2: augmented_2021, 3: cellmarker_2024
type_marker = 1;

if type_marker == 1
    markers = "PanglaoDB_Augmented_2021.txt";
elseif type_marker == 2
    markers = "CellMarker_Augmented_2021.txt";
elseif type_marker == 3
    markers = "CellMarker_2024.txt";
else
    type_marker = 3;
    markers = "CellMarker_2024.txt";
end
fprintf("Utilizing : %s \n",markers)

T = readtable(markers,'ReadVariableNames',false,'Delimiter','tab');
%T.Properties.VariableNames = {};
T = T(:,~all(ismissing(T)));
T = convertvars(T,@iscellstr,"string");
T.Var1 =  lower( string(T.Var1));

if type_marker == 3
    organism = "mouse"; % mouse | human
    T = T(contains(T.Var1, organism),: );
    tissue_types = ["colon" "intestin"];
    T = T(contains(T.Var1, tissue_types),: );
elseif type_marker == 2 
    tissue_types = ["colon" "intestin"];
    T = T(contains(T.Var1, tissue_types),: );
elseif type_marker == 1 
    tissue_types = ["stem" "epithelial" "enterocytes" "tuft" ...
                    "t cells" "t memory cells" "goblet" ...
                    "neurons" "enteroendocrine" "b cells" "beta" ...
                    "enterochromaffin"  "fibroblast" "endothelial"...
                    "smooth muscle" "dendritic" "macrophages" "plasma cells"...
                    "schwann"];
    % tissue_types = ["goblet" "enterocytes" "stem" "smooth muscle" "b cells"...
    %                 "t cells" "dendritic" "marcophages" "schwann" "neurons"...
    %                 "plasma cells" "enteroendocrine" "fibroblasts" ];
    % tissue_types = [ "t cells" "b cells" "beta" "ductal" "endothelial" "enteroendocrine" ...
    %                  "fibroblast" "erythroid" "mural" "alpha" "smooth muscle"];
    T = T(contains(T.Var1, tissue_types),: );

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

fprintf("Pre-scoring unique genes of database \n");
% ntot_genes are all the possible genes in the dataset
ntot_genes = size(genes,1);
scores = zeros( ntot_genes, ndb);
for ig = 1:ntot_genes
    for idb = 1:ndb
        bool = genes(ig)== table2array( Ttmp(idb,1:iends(idb)) );
        scores(ig,idb) = sum(bool); 
    end
    %score_sum(ig) = sum(scores(ig,:));
    val = sum(scores(ig,:));
    scores(ig,:) = scores(ig,:)./val;
end
clear Ttmp bool;
% Finished database preparation... we may store this to do not redo it


fprintf("Intersected genes within DB and SCE \n");
% Dataset look for common genes within DB 
sce.g = upper(sce.g);
[common_genes, idx, idx2 ]= intersect(sce.g, genes, 'stable');
genes_sce = sce.g(idx);
X = full( sce.X(idx,:) );
% Re-mapping information
genes = genes(idx2);
scores = scores(idx2,:);

% Scores in db x gene now
scores = scores';

clear idx idx2;

% Get cluster ids
sce = sce.embedcells('umap3d', true, false, 3);
sce = sce.clustercells([], [], true);

ncell = size(X,2);
cell_types = T.Var1;
clusters = unique(sce.c_cluster_id);
%cell_annotation = strings(ncell,1);
nclus = size(clusters,1);
fprintf("Annotating %d cells in %d clusters... \n",ncell, nclus);

score_record = zeros(nclus,5);
db_record = strings(nclus,5);
for iclus = 1:nclus
    cell_idx = find(clusters(iclus)==sce.c_cluster_id);
    Xtmp = X(:,cell_idx);
    if std_bool
        Xtmp = zscore(Xtmp);
    end
    Xtmp = scores*Xtmp;
    sum_score = sum(Xtmp,2);
    [maxval, idx] = maxk(sum_score,5);
    score_record(iclus,:) = log(maxval'+1);
    db_record(iclus,:) = cell_types(idx);
    sce.c_cell_type_tx(cell_idx) = cell_types(idx(1));
end

% 
% db_weights = scores'*X;
% cell_annotation = strings(ncell,1);
% cell_types = T.Var1;
% for icell = 1:ncell
%     [maxval, idx]=max(db_weights(:,icell));
%     cell_annotation(icell) = cell_types(idx);
% end
% 
% sce.c_cell_type_tx = cell_annotation;




