% Standarize cluster matrix to score cells?
std_bool = true;

% Marker database file
markers = "super_markers.txt";

T = readtable(markers,'ReadVariableNames',false,'Delimiter','tab');
T = convertvars(T,@iscellstr,"string");
T.Var2 =  lower( string(T.Var2));
 
% organism : "mouse" "human" "NONE"
organism = "NONE";

% Type markers: "cellmarker2024" "cellmarker_aug2021" "panglao2020" "panglao_aug2021" "NONE"
type_marker = "panglao_aug2021" ; 

% Tissue type : ... not guaranteed "NONE"
tissue_type = "NONE";

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

size(T)

target_tissue = [];

% % Immune cells 
% target_tissue = [ "macrophages" "mast" "b cells" "dendritic" ...
%                   "basophils" "gamma" "megakaryocytes" "monocytes" ...
%                   "myeloid" "nk" "neutrophils" "nuocytes" "plasma" ...
%                   "t cell" "red" "t memory"];
% 
% % Other?
% target_tissue =[ "stem" "epithelial" "enterocytes" "tuft" ...
%                     "goblet" "neurons" "enteroendocrine" "beta" ...
%                     "enterochromaffin"  "fibroblast" "endothelial"...
%                     "smooth muscle" "schwann" target_tissue];
%T = T(contains(T.Var1, target_tissue),: );


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
sce = sce.embedcells('umap2d', true, false, 3);
sce = sce.clustercells([], [], true);

ncell = size(X,2);
cell_types = T.Var2;
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




