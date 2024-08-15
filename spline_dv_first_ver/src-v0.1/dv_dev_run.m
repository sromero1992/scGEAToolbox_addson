main_path = "HFD_LFD/";
sample_id = "HFD_LFD";
path = strcat(main_path,sample_id,".mat");
data  = load(path);
sce1 = data.sce;
clear data;

CellTypeList1 = unique(sce1.c_cell_type_tx);
celltype = "Adipocytes";
celltypef = replace(celltype," ","_");
label1 = strcat(celltypef,"-", sample_id);
fname1 = strcat(main_path,label1);

% Take out 5% mitochondrial cells and corresponding batch id

cellidx = and(sce1.c_batch_id == "HFD", sce1.c_cell_type_tx == celltype);
[X, g ] = sc_selectg(sce1.X(:,cellidx), sce1.g, 1, 0.05);
cellidx2 = and(sce1.c_batch_id == "LFD", sce1.c_cell_type_tx == celltype);
[X2,g2] = sc_selectg(sce1.X(:,cellidx2), sce1.g, 1, 0.05);

X = full(X);
X2 = full(X2);

% Work only with common genes
[gl, irows, jrows] = intersect(g, g2, 'stable');
X = X(irows,:);
X2 = X2(jrows,:);

% Assing a boolean to normalize
[X] = sc_norm(X,'type','libsize');
[X2]= sc_norm(X2,'type','libsize');

% Differential variability
[Tdv, Tig, influgenes] = sc_variability(X, X2, gl, gl, fname1, false);

% Differential expression
[Tde, Tup, Tdn] = sc_deg(X, X2, gl, 1, false);
label2 = strcat("DE-", label1);
fname2 = strcat(main_path, label2);
writetable(Tde, fname2, 'FileType', 'spreadsheet', 'Sheet', 'All genes');
writetable(Tup, fname2, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
writetable(Tdn, fname2, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');

%[g_up_dv, ~, ~] = intersect(Tdv.genes(1:200), Tup.gene(1:200), 'stable');
%[g_dn_dv, ~, ~] = intersect(Tdv.genes(1:200), Tdn.gene(1:200), 'stable');
[Tde, ~] = sortrows(Tde, 2, {'ascend'});

[g_de_dv, ~, ~] = intersect(Tdv.genes(1:200), Tde.gene(1:200), 'stable');