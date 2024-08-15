home="/home/ssromerogon/";
addpath( genpath( home + "scGEAToolbox" ) );

path0 =  'mat_files';
a = dir('mat_files/*.mat');

SCE=cell(size(a,1),1);
for k=1:size(a,1)
    fprintf(' processing %d \n', k)
    fname = fullfile( a(k).folder, a(k).name);
    %disp(fname) 
    load(fname);
    grpid = strrep(a(k).name,'.mat','');
    sce.c_batch_id = string(repmat(grpid, sce.NumCells,1));
    SCE{k}=sce;
end
fprintf("Merging files... \n");
sce=sc_mergesces(SCE);
% Preprocessing...
% QC filtering stringent
sce = sce.qcfilter;
% Embedding cells
use_hvgs = false;
force_write = true;
emb_type ='umap3d';% 'tsne3d';
sce = sce.embedcells(emb_type, force_write, use_hvgs, 3);
% Clusters cells
% Method to cluster cells
method = 'kmeans'; % can be snndpc as well
%% number of cluster per cel default is k = ceil( sce.NumCells/100); approx 100 cells per cluster...
k = ceil(sce.NumCells/500); 
sce = sce.clustercells(k, 'kmeans', true);
%% Annotate cells with panglao DB 
species = 'mouse';
sce = sce.assigncelltype(species, false);
save('sce_merged','sce','-v7.3');
