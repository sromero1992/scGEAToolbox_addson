home="/home/ssromerogon/";
addpath( genpath( home + "scGEAToolbox" ) );

D=readtable('h5_files/metadata.csv');
id = D.sample_id;
grpid = D.grpid;

target_dir = 'mat_files';
mkdir(target_dir);

for k = 1:length(id)
    fprintf('Processing %d ...\n', k);
    fname_target = fullfile( target_dir, string( id(k) ) + ".mat" ); 
    if exist(fname_target,"file"), continue; end
    
    fname_source = fullfile('h5_files',num2str(id(k)),'sample_filtered_feature_bc_matrix.h5');
    [X, g, b] = sc_read10xh5file(fname_source);
    sce=SingleCellExperiment(X,g);
    sce.c_cell_id = b;
    sce.c_batch_id = string( repmat( grpid(k), sce.NumCells, 1) );
    %sce = sce.onestepraw2anno('mouse');
    save(fname_target,'sce','-v7.3');
end
