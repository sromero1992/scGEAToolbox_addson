function sce = rm_non_inf_genes(sce, species)
    % This function removes non informative genes by species, it finds
    % Biomart_${species}_gene.mat defined in the toolbox that should be
    % located under Documents that MATLAB can read.
    % AUTHOR: Selim Romero
    % INPUT: 
    % sce -----> SCE object defined by scGEAToolbox
    % species--> Species either 'mouse' or 'human'
    % OUTPUT:
    % sce -----> SCE object with removed non informative genes
    % USAGE: 
    % sce_tmp = rm_non_inf_genes(sce_tmp, 'mouse');
    
    dir_file = strcat('Biomart_', species, '_genes.mat');
    
    % Define the root folder where you want to start searching (e.g., Documents folder)
    root_folder = fullfile(getenv('USERPROFILE'), 'Documents');
    
    % Search for the file recursively in the specified root folder
    search_pattern = fullfile(root_folder, '**', dir_file); % '**' searches all subdirectories
    file_info = dir(search_pattern);
    
    if isempty(file_info)
        error('File not found: %s', dir_file);
    else
        % If the file is found, construct the full path
        ref_genes_wd = fullfile(file_info(1).folder, file_info(1).name);
        % Load the file
        load(ref_genes_wd, 'T');
    end
    
    ApprovedSymbol = string(T.GeneName);
               
    sce = sce.rmmtgenes;
    sce = sce.rmhemoglobingenes;
    sce = sce.rmribosomalgenes;
    [idx] = ~ismember(upper(sce.g), upper(ApprovedSymbol));
    if any(idx)
        sce.g(idx) = [];
        sce.X(idx, :) = [];
    end
end