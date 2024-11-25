function sce = umap_sc_wrapper(sce, ndim, use_hvgs, my_res)
    % umap_sc_wrapper computes UMAP and leiden clusters interfaced from
    % python3.11 scanpy library.
    % INPUT:
    % sce -------> SCE object 
    % ndim ------> n-dimentions of UMAP
    % OUTPUT:
    % sce ------> sce object containing Leiden clusters and UMAP
    % Usage:
    % sce = umap_sc_wrapper(sce);

    if nargin < 2 || isempty(ndim); ndim = 2; end
    if nargin < 3 || isempty(use_hvgs); use_hvgs = false; end
    if nargin < 4 || isempty(my_res); my_res = 1.0; end

    file_h5ad = 'sce_data.h5ad';
    file_out = 'leiden_umap.csv';
    npca = 50;
    mt_pct = 15.0;
    min_genes0 = 200;
    min_cells0 = 3;
    max_counts = 1.0e7;
    algo_cluster = 'leidenalg';

    %-----------------------------------------------------------------
    % Set the Python environment (Python 3.11)
    % Windows format
    %env_bin = 'F:\Anaconda\envs\scanpy_env_311\python.exe';
    env_bin = 'C:\Users\ssromerogon\.conda\envs\scanpy_env_311\python.exe';
    if ispc
        env_bin = strrep(env_bin,"\","\\");
    end
    % Linux format
    %env_bin = "/home/ssromerogon/packages/scanpy_env/bin/python3";
    %-----------------------------------------------------------------
    % Load python environment
    % Clear any existing Python environment to force reinitialization
    pe = pyenv('Version', env_bin);

    % Check if the environment is loaded
    if pe.Status ~= "Loaded"
        fprintf("Reinitializing Python environment...\n");
        pe = pyenv('Version', env_bin);
        %pause(20);  % Optional: Wait for 1 second?
        % Load the environment by executing a simple Python command
        py.exec('import sys');
    end
    % Display the environment details
    disp(pyenv);
    %-----------------------------------------------------------------   
    % Save main components to h5ad
    write_h5ad(sce);

    %-----------------------------------------------------------------
    % Execute python script
    python_executable = env_bin;  
    umap_wd = which('umap_sc_wrapper');
    umap_wd = erase(umap_wd,'umap_sc_wrapper.m');
    if ispc
        umap_wd = strrep(umap_wd,"\","\\");
    end
    python_script = strcat(umap_wd, 'umap_sc.py');

    % Call the Python script with the adjacency matrix file as argument
    % Base command
    system_command = sprintf('%s %s %s', python_executable, python_script, file_h5ad);
    
    % Append arguments
    args = sprintf('--output_csv %s --npca %d --ndim %d --my_res %.2f --mt_pct %.2f --min_genes0 %d --min_cells0 %d --max_counts %.2f --algo_cluster %s', ...
                   file_out, ...
                   npca, ...
                   ndim, ...
                   my_res, ...
                   mt_pct, ...
                   min_genes0, ...
                   min_cells0, ...
                   max_counts, ...
                   algo_cluster);
    
    % Conditionally add `use_hvg` argument
    if use_hvgs
        args = sprintf('%s --use_hvg', args);
    end

    % Combine base command and arguments
    system_command = sprintf('%s %s', system_command, args);
    
    % Execute commands
    [status, cmdout] = system(system_command);
    disp(cmdout)

    % Clean up the temporary adjacency matrix file
    if exist(file_h5ad, 'file')
        delete(file_h5ad);
    end    

    % Check for errors
    if status ~= 0
        disp('Error running the Python script:');
        disp(cmdout);
        return;
    else
        % Load the clustering results and umap
        mapping_clus = readtable(file_out);

        mapping_clus.cell_id = string(mapping_clus.cell_id);
        idx = ismember(sce.c_cell_id, mapping_clus.cell_id);
        sce = sce.selectcells(idx);
        % Display the clustering results
        clusters = mapping_clus.leiden;
        nclus = length( unique(clusters) );
        fprintf('Number Leiden clusters %d \n', nclus);
        sce.c_cluster_id = clusters + 1;

        % Overwritting UMAP
        if ndim == 2
            umap_s = [mapping_clus.UMAP_1, mapping_clus.UMAP_2];
            sce.struct_cell_embeddings.umap2d = umap_s;
        elseif ndim ==3
            umap_s = [mapping_clus.UMAP_1, mapping_clus.UMAP_2, mapping_clus.UMAP_3];
            sce.struct_cell_embeddings.umap3d = umap_s;
        end
        sce.s = umap_s;
        disp("Parsing umap function:")
        disp(cmdout);
    
    end

end
