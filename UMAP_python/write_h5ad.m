function write_h5ad(sce)    

    % Save the adjacency matrix to a text file
    [i, j, val] = find(sce.X);
    writematrix([i, j, val], 'sce_X.csv', 'Delimiter', 'tab');% Gene expression matrix
    writematrix(sce.g, 'sce_genes.csv');             % Gene names
    writematrix(sce.s, 'sce_embeddings.csv');        % Embedding (UMAP, etc.)
    writematrix(sce.c, 'sce_clusters.csv');          % Cluster IDs
    writematrix(sce.c_cell_type_tx, 'sce_celltypes.csv'); % Cell type annotations
    writematrix(sce.c_batch_id, 'sce_batch.csv');    % Batch IDs
    writematrix(sce.c_cell_id, 'sce_cell_ids.csv');  % Cell IDs

    %-----------------------------------------------------------------
    % Set the Python environment (Python 3.11)
    % Windows format
    env_bin = 'F:\Anaconda\envs\scanpy_env_311\python.exe';
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

    % Execute python script
    python_executable = env_bin;  
    %write_h5ad_wd = which('write_h5ad');
    %write_h5ad_wd = erase(write_h5ad_wd,'write_h5ad.m');
    if ispc
        write_h5ad_wd = strrep(write_h5ad_wd,"\","\\");
    end
    python_script = strcat(write_h5ad_wd, 'construct_h5ad.py');

    system_command = sprintf('%s %s', python_executable, python_script);
    [status, cmdout] = system(system_command);
    disp(cmdout)

end