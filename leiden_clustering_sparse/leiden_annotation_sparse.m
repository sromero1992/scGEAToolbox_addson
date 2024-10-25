function sce = leiden_annotation_sparse(sce, species, method)
    % leiden_annotation computes leiden clustering interfaced from
    % python3.11 with Mutual Nearest Neighbors (MNN) or K-Nearest
    % Neighbors (KNN).
    % INPUT:
    % sce -------> SCE object 
    % method ----> Method  to find neighbors (mnn or knn)
    % OUTPUT:
    % sce ------> sce object containing Leiden clusters and corresponding
    %             annotation
    % Usage:
    % sce = leiden_annotation_sparse(sce,'knn','mouse')
    % 
    % If no annotation wanted, then use
    % sce = leiden_annotation_sparse(sce,'knn', [])
    if nargin < 2; species = []; end
    if nargin < 3; method = 'knn'; end

    species = lower(species);
    method = lower(method);

    fprintf("WARNING: sce object should be prior QC if desired... \n")
    fprintf("Leiden annotation with method: %s \n", method);

    % Set the Python environment (Python 3.11)
    % Windows format
    env_bin = 'C:\Users\ssromerogon\.conda\envs\leiden_clustering\python.exe';
    if ispc
        env_bin = strrep(env_bin,"\","\\");
    end
    % Linux format
    %env_bin = "/home/ssromerogon/packages/scanpy_env/bin/python3";

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

    switch method
        case 'mnn'
            n_neighbors = 100; % 20, 30, 50 (the more neighbors, the less clusters)
            adjX = adj_mat_construct_sparse(sce, 'mnn', n_neighbors);
    
        case 'knn'
            n_neighbors = 15; % 10, 15, 20 % the less the neighbors produces more clusters
            adjX = adj_mat_construct_sparse(sce, 'knn', n_neighbors);

        otherwise
            error('Unknown method: %s. Method should be either ''mnn'' or ''knn''.', method);
    end

    disp(size(adjX));

    % Save the adjacency matrix to a text file
    adj_file = 'adjX.txt';
    [i, j, val] = find(adjX);
    writematrix([i, j, val], adj_file, 'Delimiter', 'tab');
    clear adjX i j val;
    
    % Path to the Python executable and the script
    python_executable = env_bin;  
    leiden_wd = which('leiden_annotation_sparse');
    leiden_wd = erase(leiden_wd,'leiden_annotation_sparse.m');
    if ispc
        leiden_wd = strrep(leiden_wd,"\","\\");
    end
    python_script = strcat(leiden_wd, 'run_leiden_sparse.py');

    % Call the Python script with the adjacency matrix file as argument
    system_command = sprintf('%s %s %s', python_executable, python_script, adj_file );
    [status, cmdout] = system(system_command);

    % Clean up the temporary adjacency matrix file
    if exist(adj_file, 'file')
        delete(adj_file);
    end    

    % Check for errors
    if status ~= 0
        disp('Error running the Python script:');
        disp(cmdout);
        return;
    else
        % Load the clustering results
        clusters = jsondecode(fileread('clusters.json'));
    
        delete('clusters.json')
    
        disp("Parsing Leiden:")
        disp(cmdout);
    
        % Display the clustering results
        nclus = length( unique(clusters) );
        fprintf('Number Leiden clusters %d \n', nclus);
    end

    % Assign clustering results to the SCE object
    sce.c_cluster_id = clusters + 1;

    % Embed cells and assign cell types
    tic;
    sce = sce.embedcells('umap2d', true, false, 2);
    %rng('default');
    %sce = sce.embedcells('tsne3d', true, false, 3);
    if ~isempty(species)
        fprintf("Annotating species %s \n\n", species)
        sce = sce.assigncelltype(species, false);
    end

    time_assign = toc;
    fprintf("Time for cell annotation and embedding: %f \n", time_assign);
end
