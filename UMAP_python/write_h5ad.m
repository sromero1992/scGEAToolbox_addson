function write_h5ad(sce, do_scratch_files)

    if nargin < 2; do_scratch_files=true; end
    % Load python environment
    python_executable = init_python_env_matlab();

    if do_scratch_files
        % Save the primary data matrices to text files
        [i, j, val] = find(sce.X);
        writematrix([i, j, val], 'sce_X.csv', 'Delimiter', 'tab'); % Gene expression matrix
        writematrix(sce.g, 'sce_genes.csv');                       % Gene names
        writematrix(sce.s, 'sce_embeddings.csv');                  % Embedding (UMAP, etc.)
        writematrix(sce.c, 'sce_clusters.csv');                    % Cluster IDs
        writematrix(sce.c_cell_type_tx, 'sce_celltypes.csv');       % Cell type annotations
        writematrix(sce.c_batch_id, 'sce_batch.csv');              % Batch IDs
        writematrix(sce.c_cell_id, 'sce_cell_ids.csv');            % Cell IDs
           
        % Get the total number of items in the cell attribute list
        ncell_att = length(sce.list_cell_attributes);
    
        % The list contains pairs of {name, data}, so we loop in steps of 2
        % This loop will process all the additional cell metadata attributes
        if ncell_att > 0
            fprintf('Processing %d additional cell attributes...\n', ncell_att/2);
            for idx = 1:2:ncell_att
                % The attribute name (e.g., 'SubCellType') is at the current index
                tag_name = sce.list_cell_attributes{idx};
                
                % The data for that attribute is at the next index
                data_to_write = sce.list_cell_attributes{idx + 1};
                
                % === This is the existence check you asked for ===
                % Check if the tag_name is a valid string and the data is not empty
                if (ischar(tag_name) || isstring(tag_name)) && ~isempty(tag_name) && ~isempty(data_to_write)
                    
                    % Construct the filename (e.g., 'SubCellType_sce_metadata.csv')
                    % Using [] for string concatenation is often more robust for filenames
                    filename = [tag_name, '_sce_metadata.csv'];
                    
                    % Write the data to the corresponding file
                    fprintf('Writing metadata to: %s\n', filename);
                    writematrix(data_to_write, filename);
                    
                else
                    % Optional: Warn the user if a pair is malformed or empty
                    fprintf('Skipping attribute at index %d: Invalid name or empty data.\n', idx);
                end
            end
        end
    
    end
        
    % Execute python script
    write_h5ad_wd = fileparts(which('write_h5ad')); % A more robust way to get the directory

    if ispc
        % On Windows, python may not need the escaped backslashes if called via system()
        % but it's safer to use fullfile() to construct paths.
        python_script = fullfile(write_h5ad_wd, 'construct_h5ad.py');
    else
        python_script = fullfile(write_h5ad_wd, 'construct_h5ad.py');
    end

    system_command = sprintf('"%s" "%s"', python_executable, python_script);
    fprintf('Executing command: %s\n', system_command);
    [status, cmdout] = system(system_command);
    
    if status == 0
        disp('Python script executed successfully.');
        disp(cmdout);
    else
        error('Error executing Python script:\n%s', cmdout);
    end
end