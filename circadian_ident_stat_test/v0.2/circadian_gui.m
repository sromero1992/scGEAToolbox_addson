function circadian_gui
    % Create the main figure window
    hFig = figure('Name', 'SCE-Circadian Analysis GUI', 'Position', [100, 100, 1000, 600]);

    % Define your custom plot types
    plot_types = {'Plot confident genes ', ...
                  'Plot non-condifent genes', ...
                  'Plot circadian only', ...
                  'Plot syncronized genes'};

    % Load the sce.mat file
    sce = load_sce_data();

    % Initialize GUI data for tmeta
    guiData = struct('tmeta', [], 'cust_cells', '', 'plot_type', 1, 'gene', '');

    % Create button to define tmeta interactively
    uicontrol('Style', 'pushbutton', 'Position', [50, 500, 150, 30], ...
              'String', 'Define Tmeta', 'Callback', @defineTmetaAndPlotCallback);

    % Create dropdown for selecting cell type
    cell_types = unique(sce.c_cell_type_tx); % Ensure this matches the field in the 'sce' structure
    
    % Create UI controls
    uicontrol('Style', 'text', 'Position', [35, 450, 100, 20], 'String', 'Cell Type:');
    hCells = uicontrol('Style', 'popupmenu', 'Position', [120, 450, 150, 20], ...
                       'String', cell_types, 'Value', 1); % Set default selection to the "None" option
   
    % Create dropdown for selecting plot_type
    uicontrol('Style', 'text', 'Position', [35, 400, 100, 20], 'String', 'Plot Type:');
    hPlotType = uicontrol('Style', 'popupmenu', 'Position', [120, 400, 150, 20], ...
                          'String', plot_types);

    % Create checkbox for selecting period12
    hPeriod12 = uicontrol('Style', 'checkbox', 'Position', [70, 350, 150, 20], ...
                          'String', 'Period 12 (otherwise 24)');


    % Create button for Analysis
    uicontrol('Style', 'pushbutton', 'Position', [50, 300, 150, 30], ...
              'String', 'Analyze (Faster)', 'Callback', @analyzeButtonCallback);

    % Create button for Plot All
    uicontrol('Style', 'pushbutton', 'Position', [50, 250, 150, 30], ...
              'String', 'Plot Genes & Analyze', 'Callback', @plotAllCallback);

    
    % Create text box for gene input 
    genes = sce.g; % Assuming sce.g is a cell array or similar
    uicontrol('Style', 'text', 'Position', [45, 150, 100, 20], 'String', 'Single gene plot:');
    hGene = uicontrol('Style', 'popupmenu', 'Position', [140, 150, 150, 20], ...
                       'String', genes); % Dropdown menu for genes


    % Create button for Plot Gene
    uicontrol('Style', 'pushbutton', 'Position', [50, 100, 150, 30], ...
              'String', 'Plot single Gene', 'Callback', @plotGeneCallback);

    % Create checkbox for printing sc-data
    hPrintSCdata = uicontrol('Style', 'checkbox', 'Position', [70, 70, 150, 20], ...
                          'String', 'Print sc-data (single gene)');

    % Create axes for plotting
    hPlotAxes = axes('Parent', hFig, 'Position', [0.4, 0.1, 0.55, 0.8]);

    function defineTmetaAndPlotCallback(~, ~)
        % Get unique batches from sce.c_batch_id
        batches = unique(sce.c_batch_id);
    
        % Check the format of batches and convert accordingly
        if isnumeric(batches)
            old_labels = cellstr(num2str(batches)); % Convert numeric batches to strings
        elseif iscellstr(batches)
            old_labels = batches; % Already in cell array of strings
        elseif isstring(batches)
            old_labels = cellstr(batches); % Convert string array to cell array of character vectors
        else
            error('Unexpected format for batches.');
        end
    
        % Calculate the number of cells for each batch
        num_batches = numel(old_labels);
        cell_counts = zeros(num_batches, 1); % Initialize cell counts array
        for i = 1:num_batches
            cell_counts(i) = sum(strcmp(sce.c_batch_id, old_labels{i}));
        end
    
        % Create a new figure for tmeta definition
        tmetaFig = figure('Name', 'Define tmeta', 'Position', [150, 150, 550, 400]);
    
        % Initialize table data with old_labels, times, and cell_counts
        initial_data = cell(num_batches, 3);
        initial_data(:, 1) = old_labels; % Set old_labels
        initial_data(:, 2) = num2cell(zeros(num_batches, 1)); % Set times with default values
        initial_data(:, 3) = num2cell(cell_counts); % Set cell counts for each batch
    
        % Create a table UI to input new labels and times, with cell counts as non-editable
        tmetaTable = uitable('Parent', tmetaFig, 'Position', [25, 75, 500, 250], ...
                             'Data', initial_data, ...
                             'ColumnName', {'Old Labels', 'Times', 'Cell Count'}, ...
                             'ColumnEditable', [false, true, false]);  % Cell count column is non-editable
    
        % Button to save tmeta changes and update SCE
        uicontrol('Style', 'pushbutton', 'Position', [60, 25, 100, 30], ...
                  'String', 'Save tmeta', 'Callback', @saveTmeta);
    
        % Button to save the modified SCE object
        uicontrol('Style', 'pushbutton', 'Position', [190, 25, 100, 30], ...
                  'String', 'Save SCE', 'Callback', @saveNewSCE);
    
        % Button to close the tmeta definition window
        uicontrol('Style', 'pushbutton', 'Position', [320, 25, 100, 30], ...
                  'String', 'Close', 'Callback', @(~, ~) close(tmetaFig));
    
        % Callback to save tmeta and update the SCE object
        function saveTmeta(~, ~)
            % Retrieve the updated table data
            tableData = get(tmetaTable, 'Data');
            old_labels = tableData(:, 1); % This should be cell array of strings
            times = cell2mat(tableData(:, 2)); % Convert times from cell array to numeric
            cell_counts = cell2mat(tableData(:, 3)); % Cell counts remain unchanged
    
            % Convert times to ZT labels
            new_labels = cell(numel(times), 1);
            for i = 1:numel(times)
                if times(i) >= 0 
                    new_labels{i} = sprintf('ZT%02d', times(i)); 
                else 
                    new_labels{i} = sprintf('%d', times(i)); 
                end
            end
    
            % Convert to table and sort times (if not sorted)
            tbl = table(old_labels, new_labels, times);
            tbl = sortrows(tbl, 'times');
            guiData.tmeta = tbl; % Update GUI data
    
            % Process to remove selected batches
            rm_batch = ismember(guiData.tmeta.new_labels, '-1');
            rm_batch = guiData.tmeta.old_labels(rm_batch);
        
            for ib = 1:length(rm_batch)
                % Use pre-selected batches to remove
                idx_rm = find(strcmp(sce.c_batch_id, rm_batch{ib}));
                % Rename batch 
                sce.c_batch_id(idx_rm) = "-1";
            end
    
            % Subsetting sce dynamically
            idx = ~strcmp(sce.c_batch_id, "-1");
            sce_new = SingleCellExperiment(sce.X(:, idx), sce.g);
            sce_new.c_cell_type_tx = sce.c_cell_type_tx(idx);
            sce_new.c_cell_id = sce.c_cell_id(idx);
            sce_new.c_batch_id = sce.c_batch_id(idx);
    
            % Merge the new labels that are the same
            [unique_labels, ~, label_idx] = unique(guiData.tmeta.new_labels);
            for i = 1:length(unique_labels)
                current_label = unique_labels{i};
                % Find indices in the new SCE that match this new label
                merge_idx = ismember(sce_new.c_batch_id, guiData.tmeta.old_labels(label_idx == i));
                % Update batch ID with the merged label
                sce_new.c_batch_id(merge_idx) = current_label;
            end
    
            % Update the tmeta with the new old_labels (now as new_labels)
            unique_new_labels = unique(guiData.tmeta.new_labels);  % Get unique new labels after merging
            updated_times = zeros(length(unique_new_labels), 1); % Initialize updated times
            updated_counts = zeros(length(unique_new_labels), 1); % Initialize updated cell counts
    
            for i = 1:length(unique_new_labels)
                % Get index of the old label corresponding to the unique new label
                idx_old_label = find(strcmp(guiData.tmeta.new_labels, unique_new_labels{i}), 1);
                updated_times(i) = guiData.tmeta.times(idx_old_label); % Assign corresponding times
                updated_counts(i) = cell_counts(idx_old_label); % Keep cell counts consistent
            end
            
            updatedData = [unique_new_labels, num2cell(updated_times), num2cell(updated_counts)];
            set(tmetaTable, 'Data', updatedData);  % Update the data in the table
    
            % Replace the original sce with the new sce
            clear sce; % Clear the old SCE object to save memory
            sce = sce_new; % Assign the new SCE object to 'sce'
            clear sce_new; % Clear the temporary new SCE object
    
            % Update the GUI with the new SCE object
            assignin('base', 'sce', sce);
    
            % Notify user of successful update
            disp('tmeta saved and SCE updated.');
        end
    
        % Callback to save the new SCE object
        function saveNewSCE(~, ~)
            % Prompt user to select a file to save the new SCE object
            [file, path] = uiputfile('*.mat', 'Save Modified SCE As');
            if isequal(file, 0)
                disp('User canceled the file save.');
            else
                save(fullfile(path, file), 'sce');
                disp(['Modified SCE saved to ', fullfile(path, file)]);
            end
        end
    end


    % Callback function for Plot Gene button
    function plotGeneCallback(~, ~)
        % Retrieve parameters
        cust_cells = hCells.String{hCells.Value};
        period12 = hPeriod12.Value;
        gene = hGene.String{hGene.Value}; % Retrieve selected gene
        print_scdata = hPrintSCdata.Value;
        % Check if tmeta is defined
        if isempty(guiData.tmeta)
            errordlg('Please define tmeta before plotting.');
            return;
        end
    
        % Check if gene is empty or not in the dataset
        if isempty(gene) || ~ismember(gene, sce.g)
            errordlg('Selected gene is invalid or not found in the dataset.');
            return;
        end
    
        % Create and display waitbar
        %hWaitbar = waitbar(0, 'Plotting Gene...');
        disp("Working on Plot Gene...")

        % Suppress warnings related to directory creation
        warningState = warning('query', 'MATLAB:mkdir:DirectoryExists');
        warning('off', 'MATLAB:mkdir:DirectoryExists');
    
        try
            % Clear the existing plot in the axes
            cla(hPlotAxes);
    
            % Call the gene plotting function directly with the gene
            sce_circ_plot_gene(sce, guiData.tmeta, cust_cells, period12, gene, hPlotAxes, print_scdata);
    
            % Display completion message
            msgbox('Gene plotting completed.');
        catch ME
            % Display an error message if something goes wrong
            errordlg(['Error occurred: ', ME.message]);
        end
    
        % Close the waitbar
        %close(hWaitbar);
        disp("Finished Plot Gene")

    
        % Restore the warning state
        warning(warningState);
    end

    % Callback function for Plot All button
    function plotAllCallback(~, ~)
        % Retrieve parameters
        cust_cells = hCells.String{hCells.Value};
        %if contains(cust_cells,"None"); cust_cells = []; end
        plot_type = hPlotType.Value;
        period12 = hPeriod12.Value;

        % Check if tmeta is defined
        if isempty(guiData.tmeta)
            errordlg('Please define tmeta before plotting.');
            return;
        end

        % Create and display waitbar
        %hWaitbar = waitbar(0, 'Plotting Genes...');
        disp("Working on Plot Genes...")

        % Suppress warnings related to directory creation
        warningState = warning('query', 'MATLAB:mkdir:DirectoryExists');
        warning('off', 'MATLAB:mkdir:DirectoryExists');

        try
            % Call the sce_circ_plot function for all genes
            sce_circ_plot(sce, guiData.tmeta, cust_cells, plot_type, period12);
             
            % Display completion message
            msgbox('All plotting completed.');
        catch ME
            % Display an error message if something goes wrong
            errordlg(['Error occurred: ', ME.message]);
        end

        % Close the waitbar
        %close(hWaitbar);
        disp("Finished Plot Genes!")


        % Restore the warning state
        warning(warningState);
    end

    % Callback function for Analysis button
    function analyzeButtonCallback(~, ~)
        % Retrieve parameters
        cust_cells = hCells.String{hCells.Value};
        period12 = hPeriod12.Value;

        % Check if tmeta is defined
        if isempty(guiData.tmeta)
            errordlg('Please define tmeta before analyzing.');
            return;
        end

        % Create and display waitbar
        %hWaitbar = waitbar(0, 'Analyzing...');
        disp("Working on Plot Genes...")


        try
            % Call the sce_circ_phase_estimation function      
            [T1, T2] = sce_circ_phase_estimation_stattest(sce, guiData.tmeta, ...
                                     true, period12, [], cust_cells);   
            % Display analysis results in a message box or other UI elements
            msgbox('Analysis completed.');

            % Additional code to display T1 and T2 results can be added here
        catch ME
            % Display an error message if something goes wrong
            errordlg(['Error occurred: ', ME.message]);
        end

        % Close the waitbar
        %close(hWaitbar);
        disp("Finished Plot Genes...")

    end

    % Function to load sce data (customize as needed)
    function sce = load_sce_data()
        [fileName, filePath] = uigetfile('*.mat', 'Select SCE Data File');
        if isequal(fileName, 0)
            error('No file selected.');
        end
        loadedData = load(fullfile(filePath, fileName));
        disp("Successfully loaded sce...");
    
        % Display the fields of the loaded data
        disp('Fields in loaded data:');
        disp(fieldnames(loadedData));
    
        % Assuming 'sce' is a field in the loaded data
        sce = loadedData.sce;
    
        % Display fields of the 'sce' structure
        disp('Fields in the sce structure:');
        disp(fieldnames(sce));
    end
end