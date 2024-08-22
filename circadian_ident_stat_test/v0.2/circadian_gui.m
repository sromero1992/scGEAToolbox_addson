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

    % Create axes for plotting
    hPlotAxes = axes('Parent', hFig, 'Position', [0.4, 0.1, 0.55, 0.8]);

    % Callback function to define tmeta interactively and store it in GUI data
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

        % Create a new figure for tmeta definition
        tmetaFig = figure('Name', 'Define tmeta', 'Position', [150, 150, 500, 400]);

        % Initialize table data with old_labels and empty new_labels and times
        num_batches = numel(old_labels);
        initial_data = cell(num_batches, 2);
        initial_data(:, 1) = old_labels; % Set old_labels
        initial_data(:, 2) = num2cell(zeros(num_batches, 1)); % Set times with default values

        % Create a table UI to input new labels and times
        tmetaTable = uitable('Parent', tmetaFig, 'Position', [25, 75, 450, 250], ...
                             'Data', initial_data, ...
                             'ColumnName', {'Old Labels', 'Times'}, ...
                             'ColumnEditable', [false, true]);

        % Button to confirm and save tmeta
        uicontrol('Style', 'pushbutton', 'Position', [200, 25, 100, 30], ...
                  'String', 'Save tmeta', 'Callback', @saveTmeta);

        % Callback to save tmeta and close the window
        function saveTmeta(~, ~)
            % Retrieve data from the table
            tableData = get(tmetaTable, 'Data');
            old_labels = tableData(:, 1); % This should be cell array of strings
            times = cell2mat(tableData(:, 2)); % Convert times from cell array to numeric

            % Convert times to ZT labels
            new_labels = cell(numel(times), 1);
            for i = 1:numel(times)
                new_labels{i} = sprintf('ZT%02d', times(i)); % Format times as ZT00, ZT03, etc.
            end

            % Convert to table
            guiData.tmeta = table(old_labels, new_labels, times);

            % Close the tmeta definition window
            close(tmetaFig);
        end
    end

    % Callback function for Plot Gene button
    function plotGeneCallback(~, ~)
        % Retrieve parameters
        cust_cells = hCells.String{hCells.Value};
        period12 = hPeriod12.Value;
        gene = hGene.String{hGene.Value}; % Retrieve selected gene
    
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
        hWaitbar = waitbar(0, 'Plotting Gene...');
    
        % Suppress warnings related to directory creation
        warningState = warning('query', 'MATLAB:mkdir:DirectoryExists');
        warning('off', 'MATLAB:mkdir:DirectoryExists');
    
        try
            % Clear the existing plot in the axes
            cla(hPlotAxes);
    
            % Call the gene plotting function directly with the gene
            sce_circ_plot_gene(sce, guiData.tmeta, cust_cells, period12, gene, hPlotAxes);
    
            % Display completion message
            msgbox('Gene plotting completed.');
        catch ME
            % Display an error message if something goes wrong
            errordlg(['Error occurred: ', ME.message]);
        end
    
        % Close the waitbar
        close(hWaitbar);
    
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
        hWaitbar = waitbar(0, 'Plotting Genes...');

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
        close(hWaitbar);

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
        hWaitbar = waitbar(0, 'Analyzing...');

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
        close(hWaitbar);
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