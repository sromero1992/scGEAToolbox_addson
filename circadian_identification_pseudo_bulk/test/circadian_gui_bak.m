function circadian_gui_bak
    % Create the main figure window
    hFig = figure('Name', 'Circadian Analysis GUI', 'Position', [100, 100, 1000, 600]);

    % Define your custom plot types
    plot_types = {'Plot failed (but confident circadian gene)', ...
                  'Plot all (confident circadian gene)'};

    % Load the sce.mat file
    sce = load_sce_data();

    % Initialize GUI data for tmeta
    guiData = struct('tmeta', [], 'cust_cells', '', 'plot_type', 1);

    % Populate dropdown for selecting cust_cells based on sce object
    cell_types = unique(sce.c_cell_type_tx); % Assuming sce.CellTypes is a cell array or similar
    uicontrol('Style', 'text', 'Position', [50, 350, 100, 20], 'String', 'Cell Type:');
    hCells = uicontrol('Style', 'popupmenu', 'Position', [150, 350, 150, 20], ...
                       'String', cell_types);

    % Create text box for gene input
    uicontrol('Style', 'text', 'Position', [50, 300, 100, 20], 'String', 'Gene:');
    hGene = uicontrol('Style', 'edit', 'Position', [150, 300, 150, 20]);

    % Create dropdown for selecting plot_type
    uicontrol('Style', 'text', 'Position', [50, 250, 100, 20], 'String', 'Plot Type:');
    hPlotType = uicontrol('Style', 'popupmenu', 'Position', [150, 250, 150, 20], ...
                          'String', plot_types);

    % Create checkbox for selecting period12
    hPeriod12 = uicontrol('Style', 'checkbox', 'Position', [50, 200, 150, 20], ...
                          'String', 'Period 12 (otherwise 24)');

    % Create axes for plotting
    hAxes = axes('Parent', hFig, 'Position', [0.4, 0.1, 0.55, 0.8]);

    % Create button to define tmeta interactively
    uicontrol('Style', 'pushbutton', 'Position', [50, 150, 250, 30], ...
              'String', 'Define Tmeta', 'Callback', @defineTmetaAndPlotCallback);

    % Create button for Plot (all genes)
    uicontrol('Style', 'pushbutton', 'Position', [50, 100, 100, 30], ...
              'String', 'Plot All', 'Callback', @plotButtonCallback);

    % Create button for Plot (individual gene)
    uicontrol('Style', 'pushbutton', 'Position', [50, 50, 150, 30], ...
              'String', 'Plot Gene', 'Callback', @plotGeneButtonCallback);

    % Create button for Analysis
    uicontrol('Style', 'pushbutton', 'Position', [200, 50, 100, 30], ...
              'String', 'Analyze', 'Callback', @analyzeButtonCallback);

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

    % Callback function for Plot (all genes) button
    function plotButtonCallback(~, ~)
        % Retrieve parameters
        cust_cells = hCells.String{hCells.Value};
        plot_type = hPlotType.Value;
        period12 = hPeriod12.Value;

        % Check if tmeta is defined
        if isempty(guiData.tmeta)
            errordlg('Please define tmeta before plotting.');
            return;
        end
    
        % Suppress warnings related to directory creation
        warningState = warning('query', 'MATLAB:mkdir:DirectoryExists');
        warning('off', 'MATLAB:mkdir:DirectoryExists');
    
        try
            % Call the sce_circ_plot function without gene parameter
            sce_circ_plot(sce, guiData.tmeta, cust_cells, plot_type, period12);
            
            % Display some kind of result or message in the GUI
            msgbox('Plotting completed.');
        catch ME
            % Display an error message if something goes wrong
            errordlg(['Error occurred: ', ME.message]);
        end
    
        % Restore the warning state
        warning(warningState);
    end

    % Callback function for Plot (individual gene) button
    function plotGeneButtonCallback(~, ~)
        % Retrieve parameters
        cust_cells = hCells.String{hCells.Value};
        plot_type = hPlotType.Value;
        period12 = hPeriod12.Value;
        cust_gene = get(hGene, 'String'); % Retrieve selected gene
    
        % Check if tmeta is defined
        if isempty(guiData.tmeta)
            errordlg('Please define tmeta before plotting.');
            return;
        end
    
        % Suppress warnings related to directory creation
        warningState = warning('query', 'MATLAB:mkdir:DirectoryExists');
        warning('off', 'MATLAB:mkdir:DirectoryExists');
    
        try
            % Call the sce_circ_plot function with the gene parameter
            sce_circ_plot(sce, guiData.tmeta, cust_cells, plot_type, period12, cust_gene);
            
            % Display some kind of result or message in the GUI
            msgbox('Gene plotting completed.');
        catch ME
            % Display an error message if something goes wrong
            errordlg(['Error occurred: ', ME.message]);
        end
    
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

        % Call the sce_circ_phase_estimation function
        [T1, T2] = sce_circ_phase_estimation(sce, guiData.tmeta, true, period12, [], cust_cells);

        % Display analysis results in a message box or table
        msgbox('Analysis completed.');
    end

end

function sce = load_sce_data()
    % Load sce.mat file
    [file, path] = uigetfile('*.mat', 'Select the SCE MAT-file');
    if isequal(file, 0)
        error('No file selected. Exiting...');
    else
        sce_data = load(fullfile(path, file));
        sce = sce_data.sce; % Assuming the file contains a variable 'sce'
    end
end