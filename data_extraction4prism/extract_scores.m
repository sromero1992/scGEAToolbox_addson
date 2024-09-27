% Read the CSV file into a table
T = readtable('test_PRISM.csv');

% Get column names
column_names = T.Properties.VariableNames;

% Loop through each column and convert non-numeric columns to string if necessary
for i = 1:length(column_names)
    % Check if the current column is numeric
    if isnumeric(T.(column_names{i}))
        % Skip numeric columns
        continue;
    end
    % Convert non-numeric columns to string
    T.(column_names{i}) = string(T.(column_names{i}));
end

% Get all unique sample names and unique cell names
all_samples = unique(T.samples);
all_cells = unique(T.(column_names{1}));  % Assuming the first column contains the cell names

% Initialize a new table to store all the cells and scores for each sample
combined_table = table(all_cells, 'VariableNames', {'cells'});

% Loop through each unique sample and add scores to the combined table
for isample = 1:length(all_samples)
    % Extract rows matching the current sample
    sample_rows = all_samples(isample) == T.samples;
    
    % Extract scores for the current sample
    scores_column = T.(column_names{2});
    
    % Initialize a new column with -99 for all cells (default value)
    new_scores = -99 * ones(length(all_cells), 1);
    
    % Find the corresponding cell rows in the main table
    for icell = 1:length(all_cells)
        cell_rows = (all_cells(icell) == T.(column_names{1})) & sample_rows;
        if any(cell_rows)
            % If there's a matching cell for the sample, fill with actual score
            new_scores(icell) = scores_column(cell_rows);
        end
    end
    
    % Add the new scores column to the combined table with a sample-specific name
    new_column_name = ['scores_' char(all_samples(isample))];
    combined_table.(new_column_name) = new_scores;
end

% Write the combined table to a CSV file
writetable(combined_table, 'combined_scores_by_cells.csv');



% Read the CSV file into a table
T = readtable('test_PRISM.csv');

% Get column names
column_names = T.Properties.VariableNames;

% Loop through each column and convert non-numeric columns to string if necessary
for i = 1:length(column_names)
    % Check if the current column is numeric
    if isnumeric(T.(column_names{i}))
        % Skip numeric columns
        continue;
    end
    % Convert non-numeric columns to string
    T.(column_names{i}) = string(T.(column_names{i}));
end

% Get all unique sample names and unique cell names
all_samples = unique(T.samples);
all_cells = unique(T.(column_names{1}));  % Assuming the first column contains the cell names

% Initialize a new table to store all the cells and scores for each sample
compressed_table = table(all_cells, 'VariableNames', {'cells'});

% Loop through each unique sample and add scores to the compressed table
for isample = 1:length(all_samples)
    % Extract rows matching the current sample
    sample_rows = all_samples(isample) == T.samples;
    
    % Extract scores for the current sample
    scores_column = T.(column_names{2});
    
    % Initialize a new column with -99 for all cells (default value)
    new_scores = -99 * ones(length(all_cells), 1);
    
    % Find the corresponding cell rows in the main table
    for icell = 1:length(all_cells)
        cell_rows = (all_cells(icell) == T.(column_names{1})) & sample_rows;
        if any(cell_rows)
            % If there's a matching cell for the sample, fill with actual score
            new_scores(icell) = scores_column(cell_rows);
        end
    end
    
    % Add the new scores column to the compressed table with a sample-specific name
    new_column_name = ['scores_' char(all_samples(isample))];
    compressed_table.(new_column_name) = new_scores;
end

% Write the compressed table to a CSV file
writetable(compressed_table, 'compressed_scores_table.csv');

