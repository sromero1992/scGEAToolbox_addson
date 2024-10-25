% Load the data and extract 'sce'
data = load("GSM4878207_GSM4878210.mat");
sce_tmp = data.sce;
clear data;


% Remove non-informative genes based on species
species = 'mouse';
sce_tmp = rm_non_inf_genes(sce_tmp, species);
batches = unique(sce_tmp.c_batch_id);

% Labels for plots
sub_str1 = batches(1);
sub_str2 = batches(2);
main_str = strcat(sub_str1, "-", sub_str2);

% Select cells belonging to each batch
idx1 = sce_tmp.c_batch_id == batches(1);
sce_tmp1 = sce_tmp.selectcells(idx1);
sce_tmp1 = sce_tmp1.qcfilter;

idx2 = sce_tmp.c_batch_id == batches(2);
sce_tmp2 = sce_tmp.selectcells(idx2);
sce_tmp2 = sce_tmp2.qcfilter;

% Call to differential variability analysis base
Tdv0 = sc_splinedv(sce_tmp1.X, sce_tmp2.X, sce_tmp1.g, sce_tmp2.g, [], false);
%idx_conf = comparison_result.pval <= 0.05;
%comparison_result = comparison_result(idx_conf, :);
g0_dv = Tdv0.gene;

% Define the gene list 'gl' as the first 50 entries from Tdv0
gl = Tdv0.gene(1:50);

% Sub-sample 'Tdv0' with an exact match of the 'gl' values
Tdv0_sub = Tdv0(ismember(Tdv0.gene, gl), :);

clear idx_conf sce_tmp1 sce_tmp2;

% Randomly group cells into new batches
sce = randperm_batches2groups(sce_tmp);

% Get updated batch information
batches = unique(sce.c_batch_id);
nbatch = length(batches);

% Calculate number of unique batch comparisons
nupper = nbatch * (nbatch - 1) / 2 + 1;

% Initialize a table to store results and metadata
Tdv = table('Size', [nupper, 2], ...
    'VariableTypes', {'cell', 'cell'}, ...
    'VariableNames', {'Result', 'BatchComparison'});

% Initialize a table to store results and metadata
Tdv_sub = table('Size', [nupper, 2], ...
    'VariableTypes', {'cell', 'cell'}, ...
    'VariableNames', {'Result', 'BatchComparison'});

Tdv{1, "Result"} = {Tdv0};
str_main = strcat(sub_str1,"-",sub_str2);
Tdv{1, "BatchComparison"} = {str_main};

Tdv_sub{1, "Result"} = {Tdv0_sub};
Tdv_sub{1, "BatchComparison"} = {str_main};

it = 2;
for ibatch = 1:nbatch
    % ilfdbool = contains(batches(ibatch), 'LFD');
    % ihfdbool = contains(batches(ibatch), 'HFD');

    for jbatch = ibatch + 1:nbatch
        % 
        % jlfdbool = contains(batches(jbatch), 'LFD');
        % jhfdbool = contains(batches(jbatch), 'HFD');

        % if (ilfdbool && jlfdbool) || (ihfdbool && jhfdbool)
        %     continue;
        % end

        % Select cells belonging to each batch
        idx1 = sce.c_batch_id == batches(ibatch);
        sce1 = sce.selectcells(idx1);
        sce1 = sce1.qcfilter;

        idx2 = sce.c_batch_id == batches(jbatch);
        sce2 = sce.selectcells(idx2);
        sce2 = sce2.qcfilter;

        % Constructing the comparison name
        tmp_name = strcat(batches(ibatch), '-', batches(jbatch));

        % Call to differential expression analysis (replace with actual function)
        Tdv1 = sc_splinedv(sce1.X, sce2.X, sce1.g, sce2.g, tmp_name, false);
        %idx_conf = comparison_result.pval <= 0.05;
        %comparison_result = comparison_result(idx_conf, :);

        % Storing result and meta-info in the table
        Tdv{it, 'Result'} = {Tdv1};
        tmp_name = regexprep(tmp_name, '-seed-\d+', '');  % Removes any pattern like '-seed-XXXX'
        Tdv{it, 'BatchComparison'} = {tmp_name};

        Tdv_sub{it, 'Result'} = {Tdv1( ismember(Tdv1.gene, gl), :) }; 
        Tdv_sub{it, 'BatchComparison'} = {tmp_name};

        it = it + 1;
    end
end

% Remove empty rows or rows with missing values
Tdv = Tdv(~cellfun(@isempty, Tdv.BatchComparison), :);
Tdv_sub = Tdv_sub(~cellfun(@isempty, Tdv_sub.BatchComparison), :);




% Define custom colors
greenColor0 = [0, 1, 0];  % RGB for green (for rows with twice HFD)
greenColor = [0, 0.7, 0];  % RGB for green (for rows with twice HFD)
blackColor = [0, 0, 0];      % RGB for red (for rows with twice LFD)

% Loop through each row in Tdv_sub to determine color based on BatchComparison
for i = 1:height(Tdv_sub)
    comparisonStr = Tdv_sub.BatchComparison{i};  % Extract the comparison string
    
    % Count occurrences of 'LFD' and 'HFD' in the comparison string
    count_str2 = count(comparisonStr, sub_str2);
    count_str1 = count(comparisonStr, sub_str1);
    
    % Assign colors based on the count of 'LFD' and 'HFD' 
    if count_str2 == 2 || count_str1 == 2
        % This is considered noise
        colors(i, :) = blackColor;   % Red if 'LFD' appears twice
    else
        % This is meaninful
        colors(i, :) = greenColor; % Default color for other cases
    end
    if strcmp(comparisonStr, main_str)
        colors(i, :) = greenColor0; % Default color for other cases
    end

end

% Extract and sort the gene labels based on the first table's distance differences
firstTable = sortrows(Tdv_sub.Result{1}, 'dist_diff', 'descend');
gl = firstTable.gene;  % Save the sorted gene labels based on the first table

% Create a new figure
figure;
hold on;

% Loop through each table in Tdv_sub for plotting
for i = 1:height(Tdv_sub)
    % Extract the current table
    currentTable = Tdv_sub.Result{i};
    
    % Reorder the current table based on the sorted genes in 'gl'
    [~, idx] = ismember(gl, currentTable.gene);  % Find the indices of sorted genes in the current table

    % Create a distance difference array with NaNs for missing genes
    reordered_dist_diff = nan(size(gl));  % Initialize with NaN
    validIndices = idx > 0;               % Identify valid indices in 'idx'
    reordered_dist_diff(validIndices) = currentTable.dist_diff(idx(validIndices));  % Assign distance differences to valid indices

    % X-axis indices based on the static gene order
    x = 1:length(gl);
    
    % Plot the current table's distance differences with circle markers
    plot(x, reordered_dist_diff, '-o', 'DisplayName', Tdv_sub.BatchComparison{i}, ...
        'MarkerSize', 6, 'LineWidth', 1.2, 'Color', colors(i, :));
end


% Format the plot
set(gca, 'XTick', x, 'XTickLabel', gl);
xtickangle(45);  % Rotate x-axis labels for better visibility
xlabel('Gene Labels');
ylabel('Distance Difference');
title('Gene Labels vs Distance Difference for All Tables');
legend('show', 'Location', 'northeastoutside');
grid on;
hold off;

