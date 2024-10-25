num_iterations = 100;  % Number of iterations for resampling

% Initialize figure
figure;
hold on;

% Load the data and extract 'sce'
data = load("GSM4878207_GSM4878210.mat");
sce_tmp = data.sce;
clear data;

% Remove non-informative genes based on species
species = 'mouse';
sce_tmp = rm_non_inf_genes(sce_tmp, species);
batches = unique(sce_tmp.c_batch_id);

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
g0_dv = Tdv0.gene;

% Define the gene list 'gl' as the first 50 entries from Tdv0
gl = Tdv0.gene(1:50);

% Sub-sample 'Tdv0' with an exact match of the 'gl' values
Tdv0_sub = Tdv0(ismember(Tdv0.gene, gl), :);

clear idx_conf sce_tmp1 sce_tmp2;

% Loop for multiple iterations
for iter = 1:num_iterations
    % -------------------- RE-SAMPLING AND PROCESSING --------------------
    % Randomly group cells into new batches
    sce = randperm_batches2groups(sce_tmp, iter);

    % Get updated batch information
    batches = unique(sce.c_batch_id);
    nbatch = length(batches);

    % Calculate number of unique batch comparisons
    nupper = nbatch * (nbatch - 1) / 2 + 1;

    % Initialize a table to store results and metadata
    Tdv_sub = table('Size', [nupper, 2], ...
        'VariableTypes', {'cell', 'cell'}, ...
        'VariableNames', {'Result', 'BatchComparison'});

    Tdv_sub{1, "Result"} = {Tdv0_sub};
    Tdv_sub{1, "BatchComparison"} = {main_str};

    it = 2;
    for ibatch = 1:nbatch
        for jbatch = ibatch + 1:nbatch
            % Select cells belonging to each batch
            idx1 = sce.c_batch_id == batches(ibatch);
            sce1 = sce.selectcells(idx1);
            sce1 = sce1.qcfilter;

            idx2 = sce.c_batch_id == batches(jbatch);
            sce2 = sce.selectcells(idx2);
            sce2 = sce2.qcfilter;

            % Constructing the comparison name
            tmp_name = strcat(batches(ibatch), '-', batches(jbatch));

            % Call to differential expression analysis
            Tdv1 = sc_splinedv(sce1.X, sce2.X, sce1.g, sce2.g, [], false);

            % Storing result and meta-info in the table
            Tdv_sub{it, 'Result'} = {Tdv1( ismember(Tdv1.gene, gl), :) }; 
            tmp_name = regexprep(tmp_name, '-seed-\d+', '');  % Removes any pattern like '-seed-XXXX'

            Tdv_sub{it, 'BatchComparison'} = {tmp_name};
            it = it + 1;
        end
    end

    % Remove empty rows or rows with missing values
    Tdv_sub = Tdv_sub(~cellfun(@isempty, Tdv_sub.BatchComparison), :);

    % -------------------- PLOTTING --------------------
    % Define custom colors
    greenColor0 = [0, 0.8, 0];  % RGB for green (for rows with twice HFD)
    greenColor = [0, 1, 0];  % RGB for green (for rows with twice HFD)
    redColor = [0, 0, 0];      % RGB for red (for rows with twice LFD)

    % Loop through each row in Tdv_sub to determine color based on BatchComparison
    for i = 1:height(Tdv_sub)
        comparisonStr = Tdv_sub.BatchComparison{i};  % Extract the comparison string
        
        % Count occurrences of 'LFD' and 'HFD' in the comparison string
        countstr2 = count(comparisonStr, sub_str2);
        countstr1 = count(comparisonStr, sub_str1);
        
        % Assign colors based on the count of 'LFD' and 'HFD' 
        if countstr2 == 2 || countstr1 == 2
            colors(i, :) = redColor;   % Red if 'LFD' appears twice
        else
            colors(i, :) = greenColor; % Default color for other cases
        end
        if strcmp(comparisonStr, main_str)
            colors(i, :) = greenColor0; % Default color for HFD-LFD comparison
        end
    end

    % Extract and sort the gene labels based on the first table's distance differences
    firstTable = sortrows(Tdv_sub.Result{1}, 'dist_diff', 'descend');
    gl = firstTable.gene;  % Save the sorted gene labels based on the first table

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
        if iter == 1
            % Only assign DisplayName for the first, second, and third plots
            plot(x, reordered_dist_diff, '-o', 'DisplayName', Tdv_sub.BatchComparison{i}, ...
                'MarkerSize', 6, 'LineWidth', 1.2, 'Color', colors(i, :));
        else
            % Plot without DisplayName for other plots, and set 'HandleVisibility' to 'off' to exclude them from the legend
            plot(x, reordered_dist_diff, '-o', 'MarkerSize', 6, 'LineWidth', 1.2, ...
                'Color', colors(i, :), 'HandleVisibility', 'off');
        end
    end
end

% Format the plot
set(gca, 'XTick', x, 'XTickLabel', gl);
xtickangle(45);  % Rotate x-axis labels for better visibility
xlabel('Gene Labels');
ylabel('Distance Difference');
title('Gene Labels vs Distance Difference for All Tables');

% Show the legend only for the first, second, and third plots
legend('show', 'Location', 'northeastoutside');
grid on;
hold off;
