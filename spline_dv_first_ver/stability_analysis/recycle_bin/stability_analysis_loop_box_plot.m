num_iterations = 2;  % Number of iterations for resampling

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

ngenes = 50;
% Define the gene list 'gl' as the first 50 entries from Tdv0
gl = Tdv0.gene(1:ngenes);

% Sub-sample 'Tdv0' with an exact match of the 'gl' values
Tdv0_sub = Tdv0(ismember(Tdv0.gene, gl), :);

clear idx_conf sce_tmp1 sce_tmp2;

% Initialize a results cell array to store data points
results = cell(num_iterations*ngenes, 2);  % Columns: {datapoints_green, datapoints_black}

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
            Tdv_sub{it, 'Result'} = {Tdv1(ismember(Tdv1.gene, gl), :)};  % Match genes exactly
            tmp_name = regexprep(tmp_name, '-seed-\d+', '');  % Removes any pattern like '-seed-XXXX'

            Tdv_sub{it, 'BatchComparison'} = {tmp_name};
            it = it + 1;
        end
    end

    % Remove empty rows or rows with missing values
    Tdv_sub = Tdv_sub(~cellfun(@isempty, Tdv_sub.BatchComparison), :);

    % Data containers for this iteration
    datapoints_green = []; % Container for HFD-LFD comparisons
    datapoints_black = []; % Container for other comparisons

    % Define a single green color for all green cases   
    % Loop through each row in Tdv_sub to determine color based on BatchComparison
    for i = 1:height(Tdv_sub)
        comparisonStr = Tdv_sub.BatchComparison{i};  % Extract the comparison string
    
        % Count occurrences of sub_str1 and sub_str2 in the comparison string
        countstr2 = count(comparisonStr, sub_str2);
        countstr1 = count(comparisonStr, sub_str1);
    
        % Check if the comparison string matches the main comparison (HFD-LFD)
        if strcmp(comparisonStr, main_str)
            % Store distance differences for the HFD-LFD comparison
            datapoints_green = [datapoints_green; Tdv_sub.Result{i}.dist_diff];
        elseif countstr2 == 2 || countstr1 == 2
            % Store distance differences for comparisons involving sub_str1 or sub_str2 twice
            datapoints_black = [datapoints_black; Tdv_sub.Result{i}.dist_diff];
        else
            % Store distance differences for all other comparisons
            datapoints_green = [datapoints_green; Tdv_sub.Result{i}.dist_diff];
        end
    end

    % Store the data points in the results cell array
    results{iter, 1} = datapoints_green;  % Store HFD-LFD comparison data
    results{iter, 2} = datapoints_black;  % Store other comparison data
end


% Save results to a .mat file (optional)
% save('results.mat', 'results', 'gl');

% -------------------- PLOTTING --------------------
% Gather all data points across iterations
all_datapoints_green = [];
all_datapoints_black = [];

for iter = 1:num_iterations
    all_datapoints_green = [all_datapoints_green; results{iter, 1}];
    all_datapoints_black = [all_datapoints_black; results{iter, 2}];
end

% Combine the green and black data points into a single column for boxplot
all_data = [all_datapoints_green; all_datapoints_black];

% Create a grouping variable to distinguish between the two categories
grouping = [repmat({'HFD-LFD'}, length(all_datapoints_green), 1); ...
            repmat({'Other Comparisons'}, length(all_datapoints_black), 1)];

% Create the boxplot using the combined data and grouping
figure;
boxplot(all_data, grouping, 'Labels', {'HFD-LFD', 'Other Comparisons'}, ...
    'Colors', [0, 0.8, 0; 0, 0, 0]);

title('Boxplot of Distance Differences');
ylabel('Distance Difference');
grid on;

