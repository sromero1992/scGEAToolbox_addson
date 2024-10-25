% Number of iterations for stability analysis
num_iterations = 1;
all_jaccard_ratios = zeros(num_iterations, 1);  % Store average Jaccard overlap ratios for each iteration
all_dice_ratios = zeros(num_iterations, 1);     % Store average Dice ratios for each iteration
all_overlap_ratios = zeros(num_iterations, 1);  % Store average Overlap Coefficient ratios for each iteration

all_jaccard_ratios2 = zeros(num_iterations, 1);  % Store average Jaccard overlap ratios for each iteration
all_dice_ratios2 = zeros(num_iterations, 1);     % Store average Dice ratios for each iteration
all_overlap_ratios2 = zeros(num_iterations, 1);  % Store average Overlap Coefficient ratios for each iteration

% Load the data and extract 'sce'
data = load("GSM4878207_GSM4878210.mat");
sce_tmp = data.sce;
clear data;

% Remove non-informative genes based on species
species = 'mouse';
sce_tmp = rm_non_inf_genes(sce_tmp, species);
batches = unique(sce_tmp.c_batch_id);

% Select cells belonging to each batch
idx1 = sce_tmp.c_batch_id == batches(1);
sce_tmp1 = sce_tmp.selectcells(idx1);
sce_tmp1 = sce_tmp1.qcfilter;

idx2 = sce_tmp.c_batch_id == batches(2);
sce_tmp2 = sce_tmp.selectcells(idx2);
sce_tmp2 = sce_tmp2.qcfilter;

% Call to differential variability analysis base
comparison_result = sc_splinedv(sce_tmp1.X, sce_tmp2.X, sce_tmp1.g, sce_tmp2.g, [], false);
idx_conf = comparison_result.pval <= 0.05;
comparison_result = comparison_result(idx_conf, :);
g0_dv = comparison_result.gene;

% Call to differential expression analysis base
[g_de, irows, jrows] = intersect(sce_tmp1.g, sce_tmp2.g, 'stable');
diff_exp_result = sc_deg(sce_tmp1.X(irows, :), sce_tmp2.X(jrows, :), g_de);
idx_conf = diff_exp_result.p_val <= 0.05;
diff_exp_result = diff_exp_result(idx_conf, :);
g0_de = diff_exp_result.gene;

clear diff_exp_result comparison_result idx_conf irows jrows g_de sce_tmp1 sce_tmp2;

progressbar('Initizalizing...');
progressbar(0); 
for iter = 1:num_iterations

    % Randomly group cells into new batches
    sce = randperm_batches2groups(sce_tmp, iter);

    % Get updated batch information
    batches = unique(sce.c_batch_id);
    nbatch = length(batches);

    % Calculate number of unique batch comparisons
    nupper = nbatch * (nbatch - 1) / 2;

    % Initialize a table to store results and metadata
    Tdv = table('Size', [nupper, 3], ...
                'VariableTypes', {'cell', 'cell', 'cell'}, ...
                'VariableNames', {'Result', 'BatchComparison', 'Metadata'});

    % Initialize a table to store results and metadata DE
    Tde = table('Size', [nupper, 3], ...
                'VariableTypes', {'cell', 'cell', 'cell'}, ...
                'VariableNames', {'Result', 'BatchComparison', 'Metadata'});
    it = 1;
    for ibatch = 1:nbatch
        ilfdbool = contains(batches(ibatch), 'LFD');
        ihfdbool = contains(batches(ibatch), 'HFD');

        for jbatch = ibatch + 1:nbatch

            jlfdbool = contains(batches(jbatch), 'LFD');
            jhfdbool = contains(batches(jbatch), 'HFD');

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
            tmp_name = strcat(batches(ibatch), "_", batches(jbatch));
            tmp_name = []; % for not writing file

            % Call to differential expression analysis (replace with actual function)
            comparison_result = sc_splinedv(sce1.X, sce2.X, sce1.g, sce2.g, tmp_name, false);
            %idx_conf = comparison_result.pval <= 0.05;
            %comparison_result = comparison_result(idx_conf, :);

            [g_de, irows, jrows] = intersect(sce1.g, sce2.g, 'stable');
            diff_exp_result = sc_deg(sce1.X(irows, :), sce2.X(jrows, :), g_de);
            %idx_conf = diff_exp_result.p_val <= 0.05;
            %diff_exp_result = diff_exp_result(idx_conf, :);

            % Storing result and meta-info in the table
            Tdv{it, 'Result'} = {comparison_result};
            Tde{it, 'Result'} = {diff_exp_result};
            tmp_name = strcat(batches(ibatch), "_", batches(jbatch));
            Tdv{it, 'BatchComparison'} = {tmp_name};
            Tde{it, 'BatchComparison'} = {tmp_name};
            Tdv{it, 'Metadata'} = {sprintf('Comparison between %s and %s', batches(ibatch), batches(jbatch))};
            Tde{it, 'Metadata'} = {sprintf('Comparison between %s and %s', batches(ibatch), batches(jbatch))};

            it = it + 1;
        end
    end

    % Remove empty rows or rows with missing values
    Tdv = Tdv(~cellfun(@isempty, Tdv.BatchComparison), :);
    Tde = Tde(~cellfun(@isempty, Tde.BatchComparison), :);

    % Compute overlap metrics between all pairs of batch comparisons
    nres = size(Tdv, 1);
    %nupper = nres * (nres - 1) / 2;
    nupper = nres;
    
    % Initialize a table to store overlap results
    Tovl_dv = table('Size', [nupper, 7], ...
                'VariableTypes', {'cell', 'cell', 'double', 'double', 'double', 'double', 'cell'}, ...
                'VariableNames', {'Overlap', 'Union', 'Jaccard', 'Dice', 'OverlapCoefficient', 'NumOverlap', 'BatchComparison'});

    % Initialize a table to store overlap results
    Tovl_de = table('Size', [nupper, 7], ...
                'VariableTypes', {'cell', 'cell', 'double', 'double', 'double', 'double', 'cell'}, ...
                'VariableNames', {'Overlap', 'Union', 'Jaccard', 'Dice', 'OverlapCoefficient', 'NumOverlap', 'BatchComparison'});

    it = 1;
    avg_jaccard_ratio = 0;
    avg_dice_ratio = 0;
    avg_overlap_ratio = 0;
    avg_jaccard_ratio2 = 0;
    avg_dice_ratio2 = 0;
    avg_overlap_ratio2 = 0;
    for itab = 1:nres
        g_dv = Tdv.Result{itab}.gene;
        g_de = Tde.Result{itab}.gene;

        % Calculate the intersection and union dv
        overlap = intersect(g_dv, g0_dv, 'stable');
        unionSet = union(g_dv, g0_dv);

        % Calculate the intersection and union de
        overlap2 = intersect(g_de, g0_de, 'stable');
        unionSet2 = union(g_de, g0_de);

        % Calculate Jaccard index, Dice coefficient, and Overlap coefficient dv
        novl = length(overlap);
        nunion = length(unionSet);
        n1 = length(g_dv);
        n2 = length(g0_dv);

        jaccard_ratio = 100 * novl / nunion;
        dice_ratio = 100 * (2 * novl) / (n1 + n2);
        overlap_ratio = 100 * novl / min(n1, n2);

        % Calculate Jaccard index, Dice coefficient, and Overlap coefficient de
        novl2 = length(overlap2);
        nunion2 = length(unionSet2);
        n3 = length(g_de);
        n4 = length(g0_de);

        jaccard_ratio2 = 100 * novl2 / nunion2;
        dice_ratio2 = 100 * (2 * novl2) / (n3 + n4);
        overlap_ratio2 = 100 * novl2 / min(n3, n4);

        % Storing results and meta-info in the table dv
        Tovl_dv{it, 'Overlap'} = {overlap};
        Tovl_dv{it, 'Union'} = {unionSet};
        Tovl_dv{it, 'Jaccard'} = jaccard_ratio;
        Tovl_dv{it, 'Dice'} = dice_ratio;
        Tovl_dv{it, 'OverlapCoefficient'} = overlap_ratio;
        Tovl_dv{it, 'NumOverlap'} = novl;

        % Storing results and meta-info in the table dv
        Tovl_de{it, 'Overlap'} = {overlap2};
        Tovl_de{it, 'Union'} = {unionSet2};
        Tovl_de{it, 'Jaccard'} = jaccard_ratio2;
        Tovl_de{it, 'Dice'} = dice_ratio2;
        Tovl_de{it, 'OverlapCoefficient'} = overlap_ratio2;
        Tovl_de{it, 'NumOverlap'} = novl2;

        tab1_name = strcat('Table', Tdv{it,"BatchComparison"});
        Tovl_dv{it, 'BatchComparison'} = {tab1_name};
        Tovl_de{it, 'BatchComparison'} = {tab1_name};

        it = it + 1;

        % Update averages
        avg_jaccard_ratio = avg_jaccard_ratio + jaccard_ratio;
        avg_dice_ratio = avg_dice_ratio + dice_ratio;
        avg_overlap_ratio = avg_overlap_ratio + overlap_ratio;

        avg_jaccard_ratio2 = avg_jaccard_ratio2 + jaccard_ratio2;
        avg_dice_ratio2 = avg_dice_ratio2 + dice_ratio2;
        avg_overlap_ratio2 = avg_overlap_ratio2 + overlap_ratio2;
    end

    % Calculate and store the average ratios for this iteration
    all_jaccard_ratios(iter) = avg_jaccard_ratio / nupper;
    all_dice_ratios(iter) = avg_dice_ratio / nupper;
    all_overlap_ratios(iter) = avg_overlap_ratio / nupper;
    
    all_jaccard_ratios2(iter) = avg_jaccard_ratio2 / nupper;
    all_dice_ratios2(iter) = avg_dice_ratio2 / nupper;
    all_overlap_ratios2(iter) = avg_overlap_ratio2 / nupper;
    fprintf("DV-Iteration %d Jaccard ratio: %f, Dice ratio: %f, Overlap Coefficient: %f \n", iter, all_jaccard_ratios(iter), all_dice_ratios(iter), all_overlap_ratios(iter));
    fprintf("DE-Iteration %d Jaccard ratio: %f, Dice ratio: %f, Overlap Coefficient: %f \n", iter, all_jaccard_ratios2(iter), all_dice_ratios2(iter), all_overlap_ratios2(iter));

    progressbar(iter)
end

% Statistical summary of overlap ratios across all iterations
mean_jaccard = mean(all_jaccard_ratios);
std_jaccard = std(all_jaccard_ratios);
mean_dice = mean(all_dice_ratios);
std_dice = std(all_dice_ratios);
mean_overlap = mean(all_overlap_ratios);
std_overlap = std(all_overlap_ratios);

fprintf('------ Spline-DV results -------------\n');
fprintf('Mean Jaccard ratio across iterations: %.2f\n', mean_jaccard);
fprintf('Standard deviation of Jaccard ratio: %.2f\n', std_jaccard);
fprintf('Mean Dice ratio across iterations: %.2f\n', mean_dice);
fprintf('Standard deviation of Dice ratio: %.2f\n', std_dice);
fprintf('Mean Overlap Coefficient across iterations: %.2f\n', mean_overlap);
fprintf('Standard deviation of Overlap Coefficient: %.2f\n', std_overlap);

mean_jaccard2 = mean(all_jaccard_ratios2);
std_jaccard2 = std(all_jaccard_ratios2);
mean_dice2 = mean(all_dice_ratios2);
std_dice2 = std(all_dice_ratios2);
mean_overlap2 = mean(all_overlap_ratios2);
std_overlap2 = std(all_overlap_ratios2);

fprintf('------ Differential Expression results -------------\n');
fprintf('Mean Jaccard ratio across iterations: %.2f\n', mean_jaccard2);
fprintf('Standard deviation of Jaccard ratio: %.2f\n', std_jaccard2);
fprintf('Mean Dice ratio across iterations: %.2f\n', mean_dice2);
fprintf('Standard deviation of Dice ratio: %.2f\n', std_dice2);
fprintf('Mean Overlap Coefficient across iterations: %.2f\n', mean_overlap2);
fprintf('Standard deviation of Overlap Coefficient: %.2f\n', std_overlap2);

% Visualization of overlap ratios
figure;
histogram(all_jaccard_ratios, 'Normalization', 'probability');
title('Distribution of Jaccard Ratios Across Iterations');
xlabel('Jaccard Ratio (%)');
ylabel('Frequency');

figure;
histogram(all_dice_ratios, 'Normalization', 'probability');
title('Distribution of Dice Ratios Across Iterations');
xlabel('Dice Ratio (%)');
ylabel('Frequency');

figure;
histogram(all_overlap_ratios, 'Normalization', 'probability');
title('Distribution of Overlap Coefficients Across Iterations');
xlabel('Overlap Coefficient (%)');
ylabel('Frequency');