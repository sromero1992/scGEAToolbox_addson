% Initialize the gene list and number of iterations
ngenes = 50; % Number of genes to consider
num_iterations = 100; % Number of iterations for resampling

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

% Define the gene list 'gl' as the first 50 entries from Tdv0
gl = Tdv0.gene(1:ngenes);
Tdv0_sub = Tdv0(ismember(Tdv0.gene, gl), :);

% Initialize the main table to hold all gene information and concatenated columns
mainTable = table(gl, 'VariableNames', {'gene'});

% Add the new column to 'mainTable'
mainTable.('diff_dist_good0') = Tdv0_sub.dist_diff;

% Loop for multiple iterations and comparisons
iter0Good = 1;
iter0Bad = 1;
progressbar(0);
for iter = 1:num_iterations
    % Resampling: randomly group cells into new batches
    sce = randperm_batches2groups(sce_tmp, iter);
    batches = unique(sce.c_batch_id);
    nbatch = length(batches);
    for ibatch = 1:nbatch
        for jbatch = ibatch + 1:nbatch
            % Select cells for each batch and perform quality control
            ibname = batches(ibatch);
            idx1 = sce.c_batch_id == ibname;
            sce1 = sce.selectcells(idx1);
            sce1 = sce1.qcfilter;
            
            jbname = batches(jbatch);
            idx2 = sce.c_batch_id == jbname;
            sce2 = sce.selectcells(idx2);
            sce2 = sce2.qcfilter;

            % Perform differential variability analysis
            Tdv1 = sc_splinedv(sce1.X, sce2.X, sce1.g, sce2.g, [], false);

            % Match the resulting table with the gene list 'gl'
            [~, idx] = ismember(gl, Tdv1.gene);

            % Create a new table with 'gl' as the row reference
            Tdv1_sub_new = array2table(NaN(length(gl), 1), 'VariableNames', {'dist_diff'});
            Tdv1_sub_new.gene = gl;

            % Replace rows with the matched genes
            validIdx = idx > 0;
            Tdv1_sub_new.dist_diff(validIdx) = Tdv1.dist_diff(idx(validIdx));

            % Generate column names dynamically based on comparison type
            colName = "";
            % Removes any pattern like '-seed-XXXX'
            ibname = regexprep(ibname, '-seed-\d+', '');  
            jbname = regexprep(jbname, '-seed-\d+', '');  
            ib1bool = contains(ibname, sub_str1);
            jb1bool = contains(jbname, sub_str1);
            ib2bool = contains(ibname, sub_str2);
            jb2bool = contains(jbname, sub_str2);

            if (ib1bool && jb1bool) || (ib2bool && jb2bool)
                colName = sprintf('diff_dist_bad%d', iter0Bad);
                iter0Bad = iter0Bad + 1;
            else
                colName = sprintf('diff_dist_good%d', iter0Good);
                iter0Good = iter0Good + 1;
            end

            % Add the new column to 'mainTable'
            mainTable.(colName) = Tdv1_sub_new.dist_diff;
        end
    end
    progressbar(iter);
end

% Step 1: Separate columns based on naming
good_cols = contains(mainTable.Properties.VariableNames, 'diff_dist_good');
bad_cols = contains(mainTable.Properties.VariableNames, 'diff_dist_bad');

% Keep 'gene' and good columns
good_cols(1) = true;  % Keep the 'gene' column
bad_cols(1) = true;    % Keep the 'gene' column

% Create separate tables
good_table = mainTable(:, good_cols);  % Keep 'gene' and good columns
bad_table = mainTable(:, bad_cols);    % Keep 'gene' and bad columns

% Step 2: Transpose each table
good_table_transposed = array2table(good_table{:, 2:end}', 'VariableNames', good_table.gene);
bad_table_transposed = array2table(bad_table{:, 2:end}', 'VariableNames', bad_table.gene);

% Step 3: Add labels to identify them as strings
good_table_transposed.Type = repmat("good", height(good_table_transposed), 1);  % Use double quotes for strings
bad_table_transposed.Type = repmat("bad", height(bad_table_transposed), 1);    % Use double quotes for strings

% Combine both tables
final_table = [good_table_transposed; bad_table_transposed];

% Move 'Type' column to the start of the table
final_table = final_table(:, [end, 1:end-1]);  % Move 'Type' to the first column

writetable(final_table, 'stability_analysisDV.csv')

