function [final_table, table_cross, table_within ]= analyze_stability_dv(ngenes, num_iter, ...
                                                   species, batch1, batch2, mat_file_path)
    % SplineDV analysis tool for random permutation test and subsample
    % 50%/50% of the whole batch1 and batch2, this function depends of sce
    % object defined in scGEAToolbox
    % Author: Selim Romero
    % INPUTS:
    % ngenes --------> Number of top DV genes to analyze
    % num_iter ------> Number of iterations to subsample and do the DV analysis
    % species -------> Species ('human' or 'mouse')
    % batch1 --------> SCE object should contain batch1 string to subsample
    % batch2 --------> SCE object should contain batch2 string to subsample
    % mat_file_path -> Path to the sce object to load
    % OUTPUT:
    % final_table ---> Final table containing all the SplineDV distance
    %                  difference to asses further ploting/analysis
    % table_cross----> Similar to final but contains only cross sample comparisons
    % table_within---> Similar to final but contains only within sample comparisons
    % Addditionaly writes those tables as:
    %                 stability_analysisDV_cross.csv
    %                 stability_analysisDV_within.csv
    %                 stability_analysisDV_final.csv
    % USAGE:
    % mat_file = "C:\Users\ssromerogon\Documents\Mouse_Adipose_Tissue_1\GSM4878207_GSM4878210.mat"
    % final_table = analyze_stability_dv(50, 10, 'mouse', 'HFD', 'LFD', mat_file);

    % Load the data and extract 'sce'
    data = load(mat_file_path);
    sce_tmp = data.sce;
    clear data;

    % Remove non-informative genes based on species
    sce_tmp = rm_non_inf_genes(sce_tmp, species);
    batches = unique(sce_tmp.c_batch_id);

    % Ensure there are at least two batches
    if length(batches) < 2
        error('Not enough batches for analysis.');
    end

    sub_str1 = batch1;
    sub_str2 = batch2;
    idx_batch = batches == batch1;
    jdx_batch = batches == batch2;

    % Select cells belonging to each batch
    idx1 = sce_tmp.c_batch_id == batches(idx_batch);
    sce_tmp1 = sce_tmp.selectcells(idx1);
    sce_tmp1 = sce_tmp1.qcfilter;

    idx2 = sce_tmp.c_batch_id == batches(jdx_batch);
    sce_tmp2 = sce_tmp.selectcells(idx2);
    sce_tmp2 = sce_tmp2.qcfilter;

    clear sce_tmp;
    sce_tmp = sc_mergesces({sce_tmp1, sce_tmp2},'union');
    sce_tmp.c = sce_tmp.c_batch_id;

    % Call to differential variability analysis base
    Tdv0 = sc_splinedv(sce_tmp1.X, sce_tmp2.X, sce_tmp1.g, sce_tmp2.g, [], false);

    % Define the gene list 'gl' as the first 'ngenes' entries from Tdv0
    gl = Tdv0.gene(1:ngenes);
    Tdv0_sub = Tdv0(ismember(Tdv0.gene, gl), :);

    % Initialize the main table to hold all gene information and concatenated columns
    mainTable = table(gl, 'VariableNames', {'gene'});

    % Add the new column to 'mainTable'
    mainTable.('diff_dist_cross0') = Tdv0_sub.dist_diff;

    % Loop for multiple iterations and comparisons
    iter0cross = 1;
    iter0within = 1;
    progressbar('DV stability analysis...');
    progressbar(0);
    for iter = 1:num_iter
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
                    colName = sprintf('diff_dist_within%d', iter0within);
                    iter0within = iter0within + 1;
                else
                    colName = sprintf('diff_dist_cross%d', iter0cross);
                    iter0cross = iter0cross + 1;
                end

                % Add the new column to 'mainTable'
                mainTable.(colName) = Tdv1_sub_new.dist_diff;
            end
        end
        progressbar(iter/num_iter*100);
    end

    % Step 1: Separate columns based on naming
    cross_cols = contains(mainTable.Properties.VariableNames, 'diff_dist_cross');
    within_cols = contains(mainTable.Properties.VariableNames, 'diff_dist_within');

    % Keep 'gene' and cross columns
    cross_cols(1) = true;  % Keep the 'gene' column
    within_cols(1) = true;    % Keep the 'gene' column

    % Create separate tables
    cross_table = mainTable(:, cross_cols);  % Keep 'gene' and cross columns
    within_table = mainTable(:, within_cols);    % Keep 'gene' and within columns

    % Step 2: Transpose each table
    table_cross = array2table(cross_table{:, 2:end}', 'VariableNames', cross_table.gene);
    table_within = array2table(within_table{:, 2:end}', 'VariableNames', within_table.gene);

    % Step 3: Add labels to identify them as strings
    table_cross.Type = repmat("cross", height(table_cross), 1);  % Use double quotes for strings
    table_within.Type = repmat("within", height(table_within), 1);    % Use double quotes for strings

    % Combine both tables
    final_table = [table_cross; table_within];

    % Move 'Type' column to the start of the table
    final_table = final_table(:, [end, 1:end-1]);  % Move 'Type' to the first column
    table_cross = table_cross(:, [end, 1:end-1]);  % Move 'Type' to the first column
    table_within = table_within(:, [end, 1:end-1]);  % Move 'Type' to the first column

    % Save final table to CSV
    writetable(table_cross, 'stability_analysisDV_cross.csv');
    writetable(table_within, 'stability_analysisDV_within.csv');
    writetable(final_table, 'stability_analysisDV_final.csv');

    % Return the final table
    return
end
