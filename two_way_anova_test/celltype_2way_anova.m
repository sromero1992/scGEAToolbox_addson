%fname = "AllSamplesMerged_HMOs6Removed_OldGenesRemoved_1500-20-15-500Filtered_DoubletsRemoved_AmbientRNARemoved_1500-20-15-500Filtered_AllGenesUMAP_Col_LeidenRes1_29Clusters_CellTypes_UseThisOne.mat";
fname = "simplified_batches_binf_hmo.mat";
% Load data
data = load(fname);
sce = data.sce;
clear data;

gene_set = 'all';
cell_type = "Intestinal stem cells";

% one-way and two-way anova tests
%results = run_anova_tests(sce, gene_set, cell_type);
conditions = {'B_infantis','HMOs'}; % Avoid - or + in the batch names
interaction = true; % B_infantis*HMOs
tic;
results = run_anova_tests2(sce, gene_set, cell_type, conditions, interaction);
time_end = toc;
fprintf('Results of two-way ANOVA tests for all genes %f \n', time_end);

% Summarize pvalues from anova tests
pvalues_table = export_anova_results(results, 'anova_results.csv');

% Extract gene two way anova and pvalue adjusted
gene_name = 'Mtdh';
gene_data_row = results(strcmp(results.Gene, gene_name), :);
pvalue_row =  pvalues_table(gene_name,:);

save_filename = strcat(gene_name,'_expressionViolin.png'); % Or any other filename
plotGeneExpression(gene_data_row, pvalue_row, conditions, gene_name, save_filename)


% save_filename = strcat(gene_name,'_expressionBoxplot.png'); % Or any other filename
% plotExpressionBoxplot(gene_data_row, pvalue_row, conditions, gene_name, save_filename)
