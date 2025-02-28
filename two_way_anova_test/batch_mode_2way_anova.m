%fname = "AllSamplesMerged_HMOs6Removed_OldGenesRemoved_1500-20-15-500Filtered_DoubletsRemoved_AmbientRNARemoved_1500-20-15-500Filtered_AllGenesUMAP_Col_LeidenRes1_29Clusters_CellTypes_UseThisOne.mat";
fname = "simplified_batches_binf_hmo.mat";
% Load data
data = load(fname);
sce = data.sce;
clear data;

gene_set = 'all';

% one-way and two-way anova tests
conditions = {'B_infantis','HMOs'}; % Avoid - or + in the batch names
interaction = true; % B_infantis*HMOs

cell_types = unique(sce.c_cell_type_tx)';
for cell_type = cell_types
    disp(cell_type)
    tic;
    results = run_anova_tests2(sce, gene_set, cell_type, conditions, interaction);
    time_end = toc;
    fprintf('Results of two-way ANOVA tests for all genes %f \n', time_end);
    fname = strcat( cell_type,'_2way_anova_results.csv');
    tic;
    pvalues_table = export_anova_results(results, fname);
    time_end = toc;
    fprintf('Pvalues file finished in %f \n', time_end);
end



% Extract gene two way anova and pvalue adjusted
gene_name = 'Mtdh';
gene_data_row = results(strcmp(results.Gene, gene_name), :);
pvalue_row =  pvalues_table(gene_name,:);

save_filename = strcat(gene_name,'_expressionViolin.png'); % Or any other filename
plotGeneExpression(gene_data_row, pvalue_row, conditions, gene_name, save_filename)


% save_filename = strcat(gene_name,'_expressionBoxplot.png'); % Or any other filename
% plotExpressionBoxplot(gene_data_row, pvalue_row, conditions, gene_name, save_filename)
