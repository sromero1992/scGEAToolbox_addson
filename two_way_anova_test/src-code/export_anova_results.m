function pvalues_table = export_anova_results(results, output_file)
    % Export ANOVA results from `results` into a CSV or TSV file.
    % Arguments:
    % - results: A table containing the output from run_anova_tests.
    % - output_file: String specifying the output file name (e.g., 'output.csv').

    % Preallocate for efficiency
    nGenes = size(results, 1);
    all_pvalues = cell(nGenes, 1);  % Store p-values as cell array initially
    gene_names = cell(nGenes, 1);
    all_rownames = {};
    
    colnames = string(results.Properties.VariableNames);

    % Loop over the `results` structure
    for i = 1:nGenes
        gene_pvalues = [];
        gene_rownames = {};
        gene_names{i} = results{i, 1}{1}; % Gene name (row name)
        
        for j = 2:size(results, 2) % Start from 2, gene name already extracted
            current_result = results{i, j}{1};
            if contains(class(current_result), 'anova')
                tbl_ext = current_result.stats;
                pvalues = tbl_ext.pValue;
                rownames = tbl_ext.Properties.RowNames;
                valid_rows = ~ismember(string(rownames), ["Error", "Total"]);
                pvalues = pvalues(valid_rows);
                rownames = string(rownames(valid_rows));
                rownames = strcat(rownames, "_", colnames(j));
            else
                continue;
            end
            gene_pvalues = [gene_pvalues; pvalues(:)];  % append p-values
            gene_rownames = [gene_rownames; rownames(:)]; % append row names
        end
        all_pvalues{i} = gene_pvalues;
        if i == 1 % Initialize all_rownames with the first gene's rownames. All genes should have the same row names.
          all_rownames = gene_rownames;
        elseif ~isequal(all_rownames, gene_rownames)
            error("Row names are inconsistent across genes. Check your input data.");
        end
    end

    % Convert cell array to numeric matrix
    nrows = length(all_rownames);
    all_pvalues_reshaped = zeros(nGenes, nrows);
    for i = 1:nGenes
      all_pvalues_reshaped(i,:) = all_pvalues{i}';
    end


    % Ensure that all_rownames has the correct size for nrows
    rownames_str = string(all_rownames(1:nrows));

    % Create the final table with gene names and their corresponding p-values
    pvalues_table = array2table(all_pvalues_reshaped, 'VariableNames', rownames_str);
    pvalues_table.Properties.RowNames = gene_names;


    subres = results{1,'ANOVA_2way'}{1}.stats;
    rownames = subres.Properties.RowNames;
    valid_rows = ~ismember(string(rownames), ["Error", "Total"]);
    ntests = sum(valid_rows);

    padj1 = bh_adjust_pvalues(pvalues_table{:,1});
    padj2 = bh_adjust_pvalues(pvalues_table{:,2});
    if ntests > 2
        padj3 = bh_adjust_pvalues(pvalues_table{:,3});
    end

    % Add adjusted p-values as new columns
    pvalues_table.([pvalues_table.Properties.VariableNames{1}, '_padj']) = padj1;
    pvalues_table.([pvalues_table.Properties.VariableNames{2}, '_padj']) = padj2;
    
    if ntests > 2
        pvalues_table.([pvalues_table.Properties.VariableNames{3}, '_padj']) = padj3;
    end

    pval_2way1_bool = pvalues_table{:, 4} <= 0.05;
    pval_2way2_bool = pvalues_table{:, 5} <=0.05;
    % If testing interaction or not
    if ntests > 2 
        pval_2way3_bool = pvalues_table{:, 6} <=0.05;
        idx = pval_2way1_bool | pval_2way2_bool | pval_2way3_bool;
    else
        idx = pval_2way1_bool | pval_2way2_bool;
    end

    pvalues_table = pvalues_table(idx,:);
    pvalues_table = sortrows(pvalues_table, [6, 4, 5], ["ascend", "ascend", "ascend" ]);

    % Write to the specified output file
    [~, ~, ext] = fileparts(output_file);
    if strcmpi(ext, '.csv')
        writetable(pvalues_table, output_file, 'WriteRowNames', true);
    elseif strcmpi(ext, '.tsv')
        writetable(pvalues_table, output_file, 'WriteRowNames', true, 'FileType', 'text', 'Delimiter', '\t');
    else
        error('Unsupported file format. Use .csv or .tsv.');
    end
end
