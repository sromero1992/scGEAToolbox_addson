function plotExpressionBoxplot(gene_data_row, pvalue_row, conditions, gene_name, save_filename)
    % Plots gene expression (Y) partitioned by conditions (X).
    % Optionally saves the plot to a file if save_filename is provided.

    if isempty(gene_data_row)
        fprintf('Gene "%s" not found in the results table.\n', gene_name);
        return;
    end

    gene_data0 = struct();
    gene_data0.Expression = gene_data_row.Expression{1};
    gene_data0.Treatment = gene_data_row.Treatment{1}; % Assuming treatment is present.

    cond1 = conditions{1};
    gene_data0.(cond1) = gene_data_row.(cond1){1};
    cond2 = conditions{2};
    gene_data0.(cond2) = gene_data_row.(cond2){1};
    cond12 = strcat(cond1, '+', cond2);

    interaction_group = cell(size(gene_data0.Expression, 1), 1);
    for i = 1:size(gene_data0.Expression, 1)
        if gene_data0.(cond1)(i) && gene_data0.(cond2)(i)
            interaction_group{i} = cond12;
        elseif gene_data0.(cond1)(i)
            interaction_group{i} = cond1;
        elseif gene_data0.(cond2)(i)
            interaction_group{i} = cond2;
        else
            interaction_group{i} = 'Control';
        end
    end

    new_category_order = {'Control', cond1, cond2, cond12}; % Desired order
    interaction_group = categorical(interaction_group);
    interaction_group = reordercats(interaction_group, new_category_order);

    % Extract p-values
    p_cond1 = pvalue_row{1, 1};
    p_cond1_adj = pvalue_row{1, 4};
    p_cond2 = pvalue_row{1, 2};
    p_cond2_adj = pvalue_row{1, 5};

    % Prepare data for plotting
    expression_data = gene_data0.Expression;

    % Create boxchart
    figure;
    
    % Define colors and linewidth
    colors = {[0.4940 0.1840 0.5560], [0 0 1], [1 0 0], [0 1 0]};
    lineWidth = 2;

    % Create boxplot with black outlines
    boxplot(expression_data, interaction_group, 'BoxStyle', 'outline', ...
            'Colors', 'k', 'MedianStyle', 'line', 'Widths', 0.5);

    h = findobj(gca, 'Tag', 'Box');
    for j = 1:length(h)
        set(h(j), 'LineWidth', lineWidth);
        patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{j}, 'FaceAlpha', 0.2); % Use original colors for fill
    end

    ylabel('Expression');

    % Annotations
    if length(pvalue_row.Properties.VariableTypes) == 6
        tmp_str = sprintf('Interaction pval-adj = %.3f', pvalue_row{1, 6});
    else
        tmp_str = '';
    end

    title({['Expression of ', gene_name], tmp_str}, 'HorizontalAlignment', ...
             'center', 'Units', 'normalized', 'Position', [0.5, 0.995, 0]);

    % annotation_x = [0.15, 0.35];
    % annotation_y = [0.8, 0.8];
    % annotation_text = {sprintf('padj(%s) = %.3f', cond1, p_cond1_adj), ...
    %                    sprintf('padj(%s) = %.3f', cond2, p_cond2_adj)};
    % 
    % for iAnnot = 1:2
    %     annotation('textbox', [annotation_x(iAnnot), annotation_y(iAnnot), 0.1, 0.1], ...
    %                'String', annotation_text{iAnnot}, 'FontSize', 10, ...
    %                'HorizontalAlignment', 'center');
    % end

    % Save the plot
    if nargin > 2 && ~isempty(save_filename)
        print(gcf, '-dpng', save_filename);
        close(gcf);
        fprintf('Plot saved to "%s"\n', save_filename);
    end
end

