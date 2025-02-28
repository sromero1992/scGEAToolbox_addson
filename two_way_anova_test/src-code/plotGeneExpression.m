function plotGeneExpression(gene_data, pvalue_row, conditions, gene_name, save_filename)
    % Plots gene expression data for a specified gene from a results table.
    % Optionally saves the plot to a file if save_filename is provided.
    if isempty(gene_data)
        fprintf('Gene "%s" not found in the results table.\n', gene_name);
        return; % Exit the function if the gene is not found
    end
    gene_data0 = struct();
    gene_data0.Expression = gene_data.Expression{1};
    gene_data0.Treatment = gene_data.Treatment{1};
    cond1 = conditions{1};
    gene_data0.(cond1) = gene_data{1, cond1}{1};
    cond2 = conditions{2};
    gene_data0.(cond2) = gene_data{1, cond2}{1};
    cond12 = strcat(cond1,'+',cond2); 
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
    % Extract p-values (Corrected)
    p_cond1 = pvalue_row{1,1};
    p_cond1_adj = pvalue_row{1,4};
    p_cond2 = pvalue_row{1,2};
    p_cond2_adj = pvalue_row{1,5};
    % Calculate min and max expression
    min_expression = min(gene_data0.Expression);
    max_expression = max(gene_data0.Expression);
    % 2x2 Grid of Violin Plots
    figure;
    
    % Adjust figure size and margins for better title and subplot spacing
    set(gcf, 'Units', 'inches');
    fig_pos = get(gcf, 'Position');
    fig_pos(3) = fig_pos(3) * 1.5; % Increase width
    fig_pos(4) = fig_pos(4) * 1.7; % Significantly increase height for title
    set(gcf, 'Position', fig_pos);
    
    colors = {[0.5 0.5 0], [1 0 0], [0 0 1], [0.4940 0.1840 0.5560]};
    alpha = 0.2;
    for iSubplot = 1:4
        subplot(2, 2, iSubplot);
        
        % Adjust subplot margins
        ax = gca;
        ax.Position(1) = ax.Position(1) + 0.05; % Move subplot to the right
        ax.Position(2) = ax.Position(2) + 0.05; % Move subplot up
        ax.Position(3) = ax.Position(3) * 0.9; % Reduce width
        ax.Position(4) = ax.Position(4) * 0.9; % Reduce height
        
        switch iSubplot
            case 1, groupName = 'Control';
            case 2, groupName = cond1;
            case 3, groupName = cond2;
            case 4, groupName = cond12;
        end
        violin_object = violinplot(gene_data0.Expression(interaction_group == groupName), ...
                                   'DensityScale','count');
                                    % 'DensityWidth', 1.2);
        % Remove the title
        % title(strrep(groupName, '_',' '), 'FontSize', 12); 
        xlabel(strrep(groupName, '_', ' '), 'FontSize', 14); % Add title as xlabel
        ylim([min_expression, max_expression]);
        hold on;
        median_val = median(gene_data0.Expression(interaction_group == groupName));
        line([0.5 1.5], [median_val median_val], 'Color', 'k', 'LineWidth', 2);
        hold off;
        % Customize violin colors
        if ~isempty(violin_object)
            set(violin_object, 'FaceColor', colors{iSubplot}, 'FaceAlpha', alpha);
        end
    end
    % Title and Annotations
    han = axes('Position', [0 0 1 1], 'Visible', 'off');
    set(gcf, 'CurrentAxes', han);
    if length(pvalue_row.Properties.VariableTypes) == 6
        tmp_str = sprintf('- Interaction pval-adj = %.3f', pvalue_row{1, 6});
    else
        tmp_str = '';
    end
        
    text(0.5, 0.98, ['Expression of ', gene_name, ' ',tmp_str ], 'HorizontalAlignment', 'center', 'FontSize', 16);
    
    annotation_x = [0.48, 0.18];
    annotation_y = [0.01, 0.35];
    annotation_text = {sprintf('padj(%s) = %.3f', strrep(cond1, '_',' '), p_cond1_adj), ...
                       sprintf('padj(%s) = %.3f', strrep(cond2, '_',' '), p_cond2_adj)};
    for iAnnot = 1:2
        rotation_angle = 0;
        if iAnnot == 2
            rotation_angle = 90;
        end
        annotation('textbox', [annotation_x(iAnnot), annotation_y(iAnnot), 0.1, 0.1], ...
                   'String', annotation_text{iAnnot}, 'FontSize', 12, ...
                   'HorizontalAlignment', 'center', 'Rotation', rotation_angle);
    end
    % Save the plot if save_filename is provided
    if nargin > 2 && ~isempty(save_filename)  % Check if save_filename argument exists and is not empty
        print(gcf, '-dpng', save_filename); % Save as PNG (you can change the format)
        close(gcf); % Close the figure to prevent it from being displayed
        fprintf('Plot saved to "%s"\n', save_filename);
    end
    
end % End of the function