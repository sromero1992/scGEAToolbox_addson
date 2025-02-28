function anova_results = run_anova_tests2(sce, gene_set, cell_type, conditions, interact)

    if nargin < 5; interact = false; end

    % Subset cell type
    idx = cell_type == sce.c_cell_type_tx;
    sce = sce.selectcells(idx);
    %sce = sce.qcfilter(500, 0.20, 10);

    if gene_set == 'all'
        gene_set = sce.g;
    end
    
    fprintf("Number of cells %d in cell type %s \n", sce.NumCells, cell_type);
    % Gene set to analyze
    nset = length(gene_set);
    gene_idx = zeros(nset, 1);

    % Find indices of genes in gene_set
    for ig = 1:nset
        gene_idx0 = find(strcmp(sce.g, gene_set(ig)));
        if isempty(gene_idx0)
            error('Gene "%s" not found in sce.g.', gene_set(ig));
        end
        gene_idx(ig) = gene_idx0;
    end

    % Convert to categorical variables
    gbatch_unique =  unique(sce.c_batch_id);
    batches_treat = categorical(sce.c_batch_id, gbatch_unique, 'Ordinal', true);
    
    % Normalize the data
    %Xnorm = sc_norm(sce.X);
    Xnorm = pkg.norm_libsize(sce.X, 1e4);
    Xnorm = log1p(Xnorm);

    % Initialize a cell array to store results
    all_results = cell(size(gene_idx, 1), 1);

    % Parallel loop
    if nset > 500
        %numCores = ceil(feature('numcores')/4);
        %parpool(numCores)
        parfor geneIdx = 1:nset  % Use parfor for parallel processing
            expression = Xnorm(gene_idx(geneIdx), :);
            expression = expression(:);

            cond1 = contains(string(batches_treat), conditions{1});
            cond2 = contains(string(batches_treat), conditions{2});
            tbl3 = table(expression, cond1, cond2, 'VariableNames', {'Expression', conditions{1}, conditions{2}});

            if interact
                % Model is Expression ~ condition1 + condition2 + condition1*condition2
                model = sprintf('Expression ~ %s + %s + %s:%s', conditions{1}, conditions{2}, conditions{1}, conditions{2} );
            else
                % Model is Expression ~ condition1 + condition2
                model = sprintf('Expression ~ %s + %s', conditions{1}, conditions{2});
            end
            aov2 = anova(tbl3, model);

            anova_struct = struct('Gene', sce.g{gene_idx(geneIdx)}, ...
                                    'ANOVA_2way', aov2, ...
                                    'Expression', expression, ...
                                    'Treatment', batches_treat, ...
                                    conditions{1}, cond1, ...
                                    conditions{2}, cond2);

            all_results{geneIdx} = anova_struct;
        end
    else  % Serial loop for smaller gene sets
        for geneIdx = 1:nset
            expression = Xnorm(gene_idx(geneIdx), :);
            expression = expression(:);

            cond1 = contains(string(batches_treat), conditions{1});
            cond2 = contains(string(batches_treat), conditions{2});
            tbl3 = table(expression, cond1, cond2, 'VariableNames', {'Expression', conditions{1}, conditions{2}});

            if interact
                % Model is Expression ~ condition1 + condition2 + condition1*condition2
                model = sprintf('Expression ~ %s + %s + %s:%s', conditions{1}, conditions{2}, conditions{1}, conditions{2} );
            else
                % Model is Expression ~ condition1 + condition2
                model = sprintf('Expression ~ %s + %s', conditions{1}, conditions{2});
            end
            aov2 = anova(tbl3, model);

            anova_struct = struct('Gene', sce.g{gene_idx(geneIdx)}, ...
                                    'ANOVA_2way', aov2, ...
                                    'ANOVA_1way', aov1, ...
                                    'Pairwise_1way', pairwise_comparisons, ...
                                    'Expression', expression, ...
                                    'Treatment', batches_treat, ...
                                    conditions{1}, cond1, ...
                                    conditions{2}, cond2);

            all_results{geneIdx} = anova_struct;
        end
    end

    anova_results = struct2table(cell2mat(all_results));
end
