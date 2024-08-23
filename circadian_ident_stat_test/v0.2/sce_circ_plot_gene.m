function sce_circ_plot_gene(sce, tmeta, cust_cells, period12, cust_gene, axHandle)
    % Plot identified circadian gene for the specified custom gene

    % Ensure required inputs
    if nargin < 4 || isempty(period12); period12 = false; end
    if nargin < 5 || isempty(cust_gene); error('Custom gene must be specified.'); end
    if nargin < 6 || isempty(axHandle); error('Axes handle must be specified.'); end

    % Define time variables
    tmeta.times = sortrows(tmeta.times);
    t0 = tmeta.times(1);
    tint = mean(diff(tmeta.times));
    disp("Time steps are : " + tint);
    tf = tmeta.times(end);
    t = t0:tint:tf;
    tval = t0:0.1:tf;

    % Compute circadian information for the custom gene and cell type
    [T1, T2] = sce_circ_phase_estimation_stattest(sce, tmeta, false, period12, cust_gene, cust_cells);

    % Find the specified gene in the results
    gene_idx = find(strcmp(T1.Genes, cust_gene));

    % Generate the plot for the specified gene if found
    if ~isempty(gene_idx)
        % Set the current axes to the provided axHandle
        axes(axHandle);

        % Generate sine-fitted values
        fval = T1.Amp(gene_idx) * cos(2*pi*(tval - T1.Acrophase(gene_idx)) / T1.Period(gene_idx)) + T1.Mesor(gene_idx);
        Rzts = table2array(T2(gene_idx, 2:end));

        % Plot sine-fitted values
        plot(tval, fval, 'LineWidth', 2);
        hold on;

        % Plot original expression data
        plot(t, Rzts, 'o-', 'LineWidth', 1.5);

        % Set plot limits and labels
        xlim([t0 tf]);
        title(sprintf('Gene - %s | p-value: %.3f', T1.Genes{gene_idx}, T1.pvalue(gene_idx)));
        xlabel('Time (hrs)');
        ylabel('Expression');
        legend({'Sine-fitted expression', 'Expression'}, 'Location', 'northwest');
        hold off;
    else
        error('Specified gene not found in the dataset.');
    end
end
