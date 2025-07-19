function sce_circ_plot_gene(sce, tmeta, cust_cells, period12, cust_gene, ...
                            axHandle, print_scdata, norm_str)
    % Plot identified circadian gene for the specified custom gene

    % Ensure required inputs
    if nargin < 4 || isempty(period12); period12 = false; end
    if nargin < 5 || isempty(cust_gene); error('Custom gene must be specified.'); end
    if nargin < 6 || isempty(axHandle); error('Axes handle must be specified.'); end
    if nargin < 7 || isempty(print_scdata); print_scdata = false; end
    if nargin < 8 || isempty(norm_str); norm_str = 'lib_size'; end

    % Define time variables from metadata
    tmeta.times = sortrows(tmeta.times); % Sort times if they aren't in order
    t0 = tmeta.times(1); % Initial time
    tint = mean(diff(tmeta.times)); % Calculate average time interval
    disp("Time steps are : " + tint);
    tf = tmeta.times(end); % Final time
    t = tmeta.times; % Use actual times from tmeta
    tval = t0:0.1:tf; % Finer time points for sine-fitted values

    % Compute circadian information for the custom gene and cell type
    [T1, T2] = sce_circ_phase_estimation_stattest(sce, tmeta, false, period12, ...
                                          cust_gene, cust_cells, false, norm_str);
    
    % Find the specified gene in the results
    gene_idx = find(strcmp(T1.Genes, cust_gene));
    if isempty(gene_idx)
        error('Specified gene not found in the dataset.');
    end

    if strcmp(norm_str, 'lib_size')
        % This is regular cells pipeline
        %X = sc_norm(full(sce.X));
        X = pkg.norm_libsize(sce.X, 1e4);
        X = log1p(sce.X);
    else % 'magic_impute'
        % This for cancer cells
        X = sc_impute(sce.X, 'MAGIC');
    end
    sce.X = X;
    clear X;

    % Subset and normalize cells
    ic0 = find(sce.c_cell_type_tx == cust_cells);
    sce_sub = sce.selectcells(ic0);

    ig = find(sce_sub.g == cust_gene);

    clear ic0;

    % Prepare gene expression per time point
    batch_time = unique(sce_sub.c_batch_id);
    nzts = length(batch_time);

    Xvals = zeros( length(sce_sub.c_batch_id), 1);
    tvals = Xvals;
    ncell_cummu = 0;

    % Array for the mean expression values from sce_sub at each time point
    for it = 1:nzts
        ics = find(sce_sub.c_batch_id == batch_time(it));
        nc_loc = length(ics);
        ibeg = ncell_cummu + 1;
        iend = ncell_cummu + nc_loc;
        
        % Ensure we're assigning within valid array boundaries
        if ig <= size(sce_sub.X, 1)
            Xvals(ibeg:iend) = full(sce_sub.X(ig, ics)); % Assign expression data
            tvals(ibeg:iend) = t(it); % Assign actual time from `t`
        end
        
        ncell_cummu = ncell_cummu + nc_loc;
    end
    clear sce_sub;
    
    % Generate the plot for the specified gene if found
    % Set the current axes to the provided axHandle
    axes(axHandle);

    % Generate sine-fitted values based on circadian parameters
    fval = T1.Amp(gene_idx) * cos(2 * pi * (tval - T1.Acrophase(gene_idx)) / T1.Period(gene_idx)) + T1.Mesor(gene_idx);
    Rzts = table2array(T2(gene_idx, 2:end));

    % Plot sine-fitted values
    plot(tval, fval, 'LineWidth', 2);
    hold on;

    % Plot original expression data (Rzts)
    plot(t, Rzts, 'o-', 'LineWidth', 1.5);

    if print_scdata
        % Plot scatter points with transparency to highlight overlaps
        scatter(tvals, Xvals, 20, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', ...
            'MarkerFaceAlpha', 0.12, 'MarkerEdgeAlpha', 0.12); % Set alpha for transparency
    end

    % Plot vertical line for phase
    xline(T1.Acrophase_24(gene_idx), 'LineStyle', '--', 'Color', 'r');

    % Set plot limits and labels
    xlim([t0 tf]);
    % Modify the x-axis ticks to show every 1 hour
    xticks(t0:tint:tf);
    % Titles and legends
    title_str = sprintf('Gene - %s | p-value: %.3f | Phase: %.3f | NumCells: %d' ,...
                      T1.Genes{gene_idx}, T1.pvalue(gene_idx), ...
                      T1.Acrophase_24(gene_idx), ncell_cummu );
    title(title_str, 'FontSize', 16); % Increase title font size
    xlabel('Time (hrs)', 'FontSize', 14); % Increase x-axis label font size
    ylabel('Expression', 'FontSize', 14); % Increase y-axis label font size
    set(gca, 'FontSize', 12); % Increase axis numbers font size
    if print_scdata
        legend({'Sine-fitted expression', 'Mean Expression (Rzts)', ...
            'sc-data', 'Phase'}, 'Location', 'northwest');
    else
        legend({'Sine-fitted expression', 'Mean Expression (Rzts)', ...
            'Phase'}, 'Location', 'northwest');
    end
    hold off;
end
