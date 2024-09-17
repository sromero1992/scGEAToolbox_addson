function nb_test_per_gene(sce, igene, thresholds, plotit)
    % nb_test_per_gene fits a Negative Binomial distribution to the gene counts 
    % and computes AIC and p-values.
    %
    % INPUT:
    %   - sce: SingleCellExperiment object containing gene expression matrix and metadata.
    %   - igene: The name of the gene to be tested.
    %   - thresholds: An array of thresholds for filtering gene counts.
    %   - plotit: Boolean to decide if plots should be generated.
    %
    % OUTPUT:
    %   None. AIC and p-values are printed to the console.

    if nargin < 4; plotit = false; end
    
    % Get the gene index and extract gene expression data
    idx = find(strcmp(igene, sce.g));  
    raw_gene_count = sce.X(idx, :)';  % Raw gene counts

    % Initialize arrays for AIC and p-values
    aic_values_nb = zeros(length(thresholds), 1);
    p_values_nb = zeros(length(thresholds), 1);

    for i = 1:length(thresholds)
        thrs = thresholds(i);
        filt_gene_count = raw_gene_count(raw_gene_count >= thrs);

        if length(filt_gene_count) < 2 || numel(unique(filt_gene_count)) < 2
            warning(['Not enough data to fit a Negative Binomial for threshold ', num2str(thrs)]);
            aic_values_nb(i) = NaN;
            p_values_nb(i) = NaN;
            continue;
        end
        
        % Fit Negative Binomial distribution
        pd_nb = fitdist(filt_gene_count, 'NegativeBinomial');

        % Perform Chi-Square Goodness-of-Fit Test
        [~, p_nb] = chi2gof(filt_gene_count, 'CDF', @(x) nbincdf(x, pd_nb.R, pd_nb.p));
        p_values_nb(i) = p_nb;

        % Compute AIC
        nbpdf_full = nbinpdf(filt_gene_count, pd_nb.R, pd_nb.p);
        logL_full = sum(log(nbpdf_full + eps));  
        aic_values_nb(i) = -2 * logL_full + 2 * pd_nb.NumParameters;
        
        % Optionally, plot the results
        if plotit
            figure;
            histogram(filt_gene_count, 'Normalization', 'pdf');
            hold on;
            x = min(filt_gene_count):max(filt_gene_count);
            plot(x, nbinpdf(x, pd_nb.R, pd_nb.p), 'r-', 'LineWidth', 2);
            title(['Histogram and NB PMF for ', igene, ' with threshold ', num2str(thrs)]);
            xlabel('Gene Counts');
            ylabel('Probability Density');
            legend('Histogram', 'Negative Binomial PMF');
            hold off;
        end
    end
    
    % Display AIC and p-values
    fprintf('Threshold\tAIC (NB)\tp-value (NB)\n');
    for i = 1:length(thresholds)
        fprintf('%d\t\t%.4f\t\t%.4f\n', thresholds(i), aic_values_nb(i), p_values_nb(i));
    end
end

