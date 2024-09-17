% Specify gene of interest
igene = 'Msh2';
idx = find(igene == sce.g);  % Assuming sce and idx are defined

% Use raw counts for fitting
raw_gene_count = sce.X(idx, :)';  % Raw gene counts

% Set different thresholds to compare
thresholds = [0, 1, 2, 5, 10, 50];  % Adjust threshold values as needed

% Initialize arrays to store results
aic_values_nb = zeros(length(thresholds), 1);
p_values_nb = zeros(length(thresholds), 1);
aic_values_zinb = zeros(length(thresholds), 1);
p_values_zinb = zeros(length(thresholds), 1);

for i = 1:length(thresholds)
    thrs = thresholds(i);  % Current threshold
    filt_gene_count = raw_gene_count(raw_gene_count >= thrs);  % Filtered gene counts
    
    % --- Fit a Negative Binomial Distribution ---
    try
        pd_nb = fitdist(filt_gene_count, 'NegativeBinomial');
    catch ME
        warning('Error fitting Negative Binomial: %s', ME.message);
        continue;
    end
    
    % Calculate PMF for a range of gene counts
    x = min(filt_gene_count):max(filt_gene_count);  % Range of gene counts
    nbpdf = nbinpdf(x, pd_nb.R, pd_nb.p);  % PMF values
    
    % Goodness-of-Fit Test: Chi-Squared for Negative Binomial
    [h_nb, p_nb] = chi2gof(filt_gene_count, 'CDF', @(x) nbincdf(x, pd_nb.R, pd_nb.p));
    p_values_nb(i) = p_nb;  % Store p-value
    
    % --- AIC Calculation for Negative Binomial ---
    nbpdf_full = nbinpdf(filt_gene_count, pd_nb.R, pd_nb.p);
    logL_full = sum(log(nbpdf_full + eps));  % Adding eps to avoid log(0)
    aic_values_nb(i) = -2 * logL_full + 2 * pd_nb.NumParameters;
    
    % --- Zero-Inflated Negative Binomial (ZINB) Model Fit ---
    % Prepare data for zero-inflation modeling
    zero_ind = raw_gene_count == 0;
    zero_inflation_model = fitglm(raw_gene_count, zero_ind, 'Distribution', 'binomial', 'Link', 'logit');
    
    % Fit Negative Binomial distribution to non-zero counts
    try
        pd_nb_zinb = fitdist(filt_gene_count, 'NegativeBinomial');
    catch ME
        warning('Error fitting Negative Binomial for ZINB model: %s', ME.message);
        continue;
    end
    
    % Ensure the range is consistent for plotting and analysis
    x_range = min(filt_gene_count):max(filt_gene_count);
    
    % Calculate PMF for zero-inflated model
    zero_inflated_pdf = zeros(size(x_range));  % Initialize zero-inflated PMF
    nbpdf_zinb = nbinpdf(x_range, pd_nb_zinb.R, pd_nb_zinb.p);  % PMF for non-zero counts
    zero_prob = predict(zero_inflation_model, x_range');
    
    for j = 1:length(x_range)
        if x_range(j) == 0
            zero_inflated_pdf(j) = zero_prob(j);
        else
            zero_inflated_pdf(j) = (1 - zero_prob(j)) * nbpdf_zinb(j);
        end
    end
    
    % Ensure zero_inflated_pdf is normalized
    zero_inflated_pdf = zero_inflated_pdf / sum(zero_inflated_pdf);
    
    % Calculate the CDF from the PDF
    zero_inflated_cdf = cumsum(zero_inflated_pdf);  % Compute cumulative sum for CDF
    
    % Define the CDF function for chi2gof
    zero_inflated_cdf_func = @(x) interp1(x_range, zero_inflated_cdf, x, 'linear', 'extrap');
    
    % Perform the goodness-of-fit test using the CDF
    try
        [h_zinb, p_zinb] = chi2gof(filt_gene_count, 'CDF', zero_inflated_cdf_func);
        p_values_zinb(i) = p_zinb;  % Store p-value
    catch ME
        warning('Error in chi-squared goodness-of-fit test for ZINB: %s', ME.message);
        p_values_zinb(i) = NaN;
    end

    % Plot histogram and fitted PMFs
    figure;
    histogram(filt_gene_count, 'Normalization', 'pdf');
    hold on;
    plot(x, nbpdf, 'r-', 'LineWidth', 2);
    plot(x_range, zero_inflated_pdf, 'b-', 'LineWidth', 2);  % Use consistent range
    title(['Histogram, NB and ZINB PMFs for ', igene, ' with threshold ', num2str(thrs)]);
    xlabel('Gene Counts');
    ylabel('Probability Density');
    legend('Histogram', 'Negative Binomial PMF', 'Zero-Inflated Negative Binomial PMF');
    hold off;
end

% Display AIC values and p-values
fprintf('Threshold\tAIC (NB)\tp-value (NB)\tAIC (ZINB)\tp-value (ZINB)\n');
for i = 1:length(thresholds)
    fprintf('%d\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n', thresholds(i), aic_values_nb(i), p_values_nb(i), aic_values_zinb(i), p_values_zinb(i));
end

