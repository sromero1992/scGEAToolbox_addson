% Build probability density distribution (PDF) for Msh2 and tdTomato
idx = find(igene == sce.g);
max_count = 101;
counts = zeros(max_count,1);
for icount = 1: max_count
    counts(icount) = sum( sce.X(idx,:) == icount-1);
    if icount == max_count
        counts(icount) = sum( sce.X(idx,:) >= icount-1);
    end
end
plot(log1p(counts))




% Initial threshold settings and gene selection
thrs0 = 0;
thrs1 = log1p(0.99);
thrs2 = log1p(1.99);
thrs = 100; % For raw counts

% Specify gene of interest
igene = 'tdTomato';
idx = find(igene == sce.g);

% Use raw counts for fitting
raw_gene_count = sce.X(idx, :)';  % Assuming `sce` and `idx` are defined earlier

% Filter out counts based on a threshold
filt_gene_count = raw_gene_count(raw_gene_count >= thrs);  % Apply threshold filtering on raw counts

% Plot histogram of the raw gene counts after filtering
figure;
histogram(filt_gene_count, 51, 'Normalization', 'pdf'); % Normalized histogram to match PMF
hold on;

% --- Fit a Negative Binomial Distribution ---
% Fit Negative Binomial distribution to the full data
pd_full = fitdist(filt_gene_count, 'NegativeBinomial');

% Calculate PMF for a range of gene counts
x = min(filt_gene_count):max(filt_gene_count); % Range of gene counts
nbpdf_full = nbinpdf(x, pd_full.R, pd_full.p); % PMF values

% Plot the PMF
plot(x, nbpdf_full, 'r-', 'LineWidth', 2);

% Add labels and title
title(['Histogram and Negative Binomial PMF for ' igene]);
xlabel('Gene Counts');
ylabel('Probability Density');

% Show legend
legend('Histogram', 'Negative Binomial PMF');
hold off;

% --- AIC and Goodness-of-Fit Tests ---
% Fit models for trimmed data
% Trim 0s (if count == 0) and fit Negative Binomial model
gene_count_rm_0s = filt_gene_count(filt_gene_count > 0);
pd_rm_0s = fitdist(gene_count_rm_0s, 'NegativeBinomial');

% Trim 0s and 1s and fit Negative Binomial model
gene_count_rm_0s_and_1s = filt_gene_count(filt_gene_count > thrs);
pd_rm_0s_and_1s = fitdist(gene_count_rm_0s_and_1s, 'NegativeBinomial');

% Calculate log-likelihood for each model
nbpdf_full = nbinpdf(filt_gene_count, pd_full.R, pd_full.p);
logL_full = sum(log(nbpdf_full));

nbpdf_rm_0s = nbinpdf(gene_count_rm_0s, pd_rm_0s.R, pd_rm_0s.p);
logL_rm_0s = sum(log(nbpdf_rm_0s));

nbpdf_rm_0s_and_1s = nbinpdf(gene_count_rm_0s_and_1s, pd_rm_0s_and_1s.R, pd_rm_0s_and_1s.p);
logL_rm_0s_and_1s = sum(log(nbpdf_rm_0s_and_1s));

% Calculate AIC for each model
aic_full = -2 * logL_full + 2 * pd_full.NumParameters;
aic_rm_0s = -2 * logL_rm_0s + 2 * pd_rm_0s.NumParameters;
aic_rm_0s_and_1s = -2 * logL_rm_0s_and_1s + 2 * pd_rm_0s_and_1s.NumParameters;

% Display AIC values
fprintf('AIC for full data: %.4f\n', aic_full);
fprintf('AIC for trimmed 0s: %.4f\n', aic_rm_0s);
fprintf('AIC for trimmed 0s and 1s: %.4f\n', aic_rm_0s_and_1s);

% Goodness-of-Fit Test
% Ensure integer data for goodness-of-fit test
filt_gene_count = round(filt_gene_count);
gene_count_rm_0s = round(gene_count_rm_0s);
gene_count_rm_0s_and_1s = round(gene_count_rm_0s_and_1s);

% Goodness-of-Fit Test
try
    [h_full, p_full] = chi2gof(filt_gene_count, 'CDF', @(x) nbincdf(x, pd_full.R, pd_full.p));
    [h_trimmed_0s, p_trimmed_0s] = chi2gof(gene_count_rm_0s, 'CDF', @(x) nbincdf(x, pd_rm_0s.R, pd_rm_0s.p));
    [h_trimmed_0s_and_1s, p_trimmed_0s_and_1s] = chi2gof(gene_count_rm_0s_and_1s, 'CDF', @(x) nbincdf(x, pd_rm_0s_and_1s.R, pd_rm_0s_and_1s.p));
catch ME
    warning('Error in chi-squared goodness-of-fit test: %s', ME.message);
    h_full = NaN; p_full = NaN;
    h_trimmed_0s = NaN; p_trimmed_0s = NaN;
    h_trimmed_0s_and_1s = NaN; p_trimmed_0s_and_1s = NaN;
end

% Display goodness-of-fit test results
fprintf('Goodness-of-fit p-value for full data: %.4f\n', p_full);
fprintf('Goodness-of-fit p-value for trimmed 0s: %.4f\n', p_trimmed_0s);
fprintf('Goodness-of-fit p-value for trimmed 0s and 1s: %.4f\n', p_trimmed_0s_and_1s);

% Plot histogram of filtered gene counts after trimming strategy that fits best
figure;
if aic_rm_0s < aic_full && aic_rm_0s < aic_rm_0s_and_1s
    % Trimmed 0s is the best model
    histogram(gene_count_rm_0s, 30, 'Normalization', 'pdf');
    title(['Filtered Gene Counts for ' igene ' (Trimmed 0s)']);
    xlabel('Gene Counts');
    ylabel('Probability Density');
elseif aic_rm_0s_and_1s < aic_full && aic_rm_0s_and_1s < aic_rm_0s
    % Trimmed 0s and 1s is the best model
    histogram(gene_count_rm_0s_and_1s, 30, 'Normalization', 'pdf');
    title(['Filtered Gene Counts for ' igene ' (Trimmed 0s and 1s)']);
    xlabel('Gene Counts');
    ylabel('Probability Density');
else
    % Full data is the best model
    histogram(filt_gene_count, 30, 'Normalization', 'pdf');
    title(['Filtered Gene Counts for ' igene ' (Full Data)']);
    xlabel('Gene Counts');
    ylabel('Probability Density');
end
