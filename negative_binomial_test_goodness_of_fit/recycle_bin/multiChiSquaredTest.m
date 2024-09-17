function [Chi2, pVal] = multiChiSquaredTest(filt_gene_count1, filt_gene_count2, copula_params, nBins)
    % Compute copula variables
    pd_nb1 = fitdist(filt_gene_count1, 'NegativeBinomial');
    pd_nb2 = fitdist(filt_gene_count2, 'NegativeBinomial');
    u1 = nbincdf(filt_gene_count1, pd_nb1.R, pd_nb1.p);
    u2 = nbincdf(filt_gene_count2, pd_nb2.R, pd_nb2.p);
    
    % Ensure copula inputs are strictly between 0 and 1
    u1 = min(max(u1, eps), 1-eps);
    u2 = min(max(u2, eps), 1-eps);
    
    % Create intervals for binning
    edges = linspace(0, 1, nBins + 1);  % Binning for [0,1] range
    
    % Generate all combinations of bin indices
    [bin1, bin2] = ndgrid(1:nBins, 1:nBins);
    bins = [bin1(:), bin2(:)];
    
    % Compute observed frequencies
    O = zeros(size(bins, 1), 1);
    for j = 1:size(bins, 1)
        bin_x = bins(j, 1);
        bin_y = bins(j, 2);
        O(j) = sum(u1 > edges(bin_x) & u1 <= edges(bin_x + 1) & ...
                   u2 > edges(bin_y) & u2 <= edges(bin_y + 1));
    end
    O = O / length(u1);  % Normalize observed frequencies
    
    % Compute expected frequencies under the fitted copula model
    E = zeros(size(bins, 1), 1);
    for j = 1:size(bins, 1)
        bin_x = bins(j, 1);
        bin_y = bins(j, 2);
        lb = [edges(bin_x), edges(bin_y)];
        ub = [edges(bin_x + 1), edges(bin_y + 1)];
        E(j) = copulacdf('Gaussian', ub, copula_params) - ...
               copulacdf('Gaussian', [ub(1), lb(2)], copula_params) - ...
               copulacdf('Gaussian', [lb(1), ub(2)], copula_params) + ...
               copulacdf('Gaussian', lb, copula_params);
    end
    E = max(E, eps);  % Avoid division by zero

    % Perform Chi-squared test
    Chi2 = sum(((O - E).^2) ./ E);
    df = length(O) - 1;  % Degrees of freedom
    pVal = 1 - chi2cdf(Chi2, df);  % p-value for the Chi-squared test

    % Display results
    fprintf('Chi-squared Statistic: %.4f\n', Chi2);
    fprintf('p-value: %.4f\n', pVal);
end