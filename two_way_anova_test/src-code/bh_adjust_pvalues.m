function [p_adj] = bh_adjust_pvalues(pvals)
    % Benjamini-Hochberg (BH) adjustment with NaN handling
    % INPUT:
    % pvals => Vector of p-values (can include NaNs)
    % OUTPUT:
    % p_adj => Vector of adjusted p-values (NaN for invalid input)

    % Ensure input is a column vector
    pvals = pvals(:);

    % Initialize output with NaNs
    p_adj = NaN(size(pvals));

    % Identify valid p-values (non-NaN and in range [0, 1])
    valid_idx = find(~isnan(pvals) & pvals >= 0 & pvals <= 1);

    % Extract only valid p-values
    valid_pvals = pvals(valid_idx);

    % Number of valid tests
    m = length(valid_pvals);

    if m > 0
        % Sort valid p-values and keep track of original indices
        [sorted_pvals, sort_idx] = sort(valid_pvals);
        sort_idx2 = valid_idx(sort_idx);

        % Compute BH adjusted p-values
        p_adj_sorted = sorted_pvals .* (m ./ (1:m)'); % Vectorized adjustment

        % Ensure monotonicity (non-decreasing adjusted p-values)
        p_adj_sorted = min(cummin(p_adj_sorted(end:-1:1)), 1); % Reverse and apply cummin
        p_adj_sorted = p_adj_sorted(end:-1:1); % Reverse back

        % Map adjusted p-values back to original indices
        p_adj(sort_idx2) = p_adj_sorted;
        
    end
end
