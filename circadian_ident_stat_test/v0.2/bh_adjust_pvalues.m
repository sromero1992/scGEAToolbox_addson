function [p_adj] = bh_adjust_pvalues(pvals)
    % Function to compute Benjamini-Hochberg (BH) adjusted p-values
    % INPUT:
    % pvals => Vector of p-values from F-test
    % OUTPUT:
    % p_adj => Vector of adjusted p-values
    
    % Number of tests
    m = length(pvals);
    
    % Sort p-values in ascending order and keep track of original indices
    [sorted_pvals, sort_idx] = sort(pvals);
    
    % Initialize adjusted p-values
    p_adj = zeros(m,1);
    
    % Compute adjusted p-values using BH procedure
    for i = 1:m
        % Compute the BH adjusted p-value
        p_adj(i) = sorted_pvals(i) * m / i;
    end
    
    % Ensure that the adjusted p-values do not exceed 1
    p_adj = min(p_adj, 1);
    
    % Reorder adjusted p-values to match original order
    p_adj = p_adj(sort_idx);
end
