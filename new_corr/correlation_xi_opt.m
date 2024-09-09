function correlation_coefficient = correlation_xi_opt(X, Y, mode, precomputed_ranks)
    % Computes Spearman rank correlation coefficient for conditional dependency
    % https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1758115
    %
    % INPUT: 
    % X: Independent variable
    % Y: Probably dependent variable
    % mode: mode=1 computes for ties, mode=2 computes for no ties
    % precomputed_ranks (optional): Pre-computed ranks of X for repeated use
    %                               (useful for matrix computations)
    % OUTPUT:
    % correlation_coefficient: Degree of dependence between X and Y
    %                          (0 for independence, 1 for functional relationship)

    if nargin < 4
        precomputed_ranks = [];
    end

    if isempty(precomputed_ranks)
        % Re-arrange X and Y according to X values
        [~, precomputed_ranks] = sort(X, 'ascend');
    end
    Y = Y(precomputed_ranks);
    
    % Number of observations in Y or X
    n = length(Y);
    
    % Calculate ranks of Y using tiedrank for efficiency
    r_ranks = tiedrank(Y);
    
    % Precompute rank differences
    r_plus_minus_diff = abs( diff(r_ranks) );  % r_plus - r_minus

    if mode == 2  % No ties
        % Calculate the correlation coefficient without ties
        correlation_coefficient = 1 - (3 * sum(r_plus_minus_diff)) / (n^2 - 1);
    else  % Ties present
        % Calculate ranks and counts for ties
        l_ranks = n + 1 - r_ranks;  % Reverse rank for ties (l(i) = n - r(i) + 1)
        
        % Calculate the correlation coefficient with ties
        numerator = n * sum(r_plus_minus_diff);
        denominator = 2 * sum(l_ranks .* (n - l_ranks));
        correlation_coefficient = 1 - numerator / denominator;
    end
end