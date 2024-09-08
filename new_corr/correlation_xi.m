function correlation_coefficient = correlation_xi(X, Y, mode, precomputed_ranks)
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
    %                           (0 for independence, 1 for functional relationship)
    
    if nargin < 4
        precomputed_ranks = [];
    end

    if isempty(precomputed_ranks)
        % Re-arrange X and Y according to X values
        [~, precomputed_ranks]  = sort(X, 'ascend');
    end
    Y = Y(precomputed_ranks);
    
    % Number of observations in Y or X
    n = length(Y);
    
    % Calculate ranks of Y
    r_ranks = zeros(n, 1);
    for i = 1:n
        r_ranks(i) = sum(Y(:) <= Y(i));
    end
    
    % This are already n-1 elements
    r_plus =  r_ranks(2:end);
    r_minus = r_ranks(1:end-1);
    
    if mode == 2 % No ties       
        % Calculate the correlation coefficient
        correlation_coefficient = 1 - 3 * sum( abs(r_plus - r_minus) ) / (n^2 - 1);
    else % Ties 
        % Handle ties here (using the more complex formula)
        % Calculate ranks and counts for ties
        l_ranks = zeros(n, 1);
        for i = 1:length(Y)
            l_ranks(i) = sum(Y(:) >= Y(i));
        end
        % Calculate the correlation coefficient
        numerator = n * sum( abs(r_plus - r_minus) );
        denominator = 2 * sum( l_ranks * ( n - l_ranks) ) ;
        correlation_coefficient = 1 - numerator / denominator;
    end
end