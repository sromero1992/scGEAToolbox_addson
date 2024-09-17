function correlation_coefficient = correlation_xi(X, Y, precomputed_ranks)
    % Computes rank correlation coefficient for conditional dependency
    % https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1758115
    %
    % INPUT: 
    % X ============> Independent variable
    % Y ============> Probably dependent variable
    % mode =========> mode=2 computes for no ties, mode=1 computes for ties
    % precomputed_ranks (optional): Pre-computed ranks of X for repeated use
    %                               (useful for matrix computations)
    % use_lite =====> Boolean flag to use optimized "lite" version
    %
    % OUTPUT:
    % correlation_coefficient: Degree of dependence between X and Y
    %                          (0 for independence, 1 for functional relationship)
    
    % Set default arguments
    if nargin < 3; precomputed_ranks = []; end % Default to sort X
    
    if isempty(precomputed_ranks) && isempty(X)
        error("X is void and no precomputed ranks provided.");
    end

    % Compute ranks for X if not precomputed
    if isempty(precomputed_ranks) || ~isempty(X)
        [~, precomputed_ranks] = sort(X, 'ascend');
    end

    % Compute both forward and backward rankings
    [forward_yrank, backward_yrank] = rank_max_ties(Y);

    %n = length(Y);
    %forward_yrank = tiedrank(Y);
    %backward_yrank = n + 1 - forward_yrank;  % Reverse rank for ties (l(i) = n - r(i) + 1)

    Y = Y(precomputed_ranks);
    forward_yrank = forward_yrank(precomputed_ranks);
    backward_yrank = backward_yrank(precomputed_ranks);

    % Number of observations
    n = length(Y);
    
    % Precompute rank differences
    r_plus_minus_diff = abs( diff(forward_yrank) );  % r_plus - r_minus
    
    % Compute correlation coefficient based on mode
    % Ties present
    % Calculate the correlation coefficient with ties
    numerator = n * sum(r_plus_minus_diff);
    denominator = 2 * sum( backward_yrank .* (n - backward_yrank) );
    correlation_coefficient = 1 - numerator / denominator;
    % % No ties present
    % % Calculate the correlation coefficient without ties
    %correlation_coefficient = 1 - (3 * sum(r_plus_minus_diff)) / (n^2 - 1);
end