function mat = build_corr_xi2(A, B)
    % Build correlation matrix for new correlation coefficient
    % INPUT:
    % A ========> Matrix containing features by observations
    % B ========> Matrix containing features by observations
    %             B matrix is optional for being different matrices
    % use_lite==> Boolean flag to use optimized "lite" version
    % OUTPUT:
    % mat ======> Matrix containing new correlation matrix elements ij
    %             for cross features Ai x Bj (i-th feature and j-th feat)

    if nargin < 2; B = A; end  % If only one input matrix, use A for B as well

    tic;  % Start timing

    % Transpose A for efficient column-wise access
    A = A';
    B = B';

    % Number of columns are the number of features
    n_Afeat = size(A, 2);
    n_Bfeat = size(B, 2);
    fprintf("Number of features in matrix A %d and B %d \n", n_Afeat, n_Bfeat);

    % Defining size of correlation matrix
    mat = zeros(n_Afeat, n_Bfeat);

    % Using parfor to parallelize the computation for each feature of A
    parfor ai_feat = 1:n_Afeat
        % Temporary vector to store the correlation results for this iteration
        tmp_vec = zeros(1, n_Bfeat);

        % Sorting A for the current feature (column)
        [~, idx_sort] = sort(A(:, ai_feat), 'ascend');  % Sort current column (feature) of A
        
        % Loop over the features in B
        for bi_feat = 1:n_Bfeat
            % Calculate new correlation coefficient using the provided function
            % Use idx_sort from A and B column
            tmp_vec(bi_feat) = correlation_xi([], B(:, bi_feat), idx_sort);
        end
        
        % Gather results into the main matrix
        mat(ai_feat, :) = tmp_vec;
    end

    time_corr = toc;  % End timing
    fprintf('Computational time of xi correlation: %.4f seconds\n', time_corr);  % Display computational time
end
