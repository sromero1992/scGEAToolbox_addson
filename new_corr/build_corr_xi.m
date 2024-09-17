function mat = build_corr_xi(A, B, use_lite)
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
    if nargin < 3; use_lite = true; end

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
    % Compute the correlation with each feature in B
    vec = zeros(n_Bfeat, 1);

    % Using parfor to parallelize the computation for each feature of A
    for ai_feat = 1:n_Afeat
        % Sorting A for current feature (column)
        [~, idx_sort] = sort(A(:, ai_feat), 'ascend');  % Sort current column (feature) of A
        % A sorted is not needed but only for sorting Y
        % A_sorted = A(idx_sort, ai_feat);  % Sorted observations for the current feature of A

        parfor bi_feat = 1:n_Bfeat
            % Calculate new correlation coefficient using the provided function
            % vec(bi_feat) = correlation_xi(A(:, ai_feat), B(:, bi_feat), 2, idx_sort, true);
            vec(bi_feat) = correlation_xi([], B(:, bi_feat), idx_sort);
        end
        mat(ai_feat, 1:n_Bfeat) = vec(1:n_Bfeat);
    end
    time_corr = toc;  % End timing
    fprintf('Computational time of xi correlation: %.4f seconds\n', time_corr);  % Display computational time
    %mat_final = mat(:,:);  % Return the final correlation matrix
end