function  s = sc_louvain_clustering_cells(sce0)
    % Louvain clustering from sce object
    % Classification is for cells

    % Get full count matrix representation
    X = full(sce0.X);
    % Library size normalization
    X = sc_norm(X);
    % Log1p transformation
    X = log1p(X);

    fprintf("Shape %d %d \n", size(X));

    % Cell by genes represenation
    X = X';

    fprintf('Creating similarity matrix...\n');
    % Creating similarity matrix
    % cosineSim = @(x,y) dot(x,y)/(norm(x)*norm(y));
    % S = zeros(size(X, 1));
    % for i = 1:size(X, 1)
    %     for j = 1:size(X, 1)
    %         S(i, j) = cosineSim(X(i, :), X(j, :));
    %     end
    % end
    % Normalize the rows to unit length
    row_norms = sqrt(sum(X.^2, 2));
    X_normalized = X ./ row_norms;

    % Compute the cosine similarity matrix using matrix multiplication
    S = X_normalized * X_normalized';
    clear X_normalized;

    % Symmetrize
    S = (S + S') / 2;

    fprintf('Getting k-NN graph... \n');
    % Construct the adjacency matrix (k-NN graph)
    k = 10; % 10 neighboars
    A = zeros(size(S));
    for i = 1:size(S, 1)
        [~, idx] = sort(S(i, :), 'descend');
        A(i, idx(2:k+1)) = S(i, idx(2:k+1));
    end
    A = (A + A') / 2;

    fprintf('Iterative Louvain... \n');
    % Step 4: Run Louvain clustering
    gamma = 4.0;
    [s, ~] = iterated_genlouvain(A,[],[],gamma);
    %[s, ~] = genlouvain(A);

end