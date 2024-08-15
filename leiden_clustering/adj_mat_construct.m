function adjX = adj_mat_construct(sce, method, K)
    % INPUT:
    % sce -------> SCE object 
    % method ----> Neighbor method (knn or mnn)
    % K ---------> Number of neighbors
    % OUTPUT:
    % adjX ------> Adjacency matrix from cosine similarity and neighbors
    % USAGE:
    % 

    % Input validation
    method = lower(method);
    if ~ismember(method, {'knn', 'mnn'})
        error('Method must be either ''knn'' or ''mnn''.');
    end
    if ~isnumeric(K) || K <= 0 || K ~= round(K)
        error('K must be a positive integer.');
    end

    tic;
    % X in cells by genes basis and normalize/scale
    X = full(sce.X); 
    X = sc_norm(X,'type','libsize');
    X = log1p(X)';

    % Perform PCA up to top 50 components
    %X = svdpca(X, 50, 'random'); % Actually faster but looks different...
    %X = svdpca(X, 50, 'svd');
    [U, ~, ~] = svds(X', 50);
    X = X * U;
    time_prep = toc;
    fprintf("Time for preparing data: %f \n", time_prep);

    % Compute pairwise cosine similarity
    tic;
    simX = 1 - pdist2(X, X, 'cosine');
    time_sim = toc;
    fprintf("Time for similarity: %f \n", time_sim);
    
    % Compute pairwise Jaccard similarity (for binary data)
    %simX = 1 - pdist2(X > 0, X > 0, 'jaccard');

    tic;
    % Define adjacency matrix
    sim_size = size( simX);
    adjX = zeros(sim_size);

    switch method
        case 'knn'
            % K-Nearest Neighbors (KNN):
            for i = 1:sim_size(1)       
                [~, idx] = sort( simX(i, :), 'descend');
                adjX(i, idx(1:K)) = 1;
            end
            % Ensure the adjacency matrix is symmetric
            adjX = max(adjX, adjX');
    
        case 'mnn'
            % Mutual Nearest Neighbors (MNN):
            for i = 1:sim_size(1)
                [~, idx] = sort( simX(i, :), 'descend');
                % Ith-neighbors
                neighbors_i = idx(1:K);
                for j = neighbors_i
                    [~, idx_j] = sort( simX(j, :), 'descend');
                    % Jth-neighbor
                    neighbors_j = idx_j(1:K);
                    if ismember(i, neighbors_j)
                        adjX(i, j) = 1;
                        adjX(j, i) = 1;
                    end
                end
            end
    end
    time_nei = toc;
    fprintf("Time for neighbors: %f \n", time_nei);
end
