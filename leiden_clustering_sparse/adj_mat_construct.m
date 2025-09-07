function adjX = adj_mat_construct_sparse_blocked_fixed(sce, method, K, chunk_size, use_hvgs)
    % INPUT:
    % sce --------> SCE object
    % method -----> Neighbor method ('knn' or 'mnn')
    % K ----------> Number of neighbors
    % chunk_size -> Number of rows to process in each chunk
    % OUTPUT:
    % adjX -------> Sparse adjacency matrix
    % AUTHOR: Selim Romero, Texas A&M University
    % Input validation
    method = lower(method);
    if ~ismember(method, {'knn', 'mnn'})
        error('Method must be either ''knn'' or ''mnn''.');
    end
    if nargin < 3; K = 15; end
    if ~isnumeric(K) || K <= 0 || K ~= round(K)
        error('K must be a positive integer.');
    end
    if nargin < 4; chunk_size = 1e4; end
    if ~isnumeric(chunk_size) || chunk_size <= 0 || chunk_size ~= round(chunk_size)
        error('chunk_size must be a positive integer.');
    end
    fprintf("Computing adjacency matrix! \n");
    tic;
    % Normalize and preprocess input data (genes by cells)
    X = sce.X;
    X = sc_norm(X, 'type', 'libsize');
    if use_hvgs
        [~, X] = sc_splinefit(X, sce.g);
        X = X(1:nhvgs, :);
    end
    X = log1p(X)'; % (cells by genes)
    [U, ~, ~] = svds(X', 50); % (cells by meta-genes)
    X = X * U;
    fprintf("Time for preparing data: %f \n", toc);
    num_cells = size(X, 1);
    
    % Step 1: Compute KNN for all cells
    fprintf("Finding %d nearest neighbors for all cells...\n", K);
    tic;
    % Using a searcher object is more memory-efficient
    Mdl = KDTreeSearcher(X);
    [Idx, ~] = knnsearch(Mdl, X, 'K', K + 1);
    Idx = Idx(:, 2:end); % Remove self-matches
    fprintf("Time for finding neighbors: %f \n", toc);
    
    % Step 2: Construct the adjacency matrix based on the chosen method
    adjX = sparse(num_cells, num_cells);
    switch method
        case 'knn'
            % Construct KNN adjacency matrix in a single loop
            for i = 1:num_cells
                 adjX(i, Idx(i, :)) = 1;
            end
        case 'mnn'
            % Construct MNN adjacency matrix in a chunked loop
            fprintf("Constructing MNN adjacency matrix in chunks...\n");
            tic;
            num_chunks = ceil(num_cells / chunk_size);
            for chunk_idx = 1:num_chunks
                fprintf("Processing chunk %d of %d...\n", chunk_idx, num_chunks);
                ibeg = (chunk_idx - 1) * chunk_size + 1;
                iend = min(chunk_idx * chunk_size, num_cells);
                
                for i = ibeg:iend
                    neighbors_i = Idx(i, :);
                    for j = neighbors_i
                        % Check if i is a neighbor of j using pre-computed results
                        if ismember(i, Idx(j, :))
                            adjX(i, j) = 1;
                            adjX(j, i) = 1;
                        end
                    end
                end
            end
            fprintf("Time for MNN construction: %f \n", toc);
    end
    
    % Ensure the adjacency matrix is symmetric
    adjX = max(adjX, adjX');
    fprintf("Total time for neighbors: %f \n", toc);
    % Final adjacency matrix is sparse
    adjX = sparse(adjX);
end