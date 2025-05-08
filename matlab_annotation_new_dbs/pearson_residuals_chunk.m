function z_out = pearson_residuals_chunk(X, theta, chunk_size_cells)
    %{
    Computes Pearson residuals for a count matrix X, processing in cell-wise chunks.

    The Pearson residual z_gc for gene g and cell c is computed as:
    z_gc = (X_gc - mu_gc) / sigma_gc
    (See pearson_residuals_gene_chunks for full formula details)

    INPUT:
    X                -----> count matrix (genes x cells), can be sparse or full.
    theta            -----> dispersion parameter for the NB variance (default: 100).
    chunk_size_cells -----> Number of cells to process in each chunk (default: 10000).

    OUTPUT:
    z_out            -----> Pearson residual matrix (genes x cells), typically dense.
    %}

    if nargin < 2 || isempty(theta)
        theta = 100;
    end
    if nargin < 3 || isempty(chunk_size_cells)
        chunk_size_cells = 10000; % Default chunk size for cells
    end

    [g, c] = size(X);

    % Pre-calculate global sums
    gene_sums_total = sum(X, 2); % g x 1 vector of row sums
    cell_sums_total = sum(X, 1); % 1 x c vector of column sums
    total_sum_X = sum(gene_sums_total); % Total sum of all elements in X

    if total_sum_X == 0
        warning('Total sum of X is zero. Returning zero residuals.');
        z_out = sparse(g, c); % Or use zeros(g,c)
        return;
    end

    % Initialize output matrix
    z_out = zeros(g, c, 'like', X);

    n_total_cells = c;
    sn = sqrt(n_total_cells); % Capping threshold

    fprintf('Processing %d cells in chunks of %d...\n', c, chunk_size_cells);
    for j = 1:chunk_size_cells:c
        idx_cells_start = j;
        idx_cells_end = min(j + chunk_size_cells - 1, c);
        
        current_cell_indices = idx_cells_start:idx_cells_end;
        % num_cells_in_chunk = length(current_cell_indices); % Not directly used later
        fprintf('  Processing cells %d to %d...\n', idx_cells_start, idx_cells_end);

        % Get the chunk of X. If X is sparse, X_chunk will be sparse.
        X_chunk = X(:, current_cell_indices);
        
        % Corresponding cell sums for this chunk
        cell_sums_chunk = cell_sums_total(current_cell_indices); % (1 x num_cells_in_chunk)

        % Calculate mu_chunk for the current cell chunk
        % mu_chunk will be (g x num_cells_in_chunk), typically dense.
        mu_chunk = (gene_sums_total * cell_sums_chunk) / total_sum_X;
        
        % Calculate NB variance-based sigma for the chunk
        % sigma_chunk will be (g x num_cells_in_chunk), dense.
        sigma_chunk = sqrt(mu_chunk + mu_chunk.^2 / theta + eps);
        
        % Calculate Pearson residuals for the chunk
        if issparse(X_chunk)
            z_chunk = (full(X_chunk) - mu_chunk) ./ sigma_chunk;
        else
            z_chunk = (X_chunk - mu_chunk) ./ sigma_chunk;
        end
        
        z_chunk(isnan(z_chunk)) = 0;
        
        % Cap residuals
        z_chunk(z_chunk >  sn) =  sn;
        z_chunk(z_chunk < -sn) = -sn;
        
        % Store the processed chunk
        z_out(:, current_cell_indices) = z_chunk;
    end
    fprintf('Finished processing all cell chunks.\n');
end