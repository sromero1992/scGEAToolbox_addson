function X_normalized = normalize_library_size_cell_chunks(X, scale_factor, chunk_size_cells)
    % Normalizes a count matrix X by library size, processing in cell-wise chunks.
    %
    % Formula: X_norm_ij = (X_ij / LibSize_j) * scale_factor
    %          where LibSize_j = sum_i X_ij (total counts for cell j)
    %
    % INPUT:
    % X                -----> Count matrix (genes x cells), can be sparse or full.
    % scale_factor     -----> Target sum for each cell's counts after normalization
    %                         (e.g., 1e4 for counts per 10,000). Default: 1e4.
    % chunk_size_cells -----> Number of cells to process in each chunk. Default: 1000.
    %                         (Your provided code had 10000, using 1000 as a common default).
    %
    % OUTPUT:
    % X_normalized     -----> Normalized matrix (genes x cells).
    %                         - Sparse (double precision) if input X is sparse.
    %                         - Dense if input X is dense (data type matches X if float, else double).

    if nargin < 2 || isempty(scale_factor)
        sf_param = 1e4; 
    else
        sf_param = scale_factor;
    end
    if nargin < 3 || isempty(chunk_size_cells)
        chunk_size_cells = 10000; % Default chunk size
    end

    [g, c] = size(X);
    is_input_sparse = issparse(X);

    % Determine calculation type and scale factor type
    if isa(X, 'single')
        calc_type = 'single';
        sf = single(sf_param);
    else % double, integer types
        calc_type = 'double'; % Integers will be converted to double for calculation
        sf = double(sf_param);
    end

    % Initialize output matrix
    if is_input_sparse
        % spalloc creates a double-precision sparse matrix.
        % Calculations will be performed in calc_type (single or double).
        % If calc_type is single, single results assigned to a double sparse matrix are upcast.
        X_normalized = spalloc(g, c, nnz(X)); 
        fprintf('Input is sparse. Outputting sparse (double precision) matrix.\n');
    else
        % For dense input, output type matches calc_type (which is derived from X's class)
        X_normalized = zeros(g, c, calc_type);
        fprintf('Input is dense. Outputting dense (%s) matrix.\n', calc_type);
    end

    fprintf('Normalizing %d cells in chunks of %d...\n', c, chunk_size_cells);
    % Suppress per-chunk printing if too verbose, or make it conditional
    verbose_chunks = (c > chunk_size_cells) && (c / chunk_size_cells > 5); % Example condition

    for j_start = 1:chunk_size_cells:c
        j_end = min(j_start + chunk_size_cells - 1, c);
        current_cell_indices = j_start:j_end;
        
        if verbose_chunks
            fprintf('  Processing cells %d to %d...\n', j_start, j_end);
        end
        
        % Get the chunk of X
        X_chunk_original = X(:, current_cell_indices);
        
        % Ensure X_chunk is of the desired calculation type (single or double)
        % This is important if X was integer, it needs to become float for division.
        if ~isa(X_chunk_original, calc_type)
            X_chunk = feval(calc_type, X_chunk_original);
        else
            X_chunk = X_chunk_original;
        end
        
        % Calculate library sizes for cells in the current chunk
        % sum operation on sparse 'single' matrix results in 'single' vector
        lib_sizes_chunk = sum(X_chunk, 1); 
        if ~isa(lib_sizes_chunk, calc_type) % Ensure sum is also calc_type (e.g. if X_chunk was int cast to calc_type)
            lib_sizes_chunk = feval(calc_type, lib_sizes_chunk);
        end
        
        % Prepare for division, handle zero library sizes
        lib_sizes_for_division_chunk = lib_sizes_chunk;
        zero_lib_size_cols_in_chunk = (lib_sizes_chunk == 0);
        % Use 'feval' to ensure '1' is of calc_type (e.g., single(1) or double(1))
        lib_sizes_for_division_chunk(zero_lib_size_cols_in_chunk) = feval(calc_type, 1);
        
        % Normalize the chunk.
        % If X_chunk is sparse and calc_type is 'single', result is single sparse.
        % If X_chunk is sparse and calc_type is 'double', result is double sparse.
        normalized_chunk_values = (X_chunk ./ lib_sizes_for_division_chunk) * sf;
        
        % Assign to X_normalized.
        % If X_normalized is double sparse (from spalloc) and normalized_chunk_values is single sparse,
        % MATLAB automatically upcasts the single values to double upon assignment.
        X_normalized(:, current_cell_indices) = normalized_chunk_values;
        
        % Ensure columns that originally had zero library size are truly zero in the output.
        % For sparse matrices, assigning 0 to elements effectively removes them if they were non-zero.
        if any(zero_lib_size_cols_in_chunk)
             cols_to_zero = current_cell_indices(zero_lib_size_cols_in_chunk);
             X_normalized(:, cols_to_zero) = 0; 
        end
    end
    
    fprintf('Finished processing all cell chunks.\n');
end