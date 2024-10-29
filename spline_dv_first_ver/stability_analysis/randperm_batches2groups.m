function sce = randperm_batches2groups(sce, seed)
    % This function grabs batches and splits it in 1 and 2 with ~50% of
    % random permutation to select the cell population
    % Author: Selim Romero
    % INPUT: 
    % sce ------> SCE object defined by scGEAToolbox 
    % seed -----> Random seed of preference (default 1234)
    % OUTPUT:
    % sce ------> Split of all available batches e.g. WT1, WT2, KO1 and KO2
    % USAGE:
    % do iter = 1:niter
    %     sce = randperm_batches2groups(sce, iter);
    % end
    
    batch = sce.c_batch_id;
    if  nargin < 2 || isempty(seed)
        seed = 1234;
    end
    %fprintf('Random seed is %d \n', seed);
    rng(seed);
    % Get unique batch names and number of batches
    batch_names = unique(batch);
    nbatches = length(batch_names);
    
    % Initialize cell arrays to hold indices and random permutations
    batch_indices = cell(nbatches, 1);
    shuffled_indices = cell(nbatches, 1);
    num_subsample = zeros(nbatches, 1);
    
    % Loop over each batch to perform subsampling and renaming
    for ibatch = 1:nbatches
    
        % Get the cell indices corresponding to the current batch
        batch_indices{ibatch} = find(batch == batch_names(ibatch));
        % Randomly shuffle the indices of the current batch
        indices = randperm(length(batch_indices{ibatch}));
        % Store the shuffled indices
        shuffled_indices{ibatch} = batch_indices{ibatch}(indices);
        % Calculate the number of cells to include in each new pseudo-batch (50% of total)
        num_subsample(ibatch) = round(length(indices) * 0.5);
        % Define the two ranges for pseudo-batches
        range1 = 1:num_subsample(ibatch);
        range2 = (num_subsample(ibatch) + 1):length(indices);
        % Rename the cells based on the new pseudo-batch assignment
        seed_str = strcat('seed-', string(seed));
        tag_name = strcat("1-",seed_str);
        sce.c_batch_id(shuffled_indices{ibatch}(range1)) = strcat(batch_names(ibatch), tag_name);
        tag_name = strcat("2-",seed_str);
        sce.c_batch_id(shuffled_indices{ibatch}(range2)) = strcat(batch_names(ibatch), tag_name);
    end

end