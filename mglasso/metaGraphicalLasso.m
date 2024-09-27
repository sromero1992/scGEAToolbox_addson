function [Sigma_inv, eps_k, D] = metaGraphicalLasso(X, M, rho, learning_rate, num_epochs, batch_size)
    % Meta Graphical Lasso (MG-Lasso) implementation in MATLAB.
    % Inputs:
    %   X: Cell array of datasets {X1, X2, ..., Xk} where Xk is N_k x V (sample x variables)
    %   M: Number of latent factors
    %   rho: Weight of L1 regularization term
    %   learning_rate: Learning rate for optimization
    %   num_epochs: Number of epochs for training
    %   batch_size: Batch size for training
    % Outputs:
    %   Sigma_inv: Precision matrices of latent factors
    %   eps_k: Reciprocal of epsilon values
    %   D: Projection matrix
    
    K = length(X); % Number of datasets
    N = cellfun(@(x) size(x,1), X); % Number of samples in each dataset
    V = size(X{1}, 2); % Number of observed variables
    
    % Initialize precision matrices and projection matrix
    Sigma_inv = cell(K,1);
    for k = 1:K
        Sigma_inv{k} = eye(M); % Initialize with identity matrices
    end
    
    % Initialize projection matrix D using SVD
    Xall = cell2mat(cellfun(@(x) x./size(x,1), X, 'UniformOutput', false));
    [U, ~, ~] = svd(Xall, 'econ');
    D = U(:, 1:M);
    
    % Gradient Descent optimization loop
    for epoch = 1:num_epochs
        % Shuffle dataset indices
        indices = randperm(K);
        
        for i = 1:batch_size:K
            batch_indices = indices(i:min(i+batch_size-1, K));
            D_grad = zeros(size(D));
            
            for bi = batch_indices
                Xk = X{bi};
                Nk = size(Xk, 1);
                
                % Update Sigma_inv and eps_k for this batch
                [Sigma_inv{bi}, eps_k(bi)] = optimizeLmd(D, Xk, Sigma_inv{bi}, Nk, V, M, rho);
                
                % Compute gradient for D
                DX = D' * Xk';
                DXXD = DX * DX';
                omega = Sigma_inv{bi} - eps_k(bi) * eye(M);
                grad = D * (DXXD * omega + omega' * DXXD') - 2 * (Xk' * (omega * DX));
                D_grad = D_grad - grad / Nk;
            end
            
            % Update D using the calculated gradient
            D = D - learning_rate * D_grad;
            
            % Project D onto Stiefel manifold using SVD
            [P, ~, Q] = svd(D, 'econ');
            D = P * Q';
        end
        
        % Display loss for each epoch
        total_loss = 0;
        for k = 1:K
            total_loss = total_loss + calcLoss(D, Sigma_inv{k}, eps_k(k), X{k}, N(k), V, M);
        end
        disp(['Epoch ', num2str(epoch), ' Loss: ', num2str(total_loss)]);
    end
end
