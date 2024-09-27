function [Sigma_inv, eps_k] = optimizeLmd(D, Xk, Sigma_inv, Nk, V, M, rho)
    % Graphical Lasso optimization for each dataset
    
    % Step 1: Calculate empirical covariance matrix for X_k
    emp_cov = (D' * Xk') * (Xk * D) / Nk;
    
    % Step 2: Apply graphical Lasso
    try
        % Use your custom graphicalLasso function
        [Sgm, Sigma_inv] = graphicalLasso(emp_cov, rho);  
    catch
        warning('Graphical Lasso failed. Using previous Sigma_inv.');
    end
    
    % Step 3: Calculate epsilon
    eps_k = calcEpsilon(D, Xk, Nk, V, M);
end
