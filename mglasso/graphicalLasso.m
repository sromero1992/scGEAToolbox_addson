function [Theta, S] = graphicalLasso(Sigma, lambda, max_iter, tol)
    % graphicalLasso: Solves the graphical lasso problem using coordinate descent
    % Args:
    %     Sigma: Empirical covariance matrix (V x V)
    %     lambda: Regularization parameter for L1 penalty
    %     max_iter: Maximum number of iterations (default: 100)
    %     tol: Convergence tolerance (default: 1e-4)
    % Returns:
    %     Theta: Estimated precision (inverse covariance) matrix (V x V)
    %     S: Estimated covariance matrix (V x V)

    if nargin < 3, max_iter = 100; end
    if nargin < 4, tol = 1e-4; end

    V = size(Sigma, 1);
    Theta = inv(Sigma + lambda * eye(V));  % Initial precision matrix
    W = Theta;  % Working precision matrix

    for iter = 1:max_iter
        Theta_old = Theta;

        for j = 1:V
            not_j = setdiff(1:V, j);

            % Extract submatrices and solve for j-th row/column using lasso
            S11 = Sigma(not_j, not_j);
            s12 = Sigma(not_j, j);

            % Solve Lasso using coordinate descent for the j-th row
            beta = lassoCoordDescent(S11, s12, lambda);

            % Update precision and working matrices
            W(not_j, j) = beta;
            W(j, not_j) = beta';
            W(j, j) = 1 / (Sigma(j, j) - s12' * beta);
        end

        % Symmetrize and check convergence
        Theta = (W + W') / 2;
        if norm(Theta - Theta_old, 'fro') < tol
            fprintf('Converged after %d iterations\n', iter);
            break;
        end
    end
    S = inv(Theta);  % Covariance from precision
end

function beta = lassoCoordDescent(S, s, lambda)
    % lassoCoordDescent: Solves Lasso using coordinate descent
    % Args:
    %     S: Covariance matrix for non-j columns (V-1 x V-1)
    %     s: Covariance vector for j-th column (V-1 x 1)
    %     lambda: Regularization parameter
    % Returns:
    %     beta: Coefficients of Lasso problem (V-1 x 1)

    V_minus_1 = length(s);
    beta = zeros(V_minus_1, 1);
    max_iter = 100;
    tol = 1e-4;

    for iter = 1:max_iter
        beta_old = beta;

        for j = 1:V_minus_1
            % Residual calculation for j-th coordinate
            r = s(j) - S(j, :) * beta + S(j, j) * beta(j);

            % Soft-thresholding for Lasso
            if abs(r) > lambda
                beta(j) = (r - sign(r) * lambda) / S(j, j);
            else
                beta(j) = 0;
            end
        end

        % Convergence check
        if norm(beta - beta_old, 'fro') < tol
            break;
        end
    end
end

