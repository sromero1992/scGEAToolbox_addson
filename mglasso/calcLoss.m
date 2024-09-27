function loss = calcLoss(D, Sigma_inv, eps_k, Xk, Nk, V, M)
    % Calculate loss (likelihood term)
    loss = trace((D' * Xk') * (Xk * D) * Sigma_inv) - trace(Sigma_inv) * eps_k;
    loss = loss + (sum(sum(Xk.^2)) * eps_k) / Nk;
    loss = loss - log(det(Sigma_inv)) - (V - M) * log(eps_k);
end
