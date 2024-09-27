function eps_k = calcEpsilon(D, Xk, Nk, V, M)
    % Compute epsilon (eps_k) for each dataset
    eps_k = sum(sum(Xk.^2)) - sum(sum((D' * Xk').^2));
    eps_k = (V - M) / (Nk * eps_k);
end
