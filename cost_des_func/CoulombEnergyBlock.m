function Ecoul_block = CoulombEnergyBlock(X_block1, X_block2, MI_block, eps, isDiagonal)
    % Compute the Coulomb energy for a block of genes
    ng1 = size(X_block1, 2);  % Number of genes in block 1 (columns)
    ng2 = size(X_block2, 2);  % Number of genes in block 2 (columns)
    nc = size(X_block1, 1);   % Number of cells (rows)

    % Initialize the block for Coulomb energy
    Ecoul_block = zeros(ng1, ng2);
    
    % Parallel loop over block 1 (genes in block 1)
    if isDiagonal
        parfor ig = 1:ng1
            e1 = dot(X_block1(:, ig), X_block1(:, ig)) / nc^2;  % Self-energy for ig

            % Temporary variable to avoid race conditions
            Etmp = zeros(ng2, 1);
            MItmp = MI_block(:, ig);

            % Diagonal and upper triangular part
            for jg = ig:ng2
                diffv = X_block1(:, ig) - X_block2(:, jg);

                if jg == ig
                    % Diagonal element
                    Etmp(jg) = -1.0 * (e1^2 + MI_block(ig, jg)) / (sqrt(diffv' * diffv) + eps);
                else
                    % Upper triangular part

                    e2 = dot(X_block1(:, ig), X_block2(:, jg)) / nc^2;
                    %Etmp(jg) = 1.0 * (e1 * e2 - MI_block(ig, jg)) / (sqrt(diffv' * diffv) + eps);
                    Etmp(jg) = 1.0 * (e1 * e2 - MItmp(jg)) / (sqrt(diffv' * diffv) + eps);

                end
            end

            % Assign the computed values for ig to the corresponding row
            Ecoul_block(ig, :) = Etmp;
        end
    else
        parfor ig = 1:ng1
            %e1 = dot(X_block1(:, ig), X_block1(:, ig)) / nc^2;  % Self-energy for ig
            % Temporary variable to avoid race conditions
            Etmp = zeros(ng2, 1);

            % Off-diagonal block (X_block1 and X_block2 are different)
            for jg = 1:ng2
                diffv = X_block1(:, ig) - X_block2(:, jg);
                e2 = dot(X_block1(:, ig), X_block2(:, jg)) / nc^2;
                Etmp(jg) = 1.0 * (e2^2 - MI_block(ig, jg)) / (sqrt(diffv' * diffv) + eps);
            end

            % Assign the computed values for ig to the corresponding row
            Ecoul_block(ig, :) = Etmp(:);
        end
    end
end
