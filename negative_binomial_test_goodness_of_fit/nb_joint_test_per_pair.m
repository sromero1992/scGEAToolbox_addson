function nb_joint_test_per_pair(sce, igene, jgene, thresholds, plotit, copula_type)
    % nb_joint_test_per_pair tests the joint distribution of two genes using a copula model.
    %
    % INPUT:
    %   - sce: SingleCellExperiment object with gene expression data.
    %   - igene: Name of the first gene.
    %   - jgene: Name of the second gene.
    %   - thresholds: An array of thresholds for filtering gene counts.
    %   - plotit: Boolean to determine whether to plot joint distribution histograms.
    %   - copula_type: Type of copula to fit (e.g., 'Gaussian', 't', 'Clayton', 'Frank').
    %
    % OUTPUT:
    %   None. Prints AIC and p-values for joint distributions.

    if nargin < 5 || isempty(plotit); plotit = false; end
    if nargin < 6 || isempty(copula_type); copula_type = 'Gaussian'; end

    idx1 = find(strcmp(igene, sce.g));
    idx2 = find(strcmp(jgene, sce.g));
    raw_gene_count1 = sce.X(idx1, :)';
    raw_gene_count2 = sce.X(idx2, :)';

    nthrs = length(thresholds);
    aic_values_nb1 = -1*ones(nthrs, 1);
    aic_values_nb2 = -1*ones(nthrs, 1);
    aic_values_joint = -1*ones(nthrs, 1);
    p_values_joint = -1*ones(nthrs, 1);

    for i = 1:length(thresholds)
        thrs = thresholds(i);
        idx = or(raw_gene_count1 >= thrs, raw_gene_count2 >= thrs);
        if sum(idx) <= 5; continue; end
        filt_gene_count1 = raw_gene_count1(idx);
        filt_gene_count2 = raw_gene_count2(idx);

        % Fit NB distributions for both genes
        pd_nb1 = fitdist(filt_gene_count1, 'NegativeBinomial');
        pd_nb2 = fitdist(filt_gene_count2, 'NegativeBinomial');
        
        % AIC for NB distributions
        nbpdf_full1 = nbinpdf(filt_gene_count1, pd_nb1.R, pd_nb1.p);
        logL_full1 = sum(log(nbpdf_full1 + eps));
        aic_values_nb1(i) = -2 * logL_full1 + 2 * pd_nb1.NumParameters;
        
        nbpdf_full2 = nbinpdf(filt_gene_count2, pd_nb2.R, pd_nb2.p);
        logL_full2 = sum(log(nbpdf_full2 + eps));
        aic_values_nb2(i) = -2 * logL_full2 + 2 * pd_nb2.NumParameters;
        
        % Plot joint distributions
        if plotit
            figure;
            histogram2(filt_gene_count1, filt_gene_count2, 'DisplayStyle', 'tile');
            xlabel([igene ' Counts']);
            ylabel([jgene ' Counts']);
            title(['Joint Distribution Histogram for ', igene, '-', jgene, ...
                   ' with threshold ', num2str(thrs)]);
        end

        % Fit the copula of the specified type
        u1 = nbincdf(filt_gene_count1, pd_nb1.R, pd_nb1.p);
        u2 = nbincdf(filt_gene_count2, pd_nb2.R, pd_nb2.p);
        
        if strcmp(copula_type, 't')
            % Fit t-Copula
            [copula_params, d_f] = copulafit('t', [u1, u2]);
            joint_pdf = @(x1, x2) copulapdf('t', [nbincdf(x1, pd_nb1.R, pd_nb1.p), nbincdf(x2, pd_nb2.R, pd_nb2.p)], copula_params, d_f) ...
                                  .* nbinpdf(x1, pd_nb1.R, pd_nb1.p) .* nbinpdf(x2, pd_nb2.R, pd_nb2.p);
        elseif any(strcmp(copula_type, {'Gaussian', 'Clayton', 'Frank'}))
            % Fit Gaussian, Clayton, or Frank Copula
            copula_params = copulafit(copula_type, [u1, u2]);
            joint_pdf = @(x1, x2) copulapdf(copula_type, [nbincdf(x1, pd_nb1.R, pd_nb1.p), nbincdf(x2, pd_nb2.R, pd_nb2.p)], copula_params) ...
                                  .* nbinpdf(x1, pd_nb1.R, pd_nb1.p) .* nbinpdf(x2, pd_nb2.R, pd_nb2.p);
        else
            error('Unsupported copula type. Choose from ''Gaussian'', ''t'', ''Clayton'', or ''Frank''.');
        end

        % Compute AIC for the joint model based on the chosen copula
        joint_prob = joint_pdf(filt_gene_count1, filt_gene_count2);
        logL_joint = sum(log(joint_prob + eps));
        aic_values_joint(i) = -2 * logL_joint + 4;  % 2 parameters for each NB dist + 2 for copula

        % Compute p-value using a chi-squared test statistic
        try
            % Generate empirical and model CDFs
            empirical_cdf = sum(u1 <= nbincdf(filt_gene_count1, pd_nb1.R, pd_nb1.p) & u2 <= nbincdf(filt_gene_count2, pd_nb2.R, pd_nb2.p)) / length(u1);
            model_cdf = mean(joint_pdf(filt_gene_count1, filt_gene_count2));
            test_stat = sum((empirical_cdf - model_cdf).^2);  % Basic test statistic
            p_values_joint(i) = 1 - chi2cdf(test_stat, 2);    % Assuming 2 degrees of freedom
        catch ME
            fprintf('Error with chi2gof: %s\n', ME.message);
            p_values_joint(i) = NaN; % Set to NaN if the test fails
        end
    end

    % Display results
    fprintf('Threshold\tAIC (Joint)\tp-value (Joint)\n');
    for i = 1:length(thresholds)
        fprintf('%d\t\t%.4f\t\t%.4f\n', thresholds(i), aic_values_joint(i), p_values_joint(i));
    end
end

