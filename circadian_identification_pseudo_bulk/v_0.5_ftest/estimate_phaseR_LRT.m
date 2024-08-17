function [acrophase, amp, period, mesor, p_value] = estimate_phaseR_LRT(Xg_zts, time_step, period12)
    % This function estimates the phase, amplitude, period, and mesor from gene expression data
    % using a sine function fit. It also computes the p-value for the significance of the sine fit
    % using a Likelihood Ratio Test (LRT) comparing it to a null model (mean model).
    %
    % Inputs:
    %   - Xg_zts: A cell array where each cell contains the expression data of all cells at each time point.
    %   - time_step: The time difference between successive time points.
    %   - period12: A boolean flag indicating whether to fit a 12-hour (true) or 24-hour (false) period sine function.
    %
    % Outputs:
    %   - acrophase: The estimated acrophase (time of peak expression).
    %   - amp: The estimated amplitude of the sine function.
    %   - period: The period used for fitting (12 or 24 hours).
    %   - mesor: The estimated mesor (mean level of expression).
    %   - p_value: The p-value from the Likelihood Ratio Test (LRT) indicating the significance of the sine model fit.

    % Calculate the number of time points (nzts)
    nzts = size(Xg_zts, 2);

    % Initialize arrays for gene expression values (R) and corresponding time points (time_grid)
    icells = cellfun(@length, Xg_zts, 'UniformOutput', true); % Number of cells at each time point
    num_cells = sum(icells); % Total number of cells across all time points

    R = zeros(num_cells, 1); % Flattened array for expression data
    time_grid = zeros(num_cells, 1); % Corresponding time points for each cell

    % Reshape the data from cell array to flat arrays for fitting
    ic = 0;
    for it = 1:nzts
        % Assign expression data
        R(ic + 1:ic + icells(it)) = Xg_zts{it}(:);
        % Assign corresponding time points
        time_grid(ic + 1:ic + icells(it)) = (it - 1) * time_step;
        ic = ic + icells(it);
    end
    
    % Null model (mean model)
    mean_model = mean(R); % Mean of all expression data
    SSR_null = sum((R - mean_model).^2); % Sum of squared residuals for the null model
    logL_null = -0.5 * length(R) * log(SSR_null / length(R)); % Log-likelihood of the null model

    % Define the sine model based on the selected period (12 or 24 hours)
    if period12
        ft = fittype('amp * cos(2*pi*(t - acro)/12) + mesor', ...
                     'coefficients', {'acro', 'amp', 'mesor'}, ...
                     'independent', {'t'});
        period = 12;
    else
        ft = fittype('amp * cos(2*pi*(t - acro)/24) + mesor', ...
                     'coefficients', {'acro', 'amp', 'mesor'}, ...
                     'independent', {'t'});
        period = 24;
    end
    
    % Fit options for the sine model
    options = fitoptions('Method', 'NonlinearLeastSquares', ...
                         'Algorithm', 'Trust-Region', ...
                         'Lower', [-24, -Inf, -Inf], ...
                         'Upper', [24, Inf, Inf], ...
                         'StartPoint', [0, mean_model * 2, mean_model]);
    
    % Fit the sine model to the data
    [fmdl, gof] = fit(time_grid, R, ft, options);
    SSR_sine = gof.sse; % Sum of squared errors for the sine model
    logL_sine = -0.5 * length(R) * log(SSR_sine / length(R)); % Log-likelihood of the sine model

    % Calculate the Likelihood Ratio Test (LRT) statistic
    LRT_stat = -2 * (logL_null - logL_sine);

    % Degrees of freedom difference (sine model has 3 parameters, null model has 1)
    df = 3 - 1;

    % Calculate the p-value using the chi-square distribution
    p_value = 1 - chi2cdf(LRT_stat, df);
    
    % Outputs: the estimated parameters from the sine fit
    acrophase = fmdl.acro;
    amp = fmdl.amp;
    mesor = fmdl.mesor;
end
