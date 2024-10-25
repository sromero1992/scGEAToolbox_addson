function Tdv = sc_splinedv(X1, X2, g1, g2, fname, plotit)
    %{
    sc_splinedv computes the intersection of two experiments/batches and
    computes a cubic spline for each set to compute variability.
    INPUT
        X ---------> Count matrix of set1 (Normalized)
        X2 --------> Count matrix of set2 (Normalized)
        genelist --> Gene list of set1
        genelist2 -> Gene list of set2
        fname -----> file name of tables
        plotit ----> Boolean to plot the results
    OUTPUT
        Tdv_final -> Differentially variable table containing gene info
    %}
    
    if nargin < 6; plotit = false; end
    if nargin < 5; fname = []; end

    sortit = true;
   
    % Find intersected genes and sort
    [g3, irows, jrows] = intersect(g1, g2, 'stable');
    X1 = X1(irows, :); 
    X2 = X2(jrows, :); 

    % Normalize 
    X1 = sc_norm(X1, 'type', 'libsize');
    X2 = sc_norm(X2, 'type', 'libsize');

    % Compute cubic spline polynomial for set 1
    [T1, X1, g1, xyzp1] = sc_splinefit(X1, g3, sortit, false);
    [T1, idx] = sortrows(T1, 'genes', 'ascend');
    xyz1 = [T1.lgu, T1.lgcv, T1.dropr];
    X1 = X1(idx, :);
    g1 = g1(idx);  % Use the sorted version from above

    % Compute cubic spline polynomial for set 2
    [T2, X2, g2, xyzp2] = sc_splinefit(X2, g3, sortit, false);
    [T2, idx2] = sortrows(T2, 'genes', 'ascend');
    xyz2 = [T2.lgu, T2.lgcv, T2.dropr];
    X2 = X2(idx2, :);
    g2 = g2(idx2);  % Also sorted

    % Ensure gene lists are the same after sorting
    assert(isequal(g1, g2));

    % Spline to gene vector positions
    v1 = xyz1 - xyzp1(T1.nearidx, :);
    v2 = xyz2 - xyzp2(T2.nearidx, :);

    % Distance within two relative variable points
    dist_diff = vecnorm(v1 - v2, 2, 2);
    dist_sign = sign(vecnorm(v1, 2, 2) - vecnorm(v2, 2, 2));

    % Standardize distances
    ddz = zscore(dist_diff);
    % Assume Gaussian normal distribution for errors for p-values
    pval = 1 - normcdf(ddz);

    % Rename variables in T1 and T2 to avoid duplicates
    T1.Properties.VariableNames = strcat(T1.Properties.VariableNames, '1');
    T2.Properties.VariableNames = strcat(T2.Properties.VariableNames, '2');
    
    % Merge data into final table
    T1.Properties.VariableNames{1} = 'gene';  % Assuming gene list is in first column
    Tdv = [T1, T2, table(dist_diff), table(dist_sign), table(pval)];

    % Do we need to do this? maybe not...
    % Zero out certain distances based on index conditions
    idxx = T1.(8) == 1 | T2.(8) == 1 | T1.(8) == max(T1.(8)) | T2.(8) == max(T2.(8));
    Tdv.dist_diff(idxx) = 0;

    % Sort by distance difference
    Tdv = sortrows(Tdv, "dist_diff", "descend");

    % Save the differential variability table
    if ~isempty(fname)
        filenameT = strcat(fname, '_diff_var.csv');
        writetable(Tdv, filenameT, 'Delimiter', ',');
    end
    % Plot the splines and gene statistics if plotit is true
    if plotit
        figure;
        
        % Plot for Set 1 (cyan edges with transparency, blue spline)
        scatter3(xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), 50, 'c', 'MarkerEdgeAlpha', 0.5, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'c'); % Cyan edge with transparency, no fill
        hold on;
        plot3(xyzp1(:, 1), xyzp1(:, 2), xyzp1(:, 3), 'b-', 'LineWidth', 2); % Blue solid line for spline
        
        % Plot for Set 2 (bright green edges with transparency, dark green spline)
        scatter3(xyz2(:, 1), xyz2(:, 2), xyz2(:, 3), 50, 'g', 'MarkerEdgeAlpha', 0.5, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'g'); % Green edge with transparency, no fill
        plot3(xyzp2(:, 1), xyzp2(:, 2), xyzp2(:, 3), 'r-', 'LineWidth', 2); % Dark green solid line for spline
        
        % Title and labels
        title('Comparison of Set 1 and Set 2: Gene Statistics vs. Spline');
        xlabel('Log1p Mean');
        ylabel('Log1p CV');
        zlabel('Dropout Rate');
        
        % Legend
        legend({'Set 1: Gene Stats', 'Set 1: Spline Fit', 'Set 2: Gene Stats', 'Set 2: Spline Fit'});
        
        hold off;
    end

end

