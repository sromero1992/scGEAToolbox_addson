function [T, Xsorted_completed, gsorted_completed, ...
    xyz1] = sc_splinefit_new(X, genelist, sortit, plotit, removenan)
    %SC_SPLINEFIT identify genes with a profile deviated from normal
    %
    % USAGE:
    % >> [X,genelist]=sc_readfile('example_data/GSM3044891_GeneExp.UMIs.10X1.txt');
    % >> [X]=sc_norm(X,'type','libsize');
    % >> [T]=sc_splinefit(X,genelist,true,true);
    
    % if nargin<2 || isempty(genelist)
    %     genelist=string(1:size(X,1))';
    %     genelist=strcat("gene_",genelist);
    % end

    if nargin < 5, removenan = false; end
    if nargin < 4, plotit = false; end
    if nargin < 3, sortit = true; end
    if nargin < 2 || isempty(genelist)
        genelist = string(1:size(X, 1)); 
    end
    
    % X = full(sce.X);
    % genelist = sce.g;
    % removenan = false;
    [lgu, dropr, lgcv, gsorted, Xsorted, ...
        removedgidx, removedT] = sc_genestat0(X, genelist, sortit, removenan); 

    if removenan || ~isempty(removedgidx)
        gsorted_completed = [gsorted; genelist(removedgidx)];
        Xsorted_completed = [Xsorted; X(removedgidx,:)];
        assert(isequal(size(Xsorted_completed),size(X)),'SC_SPLINEFIT')
        assert(length(gsorted_completed) == length(genelist),'SC_SPLINEFIT')
    else
        gsorted_completed = gsorted;
        Xsorted_completed = Xsorted;
    end
    
    xyz = [lgu, lgcv, dropr];

    s = cumsum([0; sqrt(   diff( lgu(:) ).^2  + ...
                           diff( lgcv(:)).^2 + ...
                           diff(dropr(:)).^2)]);

    % Old code 
    pp1 = splinefit0(s, xyz.', 15, 0.75);
    xyz1 = ppval(pp1, s)';

    % % New code from built in MATLAB 
    % % Adjust smoothing parameter for csaps to match the robustness parameter in splinefit0
    % % 0 is linear fit, 1 is nice fit to noisy data... 0.01 has bumps
    % % (cubic-spline behavior) small numbers make similar to previous code.
    % smoothing_param = 1e-6; % Adjust as needed 
    % % Define the desired number of grid points
    % num_points = 1000; %ceil(length(s)/10); % Adjust as needed
    % % Fit the cubic smoothing spline to each dimension
    % spl = csaps(s, xyz', smoothing_param);
    % % Generate new points for the fitted curves with more grid points
    % breaks_csaps = linspace(min(s), max(s), num_points); % Generate more points
    % % Evaluate the fitted spline at the specified breakpoints
    % xyz1 = fnval(spl, breaks_csaps)';
    % xyz1 = full(xyz1);

    % Remove possible spline points that should not be there
    % maxvals = full(max(xyz));
    % idx = xyz1(:,1) <= maxvals(1) & xyz1(:,2) <= maxvals(2) & xyz1(:,3) <= maxvals(3);
    % xyz1 = xyz1(idx,:);
    % minvals = full(min(xyz));
    % idx = xyz1(:,1) >= minvals(1) & xyz1(:,2) >= minvals(2) & xyz1(:,3) >= minvals(3);
    % xyz1 = xyz1(idx,:);

    %scatter3(xyz1(:,1), xyz1(:,2), xyz1(:,3))

    [nearidx, d] = dsearchn(xyz1, xyz);
    
    dx = d(d <= quantile(d, 0.9));
    
    distFit = fitdist([-dx; dx], 'Normal');
    pval = normcdf(d, 0, distFit.sigma, 'upper');
    [~, ~, ~, fdr] = pkg.fdr_bh(pval);
    
    if ~isempty(gsorted)
        genes = gsorted;
        T = table(genes, lgu, lgcv, dropr, d, pval, fdr, nearidx);
    else
        T = table(lgu, lgcv, dropr, d, pval, fdr, nearidx);
    end
    % 'variablenames',{'Genes','Log10_Mean','Dropout_Rate','Log10_CV','Deviation_3DFeature'});
    
    T.d(T.dropr > (1 - 0.05)) = 0; % ignore genes with dropout rate > 0.95
    T.d(T.dropr < (0.01)) = 0; % ignore genes with dropout rate < 0.01 (removes ribosomal and mitochondrial genes)
    
    % disp('NOTE: Genes with dropout rate > 0.95 are excluded.');
    
    if ~isempty(removedT) && istable(removedT)
        removedT.Properties.VariableNames = T.Properties.VariableNames;
        T = [T; removedT];
    end
    
    if sortit
        [T,idx] = sortrows(T, 'd', 'descend');
        gsorted_completed = gsorted_completed(idx);
        Xsorted_completed = Xsorted_completed(idx,:);
    end
    
    if length(gsorted_completed) ~= length(genelist)
        error('Output GENES are less than input GENES (some GENES are removed).');
    end
    
    if plotit
        figure;
        scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'filled', 'MarkerFaceAlpha', .1);
        hold on
        plot3(xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-', 'linewidth', 4);
        xlabel('Mean, log');
        ylabel('CV, log');
        zlabel('Dropout rate (% of zeros)');
    
        if ~isempty(gsorted)
            dt = datacursormode;
            dt.UpdateFcn = {@i_myupdatefcn3, gsorted, Xsorted};
        end
        hold off
    end

end