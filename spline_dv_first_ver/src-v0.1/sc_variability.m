function [Tdv_final, Tig, influgenes] = sc_variability( X, X2, genelist, genelist2, fname, plotit )
    %{ 
    sc_variability computes the intersection of two experimets/batches and
    computes a cubic spline for each set to compute the variability.
    INPUT
        X ---------> Count matrix of set1 (Normalized)
        X2 --------> Count matrix of set2 (Normalized)
        genelist --> Gene list of set1
        genelist2 -> Gene list of set2
        fname -----> file name of tables
    OUTPUT
        Tdv_final -> Differentially variable table containing gene info
        Tig -------> Influenciable gene table containing polyfit info
        influgenes-> Influenciable gene list from intersected datasets
    NOTE
        1.- The intersection needs to be computed to use
            Tdv_final.dvgenes_idx as follows:
            [gene,~,~] = intersect(genelist, genelist2, 'stable'))
            gene(Tdv_final.dvgenes_idx)

        2.- Gene lists and tablers are sorted by genes to have same 
            indices in same places.

        3.- List is reduced by the 10% largest absolute differences 
            in polyfit-gene distance (dvgenes_idxf).
    USAGE
        Selected genes expressed in at least %7.5
        [X,genelist]=sc_selectg(X,genelist,1,0.075);
        [X]=sc_norm(X,'type','libsize');
        [X2,genelist2]=sc_selectg(X2,genelist2,1,0.075);
        [X2]=sc_norm(X2,'type','libsize');

        [Tdv, Tig, influgenes] = sc_variability( X, X2, genelist, genelist2, false);
    %}

    if nargin < 5, plotit = false; end
    % This needs to be done, otherwise polyfit gets junk
    sortit = true;

    [genelist3, irows, jrows] = intersect(genelist, genelist2, 'stable');

    % T = table(lgu, lgcv, dropr, d, pval, fdr);
    % Compute cubic spline polynomial of noisy data
    [T,~,genelist,xyzp] = sc_splinefit( X(irows, :), genelist3, sortit);
    xyz = [T.lgu, T.lgcv, T.dropr];
    clear T;
    %[xyz, xyzp] = minmax_scale(xyz, xyzp);
    % Remove dropout rates above and below such threasholds remove MT and R
    %idx = (xyz(:,3) > 0.01) & (xyz(:,3)<0.99);
    %xyz = xyz( idx,:);
    %genelist = genelist(idx);


    % Compute cubic spline polynomial of noisy data
    [T2,~,genelist2,xyzp2]= sc_splinefit( X2(jrows, :), genelist3, sortit);
    xyz2 = [T2.lgu, T2.lgcv, T2.dropr];
    clear T2;
    %[xyz2, xyzp2] = minmax_scale(xyz2, xyzp2);
    % Remove dropout rates above and below such threasholds remove MT and R
    %idx = (xyz2(:,3) > 0.01) & (xyz2(:,3) < 0.99) ;
    %xyz2 = xyz2( idx,:);
    %genelist2 = genelist2(idx);

    if plotit 
        figure;
        scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), 'filled', 'MarkerFaceAlpha', .1);
        hold on
        plot3(xyzp(:, 1), xyzp(:, 2), xyzp(:, 3), '-', 'linewidth', 4);
        hold on
        scatter3(xyz2(:, 1), xyz2(:, 2), xyz2(:, 3), 'filled', 'MarkerFaceAlpha', .1);
        hold on
        plot3(xyzp2(:, 1), xyzp2(:, 2), xyzp2(:, 3), '-', 'linewidth', 4);
        xlabel('Mean, log');
        ylabel('CV, log');
        zlabel('Dropout rate (% of zeros)');
    end

    % Position vector referenced at polynomial
    % dsearchn(P,QP) utilizes P data points and QP query-test points 
    % xyz(idx,1:3) will sort the points to match xyzp order and size
    %[idx,dist] = dsearchn(xyz, xyzp); 
    [kpts,dist] = dsearchn(xyzp, xyz);
    % Data point-spline vector position
    xyzp = xyzp(kpts, 1:3);
    xyz_dist = xyz - xyzp;

    Tdv = table(genelist, ...
                xyz_dist(:,1), xyz_dist(:,2) , xyz_dist(:,3), ...
                xyzp(:,1), xyzp(:,2), xyzp(:,3),...
                dist);

    clear xyzp;

    Tdv.Properties.VariableNames = {'genes', 'lgmu_dist','lgcv_dist',...
                                    'dropr_dist','s1','s2','s3',...
                                    'dist'};

    [Tdv, idx_tdv] = sortrows(Tdv,'genes','descend');
    xyz = xyz(idx_tdv,:);

    % Position vector referenced at polynomial
    [kpts2,dist2] = dsearchn(xyzp2, xyz2);
    % Data point-spline vector position
    xyzp2 = xyzp2(kpts2,1:3);
    xyz_dist2 = xyz2 - xyzp2;
    
    Tdv2 = table(genelist2, ...
                xyz_dist2(:,1), xyz_dist2(:,2) , xyz_dist2(:,3), ...
                xyzp2(:,1), xyzp2(:,2), xyzp2(:,3),...
                dist2);

    clear xyzp2;

    Tdv2.Properties.VariableNames = {'genes', 'lgmu_dist','lgcv_dist',...
                                    'dropr_dist','s1','s2','s3',...
                                    'dist'};

    [Tdv2, idx_tdv2] = sortrows(Tdv2,'genes','descend');
    xyz2 = xyz2(idx_tdv2,:);

    % (D1 - S1) - (D2 - S2)  or distance difference in vector form
    % (D1 - D2) + (S2 - S1)
    ddv = table2array( Tdv(:,2:4) - Tdv2(:,2:4) );
    % % Computing ( S1 - S2 )
    ddvs = table2array( Tdv(:,5:7) - Tdv2(:,5:7) );
    % % Computing ( D1 - D2 )
    ddvp = xyz(:,1:3) - xyz2(:,1:3);
    % % Distance within splines + 1
    ddscale = vecnorm(ddvs,2,2) + 1;
    % % Distance within splines + distance within points
    % ddscale = vecnorm(ddvs,2,2) + vecnorm(ddvp,2,2);
    clear ddvs ddvp;

    dd0 = vecnorm(ddv,2,2);
    dd =  dd0./ddscale;

    % Standarize difference distance
    ddz = zscore(dd0);
    ddz2 = zscore(dd);


    % Cumulative distribution function (cdf) according to
    % a standard normal distribution in d from dx sigma
    pval = 2*normcdf(ddz,'upper');
    pval2 = 2*normcdf(ddz2,'upper');


    % Table containing distance metrics of differentially variable genes
    Tdv_final = table(Tdv.genes, ...
                      Tdv.lgmu_dist, Tdv.lgcv_dist,...
                      Tdv.dropr_dist, xyz(:,3), Tdv.dist,...         
                      Tdv2.lgmu_dist, Tdv2.lgcv_dist, ...
                      Tdv2.dropr_dist, xyz2(:,3), Tdv2.dist, ...
                      ddv(:,1), ddv(:,2), ddv(:,3),...
                      ddscale, dd0,...
                      dd, ddz, ddz2, pval, pval2);

    clear pval ddz dd dd0 ddv Tdv Tdv2;

    Tdv_final.Properties.VariableNames = {'genes', ...
                                        'lgmu_dist1', 'lgcv_dist1',...
                                        'dropr_dist1','drop1','dist_set1',...
                                        'lgmu_dist2', 'lgcv_dist2', ...
                                        'dropr_dist2','drop2','dist_set2', ...
                                        'lgmu_dd','lgcv_dd','dropr_dd',...
                                        'dd_scale', 'dd_raw',...
                                        'dd_norm', 'dd_std_raw', 'dd_std_scale'...
                                        'pval_raw', 'pval_scale'};

    % Sort together by pval acending and abs_dd decensindg
    idx = and(Tdv_final.drop1 > 0.01, Tdv_final.drop1 < 0.99); 
    Tdv_final = Tdv_final(idx,:);
    idx = and(Tdv_final.drop2 > 0.01, Tdv_final.drop2 < 0.99); 
    Tdv_final = Tdv_final(idx,:);

    [Tdv_final, ~] = sortrows(Tdv_final, 20, {'ascend'});
   
    % Writting differential variability table 
    filenameT = '_diff_var.csv';
    filenameT = strcat(fname,filenameT);
    writetable(Tdv_final,filenameT,'Delimiter',',');
    type(filenameT);

    % %%%%%%%%%%%%%%% Termporarily shuted down %%%%%%%%%%%%%%%
    active = false;
    if active 
        % Is polynomial difference saying something else? is this a scam?
        % kpolyg_idx contains the order for accessing xyzp (mapping)
        x = Tdv.p1 - Tdv2.p2;
        y = Tdv.p2 - Tdv2.p2;
        z = Tdv.p3 - Tdv2.p3;
        distp = x.^2 + y.^2 + z.^2;
        distp = distp.^0.5;
        
        % Standarize difference distance
        ddz = zscore(distp);

        % Cumulative distribution function (cdf) according to
        % a standard normal distribution in d from dx sigma
        pval = 2*normcdf(ddz,'upper');

        % Table containing influenciable genes?
        Tig = table(Tdv.genes, x, y, z, distp, ddz, pval);

        clear x y z disp ddz pval;

        Tig.Properties.VariableNames = {'genes','dist_x',...
                                        'dist_y','dist_z',...
                                        'polyfit_dist', ...
                                        'dist_sdt','pval'};

        [Tig, ~] = sortrows(Tig, [7 5], {'ascend' 'descend'});

        % Writting differential variability polynomial fit table 
        filenameT = '_diff_var_polyfit.csv';
        filenameT = strcat(fname,filenameT);
        writetable(Tig,filenameT,'Delimiter',',');
        type(filenameT);

        % Writting intersection of polyfit and differential variability
        [influgenes, ~, ~] = intersect(Tig.genes(1:1000), ...
                                       Tdv_final.genes(1:1000), 'stable');
        filenameT = '_influgenes.csv';
        filenameT = strcat(fname,filenameT);
        writematrix(influgenes, filenameT);
        type(filenameT);
    else 
        Tig = 0;
        influgenes = 0;
    end
end