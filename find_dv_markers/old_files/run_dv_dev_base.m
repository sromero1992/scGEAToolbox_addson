main_path = "HFD_LFD/";
sample_id = "HFD_LFD";
path = strcat(main_path,sample_id,".mat");
data  = load(path);
sce1 = data.sce;
clear data;

CellTypeList1 = unique(sce1.c_cell_type_tx);
celltype = "Adipocytes";
celltypef = replace(celltype," ","_");
label1 = strcat(celltypef,"-", sample_id);
fname1 = strcat(main_path,label1);

% Take out 5% mitochondrial cells and corresponding batch id
cellidx = and(sce1.c_batch_id == "HFD", sce1.c_cell_type_tx == celltype);
[X, g ] = sc_selectg(sce1.X(:,cellidx), sce1.g, 1, 0.05);
cellidx2 = and(sce1.c_batch_id == "LFD", sce1.c_cell_type_tx == celltype);
[X2,g2] = sc_selectg(sce1.X(:,cellidx2), sce1.g, 1, 0.05);


% Work only with common genes
[gl, irows, jrows] = intersect(g, g2, 'stable');
X = X(irows,:);
X2 = X2(jrows,:);

% Assing a boolean to normalize
[X] = sc_norm(X,'type','libsize');
[X2]= sc_norm(X2,'type','libsize');

genelist = gl;
genelist2 = gl;
sortit = true;

    [genelist3, irows, jrows] = intersect(genelist, genelist2, 'stable');

    % T = table(lgu, lgcv, dropr, d, pval, fdr);
    % Compute cubic spline polynomial of noisy data
    [T,X,genelist,xyzp] = sc_splinefit( X(irows, :), genelist3, sortit);
    xyz = [T.lgu, T.lgcv, T.dropr];
    [xyz, xyzp] = minmax_scale(xyz, xyzp);

    %clear T;
    % Remove dropout rates above and below such threasholds remove MT and R
    %idx = (xyz(:,3) > 0.01) & (xyz(:,3)<0.99);
    %xyz = xyz( idx,:);
    %genelist = genelist(idx);


    % Compute cubic spline polynomial of noisy data
    [T2,X2,genelist2,xyzp2]= sc_splinefit( X2(jrows, :), genelist3, sortit);
    xyz2 = [T2.lgu, T2.lgcv, T2.dropr];
    [xyz2, xyzp2] = minmax_scale(xyz2, xyzp2);

    %clear T2;
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
        xlabel('log1p(Mean)');
        ylabel('log1p(CV)');
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
    xyz = xyz(idx_tdv,1:3);

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
    xyz2 = xyz2(idx_tdv2,1:3);

    % (D1 - S1) - (D2 - S2)  or distance difference in vector form
    % (D1 - D2) + (S2 - S1)
    ddv = table2array( Tdv(:,2:4) - Tdv2(:,2:4) );
    % % Computing ( S1 - S2 )
    ddvs = table2array( Tdv(:,5:7) - Tdv2(:,5:7) );
    % % Computing + ( D1 - D2 )
    ddvp = xyz(:,1:3) - xyz2(:,1:3);
    % % Distance within splines + distance within points
    ddscale = vecnorm(ddvs,2,2) + vecnorm(ddvp,2,2);

    clear ddvs ddvp;

    dd0 = vecnorm(ddv,2,2);
    dd =  dd0./ddscale;

    % Standarize difference distance
    %ddz = zscore(dd0);
    ddz = zscore(dd);


    % Cumulative distribution function (cdf) according to
    % a standard normal distribution in d from dx sigma
    pval = 2*normcdf(ddz,'upper');


    % Table containing distance metrics of differentially variable genes
    Tdv_final = table(Tdv.genes, ...
                      Tdv.lgmu_dist, Tdv.lgcv_dist,...
                      Tdv.dropr_dist, xyz(:,3), Tdv.dist,...         
                      Tdv2.lgmu_dist, Tdv2.lgcv_dist, ...
                      Tdv2.dropr_dist, xyz2(:,3), Tdv2.dist, ...
                      ddv(:,1), ddv(:,2), ddv(:,3),...
                      ddscale, dd0,...
                      dd, ddz, pval);

    clear pval ddz dd dd0 ddv Tdv Tdv2;

    Tdv_final.Properties.VariableNames = {'genes', ...
                                        'lgmu_dist1', 'lgcv_dist1',...
                                        'dropr_dist1','drop1','dist_set1',...
                                        'lgmu_dist2', 'lgcv_dist2', ...
                                        'dropr_dist2','drop2','dist_set2', ...
                                        'lgmu_dd','lgcv_dd','dropr_dd',...
                                        'dd_scale', 'dd_raw',...
                                        'dd_norm', 'dd_std',...
                                        'pval'};

    % Sort together by pval acending and abs_dd decensindg
    idx = and(Tdv_final.drop1 > 0.01, Tdv_final.drop1 < 0.99); 
    Tdv_final = Tdv_final(idx,:);
    idx = and(Tdv_final.drop2 > 0.01, Tdv_final.drop2 < 0.99); 
    Tdv_final = Tdv_final(idx,:);

    [Tdv_final, ~] = sortrows(Tdv_final, 19, {'ascend'});
   
    % Writting differential variability table 
    filenameT = '_diff_var.csv';
    filenameT = strcat(fname,filenameT);
    writetable(Tdv_final,filenameT,'Delimiter',',');
    type(filenameT);
