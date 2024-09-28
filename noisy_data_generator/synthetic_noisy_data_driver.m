%addpath('../src-v0.2/');
rng default;
p = 10000; %observations
n = 50; %features

% Source features (lineal)
source_f_linear = [5, 11, 7];
target_f_linear = [16, 17, 18];
linear_feat = union(source_f_linear, target_f_linear);       

% Source features (nonlinear)
source_f_nlinear = [1, 14, 8];
target_f_nlinear = [19, 20, 22];
nlinear_feat = union(source_f_nlinear, target_f_nlinear);       

% Source features are highly correlated to target features 
% where target is "regulated" source
[p,n,X,Yl,Ynl] = synthetic_data2(p, n, source_f_linear, target_f_linear, ...
                                  source_f_nlinear, target_f_nlinear);

% imagesc(X(1:25,1:25))
% % Features to extract
% K = 12;
% g = 1:n;
% X = X';
% Yl = Yl';
% Ynl = Ynl';
% 
% %Tqubol = qfeatures_qubo_base(X, g, Yl, K, false);
% %Tqubonl = qfeatures_qubo_base(X, g, Ynl, K, false);
% %inter_feat_qubo_l = intersect(Tqubol.sol_genes, linear_feat)
% %inter_feat_qubo_nl = intersect(Tqubonl.sol_genes, nlinear_feat)
% 
% Tmll  = mlfeatures_base(X, g, Yl, K, 2);
% Tmlnl  = mlfeatures_base(X, g, Ynl, K, 2);
% 
% inter_feat_lasso_l = intersect(Tmll.sol_genes_lasso, linear_feat)
% inter_feat_lasso_nl = intersect(Tmlnl.sol_genes_lasso, nlinear_feat)
% 
% inter_feat_relief_l = intersect(Tmll.sol_genes_relief, linear_feat)
% inter_feat_relief_nl = intersect(Tmlnl.sol_genes_relief, nlinear_feat)
