%data = load('knockdown_experiments_qc_stringent_decountX.mat');
data = load('dim_experiments_qc_stringent_decountX.mat');
sce = data.sce;
clear data;

my_gene = ["NR4A1" "NR4A2" "NR4A3" "EHMT2" "LEF1" "AXIN2" "VIM" ...
    "ZEB1"  "CDH1"  "CDH3"];

X = full(sce.X);
X = sc_norm(X);
X = log1p(X);

batches = unique(sce.c_batch_id);
% % Knockdown experiments array
% tmp = batches(2);
% batches(2) = batches(3);
% batches(3) = batches(4);
% batches(4) = tmp;

% Dim treatment experiments array 
tmp = batches(2);
batches(2) = batches(3);
batches(3) = tmp;

all_cells = size(X,2);
ng = length(my_gene);
nbatch = length(batches);

% Get gene index from sce
geneidx = zeros(ng,1);
for ig = 1:ng
    geneidx(ig) = find( my_gene(ig) == sce.g);
end

Xheat = zeros(ng, nbatch);

for ib = 1:nbatch
    cellidx = find( sce.c_batch_id == batches(ib) );
    Xtmp = X(geneidx, cellidx);
    % Mean of each row/gene
    Xheat(1:ng,ib) = mean(Xtmp,2);
end

%-------- MI heatmap-------------
h = imagesc(Xheat);

rangex = 1:nbatch;
rangey = 1:ng;
% Define string labels
xLabels = string(batches(rangex));
yLabels = string(my_gene(rangey));
axesHandle = get(h, 'Parent');
% Increse X ticks
axesHandle.XTick = rangex; 
axesHandle.YTick = rangey; 
% Set the x-axis tick labels
axesHandle.XTickLabels = xLabels;
axesHandle.YTickLabels = yLabels;

title('Gene average expression heatmap');
colorbar;
colormap('parula');
%------------------------------
