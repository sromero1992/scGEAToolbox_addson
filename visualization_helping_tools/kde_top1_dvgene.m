% Initialization
target = ["ANPEP", "DDX60", "NRXN3", "LGR5"]; % Target genes
cell_types = unique(sce.c_cell_type_tx);

% Epithelial cells
cell_select = cell_types(2); % Select epithelial cell type
idx = strcmp(sce.c_cell_type_tx, cell_select); % Use strcmp for string comparison
sce_sub = sce.selectcells(idx);

% Separate by batches
batches = unique(sce.c_batch_id);
sce_sub_b = {};
for ib = 1:length(batches)
    idx = sce_sub.c_batch_id == batches(ib); % Logical indexing for batch selection
    sce_sub_b{ib} = sce_sub.selectcells(idx);
    sce_sub_b{ib} = sce_sub_b{ib}.qcfilter; % Apply quality control filter
    sce_sub_b{ib}.X = log1p( sc_norm(sce_sub_b{ib}.X) ) ;
    %sce_sub_b{ib}.X = sc_norm(sce_sub_b{ib}.X);

end

% Collect Xg with metadata
Xg = {};
jg = 1;
for ib = 1:length(batches)
    for ig = 1:length(target)
        idx = find(strcmp(sce_sub_b{ib}.g, target(ig))); % Use strcmp for string comparison
        if ~isempty(idx) % Ensure the gene exists in the subset
            Xg{jg}.data = sce_sub_b{ib}.X(idx, :); % Gene expression matrix for the target gene
            Xg{jg}.batch = batches(ib); % Store the batch metadata
            Xg{jg}.gene = target(ig);   % Store the target gene metadata
            jg = jg + 1;
        end
    end
end

ig1 = 1;
ig2 = 5;
ncells1 = sce_sub_b{1}.NumCells;
ncells2 = sce_sub_b{2}.NumCells;

batch1 = Xg{ig1}.batch;
batch2 = Xg{ig2}.batch;
gene1 = Xg{ig1}.gene;
gene2 = Xg{ig2}.gene;

assert(gene2 == gene1) % Ensure gene names match

rng(0, "twister"); % For reproducibility

% Assuming Xg{ig1}.data and Xg{ig2}.data are the two datasets
data1 = full(Xg{ig1}.data); % Dataset 1 (example from Xg{ig1})
data2 = full(Xg{ig2}.data); % Dataset 2 (example from Xg{ig2})

% Perform KDE estimation for both datasets
[fp1, xfp1] = kde(data1); % KDE for data1 (pdf estimate)
[fp2, xfp2] = kde(data2); % KDE for data2 (pdf estimate)

% Calculate the mean of both datasets
mean1 = mean(data1); % Mean for data1
mean2 = mean(data2); % Mean for data2

% Plot the estimated pdf for both datasets
figure;
hold on;

% Plot KDE for data1
h1 = area(xfp1, fp1, 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'LineWidth', 2); % Blue area for data1
% Plot KDE for data2
h2 = area(xfp2, fp2, 'FaceColor', 'red', 'FaceAlpha', 0.5, 'LineWidth', 2); % Red area for data2

% Plot the mean values as vertical lines with darker shades of blue and red
h3 = plot([mean1 mean1], ylim, 'Color', [0, 0, 0.5], 'LineWidth', 2); % Dark blue line for data1
h4 = plot([mean2 mean2], ylim, 'Color', [0.5, 0, 0], 'LineWidth', 2); % Dark red line for data2

% Labels and title with dynamic gene name
xlabel('X-axis (Data Values)');
ylabel('Density');
title(['Kernel Density Estimate (KDE) for Gene: ' gene1]); % Include gene name in title

% Adjust the legend to handle area and plot objects
legend([h1, h2, h3, h4], ...
       {strcat('Data1 KDE (Batch: ', batch1, ')'), ...
       strcat('Data2 KDE (Batch: ', batch2, ')'), ...
       'Mean Data1', ...
       'Mean Data2', 'Location', 'Best'});

hold off;
