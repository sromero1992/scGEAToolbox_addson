function generateHeatmap_circ_simple(sce, celltype, strict, customName, circ)
    % Specify the file name for gene expression data
    fname = strcat(celltype, "_period_24__macro_circadian_analysis.csv");
    D = readtable(fname, 'ReadVariableNames', true);

    % Sort rows based on acrophase and amplitude
    D = sortrows(D, ["Acrophase_24", "Abs_Amp"], ["ascend", "descend"]);

    % Filter adjusted p-values and amplitude
    if strict
        Dwork = D(D.pvalue_adj <= 0.05, :);
    else
        Dwork = D(D.pvalue <= 0.05, :);
    end
    Dwork = Dwork(Dwork.Amp >= 0, :);

    % Filter for circadian genes if specified
    if circ 
        classic_circ = ["arn" "bhlh" "clock" "cry" "dbp" "tef" "hlf" "raf" "erk" ...
                        "mek" "ras" "mtor" "map" "ral" "akt" "hif" "kras" "myc" ...
                        "nfkb" "per" "wnt" "nrd" "rev" "pik"];
        glcirc = contains(lower(Dwork.Genes), classic_circ, 'IgnoreCase', true);
        Dwork = Dwork(glcirc, :);
    end
    if height(Dwork) <= 1
        disp('Small number of predicted circadian genes...')
        disp(height(Dwork))
        return;
    end

    % Extract circadian gene list for specified cell type
    gl = string(Dwork.Genes);
    ic = strcmpi(sce.c_cell_type_tx, celltype);  
    ig = ismember(lower(sce.g), lower(gl));  

    % Extract expression data for circadian genes and specified cell type
    X = sce.X(ig, ic);
    g = sce.g(ig);

    % Normalize data and create new SingleCellExperiment object
    X = sc_norm(X);
    sce_new = SingleCellExperiment(X, g);
    sce_new.c_batch_id = sce.c_batch_id(ic);
    sce_new.c_cell_type_tx = sce.c_cell_type_tx(ic);
    sce_new.c_cell_id = sce.c_cell_id(ic);

    % Calculate mean expression across batches
    time_batches = unique(sce_new.c_batch_id); 
    nt = length(time_batches);
    ng = length(sce_new.g);
    gzts = zeros(nt, ng);

    for igs = 1:ng
        for ib = 1:nt
            ic = sce_new.c_batch_id == time_batches(ib);
            gzts(ib, igs) = mean(sce_new.X(igs, ic), 'omitnan');
        end
    end

    % Scale data with chosen normalization method
    norm_type = 'zscore';  % Choose normalization type: 'zscore', 'norm', 'scale', 'center'
    switch norm_type
        case 'zscore'
            data_scaled = normalize(gzts, "zscore", "std");
        case 'norm'
            data_scaled = normalize(gzts, "norm", 2);
        case 'scale'
            data_scaled = normalize(gzts, "range", [-1, 1]);
        case 'center'
            data_scaled = normalize(gzts, "center", "mean");
        otherwise
            data_scaled = normalize(gzts, "zscore", "std");
    end
    data_scaled = data_scaled';
  
    % Set clustering method
    clus_method = 'none';  % Fall back to no clustering

    % Get index order based on final sort by "Acrophase_24"
    sort_idx = NaN(length(Dwork.Genes), 1);  % Initialize sort_idx with NaN
    gl = string(Dwork.Genes);  % Ensure gl contains the gene names from Dwork
    
    for ig = 1:length(gl)
        idx = find(sce_new.g == gl(ig), 1);  % Find the index of the gene
        if ~isempty(idx)  % Check if the gene was found
            sort_idx(ig) = idx;  % Assign index if found
        else
            % Warn if not found
            warning('Gene %s not found in the SingleCellExperiment object.', gl(ig));
            sort_idx(ig) = NaN;  % Assign NaN or handle appropriately
        end
    end
    
    % Remove NaN values from sort_idx
    sort_idx(isnan(sort_idx)) = [];
    
    % Use sort_idx to extract data_sorted
    data_sorted = data_scaled(sort_idx, :);  
    
    % Generate heatmap
    figureHandle = figure;
    h = imagesc(data_sorted);
    set(h, 'AlphaData', ~isnan(data_sorted));  % Handle NaN values gracefully
    
    % Set axis labels and customize appearance
    axesHandle = get(h, 'Parent');
    axesHandle.XTick = 1:length(time_batches);
    axesHandle.XTickLabels = string(time_batches);
    axesHandle.YTick = [];               % Remove y-axis tick marks
    axesHandle.YTickLabels = [];          % Remove y-axis tick labels
    
    title(sprintf('Mean Gene Expression Heatmap for %s', celltype));
    xlabel('Time Batches');
    ylabel('Gene Number');
    
    % Define custom blue-white-red colormap
    n = 256;
    blueToWhite = [linspace(0, 1, n/2)', linspace(0, 1, n/2)', linspace(1, 1, n/2)'];
    whiteToRed = [linspace(1, 1, n/2)', linspace(1, 0, n/2)', linspace(1, 0, n/2)'];
    customColormap = [blueToWhite; whiteToRed];
    
    colormap(customColormap);
    caxis([min(data_sorted(:)), max(data_sorted(:))]);
    colorbar;
    set(gca, 'FontSize', 10);
    axis tight;

    % Construct the filename based on `amp_pos` and `circ` parameters
    customName = strcat(customName, '_', clus_method);
    
    if strict
        fname = strcat(celltype, '_', customName, '_strict_filter_');
    else
        fname = strcat(celltype, '_', customName, '_not_strict_filter_');
    end
    
    if circ
        fname = strcat(fname, '_core_circ_');
    end
    
    % Specify the output file name for the table
    output_fname = strcat(fname, "_with_clusters.csv");
    writetable(Dwork, output_fname);
    
    % Maximize and save figure
    screenSize = get(0, 'ScreenSize');
    figureWidth = screenSize(3) * 0.2;  
    figureHeight = screenSize(4) * 0.9; 
    
    set(figureHandle, 'Position', [(screenSize(3) - figureWidth) / 2, screenSize(4) * 0.05, figureWidth, figureHeight]);
    saveas(figureHandle, strcat(fname, '.svg'));
    close(figureHandle);
end
