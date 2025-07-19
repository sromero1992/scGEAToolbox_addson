 % Include the heatmap generation function here (or ensure it's in your MATLAB path)
    function generateHeatmap_circ_simple(sce, celltype, strict, customName, circ, norm_str)
        if nargin < 6 || isempty(norm_str); norm_str = 'lib_size'; end

        % Specify the file name for gene expression data
        fname = strcat(celltype, "_period_24__macro_circadian_analysis.csv");
        
        % Check if the analysis file exists before trying to read it
        if ~exist(fname, 'file')
            warning('Analysis file "%s" not found. Cannot generate heatmap.', fname);
            return; % Exit the function if the file doesn't exist
        end
        
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
                            "nfkb" "per" "wnt" "nrd" "rev" "pik"]; % Added "nrd" and "rev" as they appear in your original code
            glcirc = contains(lower(Dwork.Genes), classic_circ, 'IgnoreCase', true);
            Dwork = Dwork(glcirc, :);
        end
        if height(Dwork) <= 1
            disp('Small number of predicted circadian genes...')
            disp(height(Dwork))
            return;
        end

        % Normalize full set
        %X = full(sce.X);
        if strcmp(norm_str, 'lib_size')
            % This is regular cells pipeline
            %X = sc_norm(full(sce.X));
            X = pkg.norm_libsize(sce.X, 1e4);
            X = log1p(sce.X);
        else % 'magic_impute'
            % This for cancer cells
            X = sc_impute(sce.X, 'MAGIC');
        end
        X = sparse(X);
        sce.X = X;
        clear X;

        % Extract circadian gene list for specified cell type
        gl = string(Dwork.Genes);
        ic = strcmpi(sce.c_cell_type_tx, celltype);
        ig = ismember(lower(sce.g), lower(gl));
        % Extract expression data for circadian genes and specified cell type
        X = sce.X(ig, ic);
        g = sce.g(ig);

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
                % Ensure batch ID comparison is robust to type (char/string)
                ic = ismember(sce_new.c_batch_id, time_batches(ib));
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
        % We need to match the genes in Dwork with the genes in sce_new
        [~, sort_idx] = ismember(lower(Dwork.Genes), lower(sce_new.g));

        % Remove indices that were not found (where ismember returned 0)
        found_genes_idx = sort_idx > 0;
        sort_idx = sort_idx(found_genes_idx);
        sorted_gene_names = Dwork.Genes(found_genes_idx); % Keep track of the names

        if isempty(sort_idx)
            disp('No matching genes found in SCE object after filtering.');
            return;
        end

        % Use sort_idx to extract data_sorted and corresponding gene names
        data_sorted = data_scaled(sort_idx, :);
        sorted_gene_names_sce = sce_new.g(sort_idx); % Get the gene names from sce_new based on sort_idx

        % Double check if the gene names match after indexing (should if ismember worked)
        if ~isequal(lower(sorted_gene_names), lower(sorted_gene_names_sce))
             warning('Gene name mismatch after indexing. Proceeding with indexing results.');
             % You might want to inspect sorted_gene_names and sorted_gene_names_sce here
        end

        % Generate heatmap
        figureHandle = figure;
        h = imagesc(data_sorted);
        set(h, 'AlphaData', ~isnan(data_sorted));  % Handle NaN values gracefully

        % Set axis labels and customize appearance
        axesHandle = get(h, 'Parent');
        axesHandle.XTick = 1:length(time_batches);
        axesHandle.XTickLabels = string(time_batches);
        axesHandle.YTick = 1:length(sorted_gene_names_sce); % Set YTicks for each gene
        axesHandle.YTickLabels = sorted_gene_names_sce;    % Set YTickLabels to sorted gene names
        
        % Set axis labels and customize appearance
        if size(data_sorted, 1) > 20
            axesHandle.YTick = [];
            axesHandle.YTickLabels = [];
        end
        
        title(sprintf('Mean Gene Expression Heatmap for %s', celltype));
        xlabel('Time Batches');
        ylabel('Genes (Sorted by Acrophase)');

        % Define custom blue-white-red colormap
        n = 256;
        blueToWhite = [linspace(0, 1, n/2)', linspace(0, 1, n/2)', linspace(1, 1, n/2)'];
        whiteToRed = [linspace(1, 1, n/2)', linspace(1, 0, n/2)', linspace(1, 0, n/2)'];
        customColormap = [blueToWhite; whiteToRed];

        colormap(customColormap);
        % Set caxis based on the actual data range
        data_range = [min(data_sorted(:)), max(data_sorted(:))];
        if diff(data_range) > 0 % Avoid error if all values are the same
             caxis(data_range);
        end
        colorbar;
        set(gca, 'FontSize', 10);
        axis tight;
        % Construct the filename based on `customName`, `clus_method`, `strict`, and `circ`
        fname_prefix = strcat(celltype);
        if ~isempty(customName)
            fname_prefix = strcat(fname_prefix, '_', customName);
        end
        fname_prefix = strcat(fname_prefix, '_', clus_method);

        if strict
            fname_prefix = strcat(fname_prefix, '_strict');
        else
            fname_prefix = strcat(fname_prefix, '_not_strict');
        end

        if circ
            fname_prefix = strcat(fname_prefix, '_core_circ');
        end

        % Specify the output file name for the table (using Dwork which is already filtered)
        output_table_fname = strcat(fname_prefix, "_Dwork_analysis_results.csv");
        writetable(Dwork, output_table_fname);

        % Save figure
        saveas(figureHandle, strcat(fname_prefix, '.svg'));
        close(figureHandle);
    end