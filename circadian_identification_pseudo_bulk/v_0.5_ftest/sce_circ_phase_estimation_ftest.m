function [T1, T2] = sce_circ_phase_estimation_ftest(sce, tmeta,  rm_low_conf, period12, ...
                                    custom_genelist, custom_celltype)
    tic;
    rng('default');
    % old_labels = ["1La" "2La" "3La" "4La" "5Da" "6Da" "7Da" "8Da"]';
    % new_labels = ["ZT00" "ZT03" "ZT06" "ZT09" "ZT12" "ZT15" "ZT18" "ZT21"]';
    % times = [0 3 6 9 12 15 18 21]';
    % tmeta = table( old_labels, new_labels, times);
    %custom_genelist = ["Clock","Arntl","Per1","Per2","Per3","Cry1","Cry2","Nr1d1","Nr1d2","Rora",...
    %              "Rorc","Sirt1","Bhlhe40","Bhlhe41","Timeless", "Xbp1", "Atf4", "Atf6", "Hif1a"];
    if nargin < 3 || isempty(rm_low_conf); rm_low_conf = true; end
    if nargin < 4 || isempty(period12); period12 = false; end
    if nargin < 5 || isempty(custom_genelist); custom_genelist = []; end
    if nargin < 6 || isempty(custom_celltype); custom_celltype = []; end

    if period12 
        disp( "Circadian identification with 12 hrs period...")
    else
        disp( "Circadian identification with 24 hrs period...")
    end
   

    % Merge replicates for this analysis or just re-label
    batches = unique(sce.c_batch_id);

    time_cycle = max(tmeta.times);
    time_step = mean( diff(tmeta.times) ); % Based on Sato data
    % Rename batches 
    for ib = 1:length(batches)
        str_idx = find( batches(ib) == tmeta.old_labels );
        idx = find(sce.c_batch_id == batches(ib));
        sce.c_batch_id(idx) = tmeta.new_labels(str_idx);
    end
    batches = unique(sce.c_batch_id);
    
    disp("New batches")
    disp(batches')

    % All cell types available (Make this optional any other time)
    cell_type_list = unique(sce.c_cell_type_tx);
    ncell_types = length(cell_type_list);

    % Number of time points
    nztps = length(unique(tmeta.new_labels));
    
    info_p_type = zeros(ncell_types, 5);
    for icell_type = 1:ncell_types   
        % Extract count matrix for ith cell type
        cell_type = cell_type_list(icell_type);
        if ~isempty(custom_celltype)
            if ~ismember(cell_type, custom_celltype); continue; end
        end
        % Count matrix only for cell_type
        idx = find(sce.c_cell_type_tx == cell_type);
        X = sce.X(:,idx);
        X = full(X);
        fprintf("Processing cell type %s \n", cell_type);
        % Normalizing count data for that cell type / (pearson residuals?)
        X = sc_norm(X);
        X = sparse(X);
        %X = sc_transform(X,"type","PearsonResiduals"); % This is not best for Relative error since gets 0 fit values
        
        % Gene genes and cells information
        g = sce.g;
        cell_batch = sce.c_batch_id(idx);
        
        sce_sub = SingleCellExperiment(X,g);
        sce_sub.c_batch_id = cell_batch;
        sce_sub.c_cell_type_tx = sce.c_cell_type_tx(idx);
        % non-stringent filtering (not filtering anymore)
        %sce_sub.qcfilter(500, 0.20, 10);
        clear X g;
        
        %gene_list = ["Clock","Arntl","Per1","Per2","Per3","Cry1","Cry2","Nr1d1","Nr1d2","Rora",...
        %              "Rorc","Sirt1","Bhlhe40","Bhlhe41","Timeless", "Xbp1", "Atf4", "Atf6", "Hif1a"];
        if isempty( custom_genelist)
            disp("Circadian analysis for all genes");
            gene_list = sce_sub.g;
        else
            disp("Circadian analysis for custom genes");
            nx = size(custom_genelist, 1);
            if nx == 1; custom_genelist = custom_genelist'; end
            [lic, ~] = ismember(custom_genelist, sce_sub.g);
            custom_genelist = custom_genelist(lic);
            gene_list = custom_genelist;

        end
        ngene = length(gene_list);

        % Batch name for each time point 
        batch_time = unique(sce_sub.c_batch_id);
        % Number of time points
        nzts = length(batch_time);
        fprintf("Number of NZTS after sub-sample: %d \n", nzts);
        
        % If number of time points do not match experimental time points, dump cell type
        if nzts ~= nztps
            disp("Number of time points does not match to input metadata table")
            continue; 
        end
   
        % Initialize temporary arrays to store intermediate results
        tmp_acro = zeros(ngene, 1);
        tmp_amp = zeros(ngene, 1);
        tmp_T = zeros(ngene, 1);
        tmp_p_value = zeros(ngene, 1);
        tmp_R0 = zeros(ngene, nzts);

        tmp_mesor = zeros(ngene, 1);
        for igene = 1:ngene
             prog = floor(igene/ngene*100);
             textprogressbar(prog);
        %parfor igene = 1:ngene
            % Gene index to work on from list
            ig = find(sce_sub.g == gene_list(igene));
        
            Xg_zts = {};
            % Prepare gene expression per time point
            for it = 1:nzts
                ics = find(sce_sub.c_batch_id == batch_time(it));
                % Gene expression for a time point
                Xg_zts{it} = full(sce_sub.X(ig, ics));
                tmp_R0(igene, it) = mean(Xg_zts{it});
            end 

            [tmp_acro(igene), tmp_amp(igene), tmp_T(igene), tmp_mesor(igene), tmp_p_value(igene)] = ...
                            estimate_phaseR_Ftest(Xg_zts, time_step, period12);       

        end
        
        % Aggregate results from temporary arrays to the main arrays
        acro = tmp_acro;
        amp = tmp_amp;
        T = tmp_T;
        mesor = tmp_mesor;
        R0 = tmp_R0;
        p_value = tmp_p_value;
        
        clear tmp_mesor tmp_T tmp_amp tmp_acro;

        acro_formatted = acro;
        acro_formatted(acro < 0) = acro(acro < 0) + 24;
        acro_formatted(acro > 24) = acro(acro > 24) - 24;
    
        T1 = table(gene_list, amp, abs(amp), mesor, acro, ...
                                      acro_formatted, T, p_value);
        
        %clear amp mesor acro T mae rmse mae_rel rmse_rel failed;
        
        T1.Properties.VariableNames = ["Genes","Amp", "Abs_Amp", "Mesor","Acrophase", ...
                                       "Acrophase_24", "Period", "pvalue"];
        
        T2 = table(gene_list, R0(:,1), R0(:,2), R0(:,3), R0(:,4), ...
                               R0(:,5), R0(:,6), R0(:,7), R0(:,8));
        
        T2.Properties.VariableNames = ["Genes","ZT00","ZT03","ZT06","ZT09","ZT12","ZT15","ZT18","ZT21"];  


        % Remove NaNs (not optional)
        rm_nans_idx = ~isnan(T1.pvalue) ;
        T1 = T1(rm_nans_idx,:);
        T2 = T2(rm_nans_idx,:);

        % remove low confidence genes?
        rm_nconfs_idx = T1.pvalue < 0.05;
        num_conf_g = sum(rm_nconfs_idx);
        num_n_conf_g = sum(~rm_nconfs_idx);
        if rm_low_conf
            T1 = T1(rm_nconfs_idx, :);
            T2 = T2(rm_nconfs_idx, :);
        end
       
        % Sort by acrophase and amplitude to classify better
        [T1, idx] = sortrows(T1,["pvalue","Acrophase","Abs_Amp"],...
                                {'ascend','ascend','descend'});
        
        T2 = T2(idx,:);
        if period12
            per_lab = "_period_12_";
        else
            per_lab = "_period_24_";
        end
    
        fname = strcat(cell_type, per_lab); % Run this block if you want a new file
        ftable_name = strcat(fname,"_macro_circadian_analysis.csv");
        writetable(T1,ftable_name);
        
        ftable_name = strcat(fname,"_macro_circadian_ZTs.csv");
        writetable(T2,ftable_name);
    
        info_p_type(icell_type,:) = [ sce_sub.NumCells, sce_sub.NumGenes, ...
                                             length(T1.Genes), num_conf_g, num_n_conf_g];
        
    end


    T0 = table( cell_type_list, info_p_type(:,1), info_p_type(:,2), ...
                info_p_type(:,3), info_p_type(:,4), info_p_type(:,5));

    T0.Properties.VariableNames = ["CellType", "Num. cells", "Num. genes", ...
                                   "Num. circadian genes", ...
                                   "Num. confident genes", ...
                                   "Num. not confident genes"];

    ftable_name = "cell_genes_info";
    if ~isempty(custom_celltype)
        ftable_name = strcat(ftable_name, custom_celltype);
    end
    ftable_name = strcat(ftable_name,per_lab);
    ftable_name = strcat(ftable_name,"analysis.csv");
    writetable(T0, ftable_name);
    
    time_end = toc;
    disp("Computational time: " + time_end)

end