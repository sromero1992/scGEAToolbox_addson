function sce_circ_phase_estimation(sce, tmeta,  rm_low_conf, period12, ...
                                    custom_genelist, custom_celltype)
    tic;
    % old_labels = ["1La" "2La" "3La" "4La" "5Da" "6Da" "7Da" "8Da"];
    % new_labels = ["ZT00" "ZT03" "ZT06" "ZT09" "ZT12" "ZT15" "ZT18" "ZT21"]';
    % tmeta = table( old_labels', new_labels');
    %custom_genelist = ["Clock","Arntl","Per1","Per2","Per3","Cry1","Cry2","Nr1d1","Nr1d2","Rora",...
    %              "Rorc","Sirt1","Bhlhe40","Bhlhe41","Timeless", "Xbp1", "Atf4", "Atf6", "Hif1a"];
    if nargin < 3; rm_low_conf = true; end
    if nargin < 4; period12 = false; end
    if nargin < 5; custom_genelist = []; end
    if nargin < 6; custom_celltype = []; end

    if period12 
        disp( "Circadian identification with 12 hrs period...")
    else
        disp( "Circadian identification with 24 hrs period...")
    end
   
    % Merge replicates for this analysis or just re-label
    batches = unique(sce.c_batch_id);

    % % Check if batches in sce match the old labes in meta table
    % right_input = all( batches == tmeta.old_labels );
    % if ~right_input 
    %     fprintf('Your metadata table does not match sce object batches \n' );
    %     fprintf('Please, ensure that same number of batches match old_input');
    %     return;
    % end

    % Rename batches 
    for ib = 1:length(batches)
        str_idx = find( batches(ib) == tmeta.old_labels );
        idx = find(sce.c_batch_id == batches(ib));
        sce.c_batch_id(idx) = tmeta.new_labels(str_idx);
    end
    batches = unique(sce.c_batch_id);
    
    disp("New batches ")
    disp(batches')

    % All cell types available (Make this optional any other time)
    cell_type_list = unique(sce.c_cell_type_tx);
    ncell_types = length(cell_type_list);

    % Compute R with cell percentage (true)
    % Compute R with gene expression (false)
    cell_pct = false;
    % N pseudo-bulk chunks
    n_bulk = 200; % 100 or 50
    % minimum number of cells to bootstrap pseudo-bulk
    ncellb = 50;% 100 or 500
    % Number of time points
    nztps = length(unique(tmeta.new_labels));
    
    info_per_cell_type = zeros(ncell_types, 3);
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
            gene_list = custom_genelist;
        end
        ngene = length(gene_list);
    
        % Gene names stored 
        gene_names = strings(ngene,1);
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
    
        R0 = zeros(ngene,nzts);
        acro = zeros(ngene,1);
    
        amp =zeros(ngene,1);
        T = zeros(ngene,1);
        mesor = zeros(ngene,1);
        failed = zeros(ngene,1);
        
        mae = zeros(ngene,1);
        mae_rel = zeros(ngene,1);
    
        parfor igene = 1:ngene
        %for igene = 1:ngene

            %prog = floor(igene/ngene*100);
            %textprogressbar(prog);
            % Gene index to work on from list
            gene_names(igene) = gene_list(igene);
            ig = find(sce_sub.g == gene_names(igene));
            % Set the same gene name for pseudo-bulks
            range = 1:n_bulk;
            range = range';
            R = zeros(n_bulk, nzts);

            % Pre-screening and computing the Gene exp for all time points
            % with all the cells
            for it = 1:nzts
                ic = find(sce_sub.c_batch_id == batch_time(it));
                % Gene expression for a time point
                Xg_tp = full(sce_sub.X(ig,ic));
                R0(igene, it) = compute_pseudoB_R( Xg_tp, 1, size(Xg_tp,2), cell_pct);
            end
    
            % Failure type 1, mean expression /cell pct 
            % is zero in all time points and don't compute anything
            if ~all(R0(igene, :)) 
                failed(igene) = 1;
                continue;
            end
            % Compute Pseudo-bulk gene expression by time point in batches of
            % ncellb if possible with Bootstrap method, else take half of
            % possible to make it noisy prediction and being filtered
            for it = 1:nzts
                ic = find(sce_sub.c_batch_id == batch_time(it));
                % Gene expression for a time point
                Xg_tp = full(sce_sub.X(ig,ic));
                % Bootstrap subsampling in batches of ncellb cells if possible
                % otherwise take half of the cell population to make it noisy
                if ncellb > size( Xg_tp, 2)
                    ncellb0 = ceil(size(Xg_tp, 2)/2);
                    % Consider one time point being unreliable: failure type 2
                    failed(igene) = 2;
                else
                    ncellb0 = ncellb;
                end
                R(range,it) = compute_pseudoB_R( Xg_tp, n_bulk, ncellb0, cell_pct);
            end
            % if failed skip this gene
            time_cycle = 21; % Based on Sato data
            time_step = 3; % Based on Sato data
        
            % Sine function allocations R is nbulks * nzts
            [acro(igene), amp(igene), T(igene), mesor(igene)] = ...
                    estimate_phaseR( R, time_cycle, time_step, period12);
        
            % Allocation for errors
            mae(igene) = 10e10;
            mae_rel(igene) = 10e10;
        
            t = 0:3:21;
            % Evaluate cos function for every time point cell percentage
            fval = zeros(1, nzts);
            fval(:) = amp(igene).*cos( 2*pi*( t - acro(igene))./T(igene) ) + mesor(igene);
        
            % Measure error
            abs_err = abs(R - fval);
            abs_err_rel = abs_err./abs(fval);
            abs_err = mean(abs_err, 1);
            mae(igene) = mean(abs_err, 2);
            
            % Measure relative error
            abs_err_rel = mean(abs_err_rel, 1);
            mae_rel(igene) = mean(abs_err_rel, 2);

            % Compare the derivative direction of fval and R0 if no match, true
            % also test if there is more than 30% of different directions in
            % slopes
            diff_bool = diff_test(fval,R0(igene, :));
            if diff_bool
                failed(igene) = 3;
                continue;
            end
        end
    
        acro_formatted = acro;
        acro_formatted(acro < 0) = acro(acro < 0) + 24;
        acro_formatted(acro > 24) = acro(acro > 24) - 24;
    
        T1 = table(gene_names, amp, abs(amp), mesor, acro, acro_formatted, T, mae, ...
                   mae_rel, failed );
        
        %clear amp mesor acro T mae rmse mae_rel rmse_rel failed;
        
        T1.Properties.VariableNames = ["Genes","Amp", "Abs_Amp", "Mesor","Acrophase", ...
                                       "Acrophase 24", "Period", ...
                                       "MAE", "MAE_rel" ,"Failed"];
        
        T2 = table(gene_names, R0(:,1), R0(:,2), R0(:,3), R0(:,4), ...
                               R0(:,5), R0(:,6), R0(:,7), R0(:,8));
        
        T2.Properties.VariableNames = ["Genes","ZT00","ZT03","ZT06","ZT09","ZT12","ZT15","ZT18","ZT21"];  
        
        % Remove failed with zero expression
        idx = T1.Failed ~= 1;
        T1 = T1(idx,:);
        T2 = T2(idx,:);
    
        % Remove non- significant pvalue predictions 
        [T1, conf_idx] = gen_pvalues( T1, 1);
        % remove low confidence genes?
        if rm_low_conf
            T1 = T1(conf_idx, :);
            T2 = T2(conf_idx, :);
        end
       
        % Remove 30% errors in direction or not same direction in end points
        idx = T1.Failed ~= 3;
        T1 = T1(idx,:);
        T2 = T2(idx,:);
    
        % Sort by acrophase and amplitude to classify better
        [T1, idx] = sortrows(T1,["pvalues","Acrophase","Abs_Amp"],...
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
    
        info_per_cell_type(icell_type,:) = [ sce_sub.NumCells, sce_sub.NumGenes, ...
                                             length(T1.Genes)];
        
    end
    T0 = table( cell_type_list, info_per_cell_type(:,1),info_per_cell_type(:,2), ...
                               info_per_cell_type(:,3));
    
    T0.Properties.VariableNames = ["CellType", "No. cells", "No. genes", "No. Circadian genes"];
    ftable_name = strcat("cell_genes_info",per_lab);
    ftable_name = strcat(ftable_name,"analysis.csv");
    writetable(T0, ftable_name);
    
    % %
    % t0 = 0:0.1:21;
    % t = 0:3:21;
    % igene = find(T1.Genes=="Tmsb10");
    % for i = igene
    %     fval = T1.Amp(i).*cos( 2*pi*( t0 - T1.Acrophase(i))./T1.Period(i) ) + T1.Mesor(i);
    %     Rzts = table2array( T2(i,2:end) );
    %     data_range = abs( max(Rzts) - min(Rzts) );
    %     plot(t0, fval);
    %     hold on;
    %     plot(t,Rzts);
    %     hold off;
    %     % scatter(t,R(1:100,:));
    %     % hold on;
    % end
    
    % inv_genes = readtable( "invasion_genes.txt", 'readvariablenames', false);
    % inv_genes = string(table2array(inv_genes));
    
    % for i = 1:length(inv_genes)
    %     targets = find(T1.Genes==inv_genes(i));
    %     targets
    % end
    time_end = toc;
    disp("Computational time: " + time_end)

end