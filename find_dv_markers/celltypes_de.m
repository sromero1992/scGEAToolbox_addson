function celltypes_de(sce, comparisons, fname)
    % INPUTS:
    % sce -----------------> SCE object
    % fname ---------------> file name to save processed SCEs
    % comparisions --------> cell object with two class labes
    %                         (cell comparisons = {'08w', '14w'};)
    % qc_true -------------> boolean to subsample and do light qc
    % OUTPUT:
    % Files...
    % Example: find_dv_markers(sce,'example_mouse',500,true)
    if nargin <3 || isempty(fname)
        fname = 'results';
    end
    cell_types = unique(sce.c_cell_type_tx);
    % Assign a boolean to normalize
    X = sc_norm(sce.X,'type','libsize');
    X = log1p(X);
    tic
    fprintf("Starting DE analysis for clusters...\n");
    
    for icells = cell_types'
        disp(["Working on: " icells])
        
        % Clean the cell type name to create a safe file and sheet name
        icells_safe = replace(icells, {' ', ':', '\', '/', '?', '*', '[', ']'}, '_');
        
        % Define the unique Excel filename for this specific cell type here, inside the loop
        fname_excel_celltype = sprintf('%s_%s_DE.xlsx', fname, icells_safe);

        % Get ith-cluster cells
        cell_idx = icells == sce.c_cell_type_tx;
        idx1 = comparisons{1} == sce.c_batch_id;
        idx2 = comparisons{2} == sce.c_batch_id;
    
        X1 = X(:, and(idx1, cell_idx) );
        X2 = X(:, and(idx2, cell_idx) );
        
        % Call the DE function 
        [Tde, Tup, Tdn] = sc_deg(X1, X2, sce.g, 1, false);
        
        % Save DE results to different sheets in the unique Excel file
        % We use 'fname_excel_celltype' consistently for all calls.
        sheetName_all = 'All_Genes';
        writetable(Tde, fname_excel_celltype, 'Sheet', sheetName_all);
        
        sheetName_up = 'Up_Genes';
        writetable(Tup, fname_excel_celltype, 'Sheet', sheetName_up);
        
        sheetName_dn = 'Down_Genes';
        writetable(Tdn, fname_excel_celltype, 'Sheet', sheetName_dn);
        
        % Pass the same filename to the enrichment function
        gui.e_enrichrxlsx(Tup, Tdn, Tde, fname_excel_celltype);
    end
    time_end_de = toc;
    fprintf("DE analysis for clusters time: %f sec \n", time_end_de);
end