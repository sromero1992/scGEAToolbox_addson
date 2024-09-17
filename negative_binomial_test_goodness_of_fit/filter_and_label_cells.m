function [sce, cell_labels]= filter_and_label_cells(sce, gene1, gene2, ntests)
    % filter_and_label_cells filters cells based on gene expression thresholds,
    % calculates cell percentages, and labels cells based on gene expression patterns.
    %
    % INPUT:
    %   - sce: A SingleCellExperiment object with gene expression data.
    %   - gene1: The name of the first gene (e.g., 'tdTomato').
    %   - gene2: The name of the second gene (e.g., 'Msh2').
    %   - ntests: The number of threshold tests to run.
    %
    % OUTPUT:
    %   - cell_labels: A 1D array of strings containing labels 'mode1_Msh2_tdTomato'
    %                  and 'negative_model1_Msh2_tdTomato' for filtered cells.
    % USAGE:
    % cells = filter_and_label_cells(sce,'Msh2', 'tdTomato', 10);

    % Labes for identified cells
    test_label0 = 'mode1_Msh2_tdTomato';
    test_label1 = 'negative_model1_Msh2_tdTomato';

    % Get indices for the two genes
    igene = find(strcmp(sce.g, gene1));
    jgene = find(strcmp(sce.g, gene2));

    % Initialize array to store cell percentages
    cell_pct = zeros(ntests, 8);
    Xg1 = sce.X(igene, :);
    Xg2 = sce.X(jgene, :);
    ntot_cells = length(Xg1);

    % Loop through thresholds and compute percentages
    for ithrs = 1:ntests + 1
        idx = Xg1 >= ithrs - 1;
        cell_pct(ithrs, 1) = sum(idx) / ntot_cells;
        cell_pct(ithrs, 2) = sum(~idx) / ntot_cells;

        jdx = Xg2 >= ithrs - 1;
        cell_pct(ithrs, 3) = sum(jdx) / ntot_cells;
        cell_pct(ithrs, 4) = sum(~jdx) / ntot_cells;

        % + +
        kdx = and(idx, jdx);
        cell_pct(ithrs, 5) = sum(kdx) / ntot_cells;
        cell_pct(ithrs, 6) = sum(~kdx) / ntot_cells;

        % + / -
        mdx = or(idx, jdx);
        cell_pct(ithrs, 7) = sum(mdx) / ntot_cells;
        cell_pct(ithrs, 8) = sum(~mdx) / ntot_cells;
    end
    disp(cell_pct)
    % Apply specific threshold conditions to label cells
    idx = Xg1 >= 1;
    jdx = Xg2 <= 1;
    mdx = and(idx, jdx);

    % Initialize cell_labels as a string array
    cell_labels = strings(ntot_cells, 1);

    % Label cells based on mdx
    cell_labels(mdx) = test_label0;
    cell_labels(~mdx) = test_label1;

    % Output the cell labels as a 1D array
    cell_labels = cell_labels(mdx | ~mdx);

    % Creating new attribute in sce object
    ncell_att = length(sce.list_cell_attributes);
    sce.list_cell_attributes{ncell_att + 1} = 'Expression-activation';
    sce.list_cell_attributes{ncell_att + 2} = cell_labels;

end
