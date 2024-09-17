function gene_qc_test(sce, gene_list, plotit, thresholds)
    % gene_qc_test tests the distribution of gene expression and joint 
    % distributions of gene pairs in single-cell data.
    %
    % INPUT:
    %   - sce: A SingleCellExperiment object containing the data matrix X (genes x cells) 
    %          and metadata fields like sce.g (gene names).
    %   - gene_list: A cell array containing the gene names to be tested.
    %                ( e.g. gene_list = {'Msh2' 'tdTomato'}; )
    %   - plotit: A boolean to determine whether to plot the results (default = false).
    %   - thresholds: An array of thresholds for filtering gene counts (default = 
    %                 [0, 1, 2, 3, 4, 5, 10, 15, 20, 30, 50, 100, 1000]).
    %
    % OUTPUT:
    %   None. Prints AIC and p-values for each gene and joint gene pair in the console.
    % USAGE: 
    % gene_list = {'Msh2' 'tdTomato'}; 
    % gene_qc_test(sce, gene_list, true, 0:3)
    if nargin < 3; plotit = false; end
    if nargin < 4
        thresholds = [0, 1, 2, 3, 4, 5, 10, 15, 20, 30, 50, 100, 1000];
    end

    ngenes = length(gene_list);
    if ngenes < 1
        error("Number of genes is less than 1...");
    end
    fprintf("Number of genes to test: %d \n", ngenes);

    % Test individual genes
    for igene = 1:ngenes
        loc_gene = gene_list{igene};
        disp("--------------------------------------------");
        disp("Testing " + loc_gene);
        nb_test_per_gene(sce, loc_gene, thresholds, plotit);
    end 

    % Test joint distribution of gene pairs
    disp("--------------------------------------------");
    disp("------------- Joint gene test --------------");
    disp("--------------------------------------------");

    for igene = 1:ngenes
        for jgene = igene + 1:ngenes
            iloc_gene = gene_list{igene};
            jloc_gene = gene_list{jgene};
            disp("--------------------------------------------");
            fprintf("Testing pair %s - %s \n", iloc_gene, jloc_gene);
            nb_joint_test_per_pair(sce, iloc_gene, jloc_gene, thresholds, plotit);
        end
    end
end
