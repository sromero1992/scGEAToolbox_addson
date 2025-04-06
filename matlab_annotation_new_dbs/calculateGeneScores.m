function scores = calculateGeneScores(Tsub, genes)
    % calculateGeneScoresCVersion Calculates gene scores based on gene presence in table Tsub.
    %
    %   scores = calculateGeneScoresCVersion(Tsub, genes) calculates scores for
    %   each gene in 'genes' based on its presence in the columns of table 'Tsub'.
    %
    %   Input:
    %       Tsub : Table containing gene information.
    %       genes : Cell array or string array of unique genes.
    %
    %   Output:
    %       scores : Matrix of gene scores.

    fprintf("Pre-scoring unique genes of database... \n");

    ntot_genes = length(genes);
    ndb = size(Tsub, 1);
    scores = zeros(ntot_genes, ndb);

    for idb = 1:ndb
        rowGenes = table2array(Tsub(idb, :)); % Get all genes in the current row
        bool = ismember(genes, rowGenes);
        scores(:, idb) = double(bool); % Convert logical to double
    end

    for ig = 1:ntot_genes
        val = sum(scores(ig, :));
        if val ~= 0
            scores(ig, :) = scores(ig, :) ./ val;
        end
    end

    %clear Ttmp bool; %Ttmp not created in this version.
end