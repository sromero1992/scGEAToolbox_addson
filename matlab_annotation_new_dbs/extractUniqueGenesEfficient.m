function uniqueGenes = extractUniqueGenesEfficient(T)
    % Efficiently extract unique genes from table T.

    Ttmp = T(:, 2:end); % Exclude the first column
    geneCells = table2cell(Ttmp); %Convert to cell array.
    allGenes = [geneCells{:}];%Flatten.
    allGenes(allGenes == "") = [];%remove empty strings.
    uniqueGenes = unique(allGenes);% Find unique genes
end