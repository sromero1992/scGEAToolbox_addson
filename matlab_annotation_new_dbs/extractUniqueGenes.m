function myUniqueGenes =  extractUniqueGenes(T)
    % Pre-weight stage
    ndb = size(T, 1);
    ngc = size(T, 2);
    % Find all unique genes first
    genes = strings(ngc, 1);
    iends = zeros(ndb, 1);
    fprintf("DB Size : %d  %d \n", size(T));

    % Corrected logic: Compare all rows with each other
    for idb = 1:ndb
        iends(idb) = sum(table2array(T(idb, :)) ~= "");
        for jdb = 1:ndb
            if idb ~= jdb % Avoid comparing a row with itself
                iends(jdb) = sum(table2array(T(jdb, :)) ~= "");
                tmp1 = table2array(T(idb, 1:iends(idb)));
                tmp2 = table2array(T(jdb, 1:iends(jdb)));
                tmp = union(tmp1, tmp2);
                genes = union(tmp, genes);
            end
        end
    end

    clear tmp1 tmp2 tmp;

    % Remove voids
    genes(genes == "") = [];
    myUniqueGenes = genes; % return the genes.
end