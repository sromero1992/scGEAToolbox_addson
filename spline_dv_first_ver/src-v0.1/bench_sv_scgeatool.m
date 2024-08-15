T=readtable('Anti-DE datasets2.csv','ReadVariableNames',false);
T.Properties.VariableNames = {'sample_id', 'species','tissue','geo_acc','geo_mat','batch_id'};

nsamples = length( unique(T.sample_id));
for ids = 1:nsamples
    idx = T.sample_id == ids;
    T1 = T(idx,:);
    gse = string( strtrim(T1.geo_acc(1)) );
    
    % Read all files in current path
    main_path = strcat(gse,"/");
    nsamples = size(T1,1);
    samples = string( strtrim(T1.batch_id) );

    fprintf("Processing sample_id %d \n",ids);

    % If the directory has finished successfully, skip it
    fname3 = strcat(main_path,"finished.txt");
    %if isfile(fname3)
    %    continue;
    %end

    % We can merge more than these and do the sce1 and sce2 comparison
    % It's up to tester 
    path1 = strcat(main_path,samples(1),".mat");
    path2 = strcat(main_path,samples(2),".mat");
    data  = load(path1);
    data2 = load(path2);
    sce1 = data.sce;
    sce2 = data2.sce;
    clear data data2;

    label0 = strcat( samples(1), "_", samples(2) ); 

    % Cell types of dataset1
    celltype1 = sce1.c_cell_type_tx;
    % Cell types of dataset1
    celltype2 = sce2.c_cell_type_tx;
    % Different cells along dataset1
    CellTypeList1 = unique(celltype1);
    % Different cells along dataset1
    CellTypeList2 = unique(celltype2);
    
    % Mapping the ideces of cell occurences 
    n = size(CellTypeList1,1);
    m = size(CellTypeList2,1);
    maxval = max(m,n);
    idx = zeros(maxval,2);
    k = 1;
    for i = 1:n
        for j = 1:m
            if CellTypeList1(i) == CellTypeList2(j)
                newStr = eraseBetween(CellTypeList1(i),"[","]");
                newStr = erase(newStr,"[");
                newStr = erase(newStr,"]");
                newStr = strrep(newStr,"/","-");
                newStr = strrep(newStr,"\","-");
                CellTypeList1(i) = newStr;
                CellTypeList2(j) = newStr;
                fprintf("Cell %s in %d , %d\n", CellTypeList1(i),i,j);
                idx(k,:) = [i, j] ;
                k = k + 1;
            end
        end
    end
    
    maxval1 = nnz(idx(:,1));
    maxval2 = nnz(idx(:,2));
    maxval = min(maxval1, maxval2);
    % Taking out zero occurences in mapped idx
    idx = idx(1:maxval,:);
    for k = 1:maxval
        fprintf("-----Processing %s %s %d  ... \n", CellTypeList1(idx(k,1)), k);
        i = idx(k,1);
        j = idx(k,2);
        % Getting cell subsample 
        cellidx  = celltype1 == CellTypeList1(i);
        cellidx2 = celltype2 == CellTypeList2(j);
        label1 = strcat(CellTypeList1(i),"-", label0);
        fname1 = strcat(main_path,label1);
    
        % Take out 5% mitochondrial cells
        [X, g ] = sc_selectg(sce1.X(:,cellidx), sce1.g, 1, 0.05);
        [X2,g2] = sc_selectg(sce2.X(:,cellidx2), sce2.g, 1, 0.05);
    
        % Work only with common genes
        [gl, irows, jrows] = intersect(g, g2, 'stable');
        X = X(irows,:);
        X2 = X2(jrows,:); 

        if isempty(gl)
            continue;
        end

        % Assing a boolean to normalize
        [X] = sc_norm(X,'type','libsize');
        [X2]= sc_norm(X2,'type','libsize');
    
        % Differential variability
        [Tdv, Tig, influgenes] = sc_variability(X, X2, gl, gl, fname1, false);
      
        % Differential expression
        [Tde, Tup, Tdn] = sc_deg(X, X2, gl, 1, false);
        label2 = strcat("DE-", label1);
        fname2 = strcat(main_path, label2);
        writetable(Tde, fname2, 'FileType', 'spreadsheet', 'Sheet', 'All genes');
        writetable(Tup, fname2, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
        writetable(Tdn, fname2, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');

    end
    % Writting this if finished the loop for checking up
    fname3 = strcat(main_path,"finished.txt");
    writelines("Finished successfully the pipeline...",fname3)
    type(fname3);
end


