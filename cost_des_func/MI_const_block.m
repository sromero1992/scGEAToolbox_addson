function MI_mat = MI_const_block(data, data2, mode)
    % MI_construction computes the mutual information from data
    % INPUT:
    % data =====> Contains X count matrix and y target as follows
    %             data = [X; y]; (gene by (cell basis + target) )
    % mode =====> mode 1 is upper triangular, mode 2 is full
    % USAGE:
    % R0 = MI_construction(data);
    % save('R0.mat', 'R0', '-v7.3');

    if nargin < 2 || isempty(data2) ; data2 = data; end
    if nargin < 3 || isempty(mode) ; mode = 1; end

    % Data is count matrix with genes in rows and cells in columns
    data = sparse(data);
    data2 = sparse(data2);

    % transpose for efficient access (vectorization) across observations
    data = transpose(data);
    data2 = transpose(data2);
    nobs = size(data,1);
    nobs2 = size(data2, 1);
    if nobs ~= nobs2
        error("Not same number of observations in input matrices");
    end
    ngene = size(data,2);
    ngene2 = size(data2, 2);
    MI_mat = zeros(ngene, ngene2);
    % MI across data's rows (Computing upper triangular) in parallel 
    if mode == 1
        parfor jg = 1:ngene
            tmp = zeros(ngene,1);
            for ig = jg:ngene2
                % Computing pair MI with binning distribution
                tmp(ig) = BinPairMI( full(data(1:nobs,jg)), ...
                                     full(data2(1:nobs,ig)) );
            end
            MI_mat(jg,:) = tmp;
        end
    else
        parfor jg = 1:ngene
            tmp = zeros(ngene,1);
            for ig = 1:ngene2
                % Computing pair MI with binning distribution
                tmp(ig) = BinPairMI( full(data(1:nobs,jg)), ...
                                     full(data2(1:nobs,ig)) );
            end
            MI_mat(jg,:) = tmp;
        end
    end
    % Copy upper triangular to lower triangular
    %MI_mat = MI_mat + triu(MI_mat, 1)';
end