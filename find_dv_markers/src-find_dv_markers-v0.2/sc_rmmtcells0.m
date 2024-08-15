function [X, keptidx] = sc_rmmtcells0(X, genelist, mtratio, mtgenenamepat, vebrose)
%Remove cells with high mtDNA ratio
if nargin < 3, mtratio = 0.1; end
if nargin < 4, mtgenenamepat = "mt-"; end
if nargin < 5, vebrose = false; end

assert(size(X, 1) == length(genelist), ...
    sprintf('%s: size(X,1) is not equal to length(g).', upper(mfilename)))
idx = startsWith(genelist, mtgenenamepat, 'IgnoreCase', true);
if sum(idx) > 0
    if vebrose
        fprintf('%d mt-genes found.\n', sum(idx));
    end
    lbsz = sum(X, 1);
    lbsz_mt = sum(X(idx, :), 1);
    f_mtreads = lbsz_mt ./ lbsz;
    keptidx = f_mtreads < mtratio;
    if sum(~keptidx) > 0
        X = X(:, keptidx);
        %X(:,~keptidx)=[];
        if vebrose
            fprintf('%d cells with mt-read ratio >=%f (or %f%%) are removed.\n', ...
                sum(~keptidx), mtratio, mtratio*100);
        end
    else
        if vebrose
            fprintf('No cells with mt-read ratio >=%f (or %f%%) are removed.\n', ...
                mtratio, mtratio*100);
        end
    end
else
    if vebrose
        fprintf('No mt-genes found.\n');
        fprintf('No cells with mt-read ratio >=%f (or %f%%) are removed.\n', ...
            mtratio, mtratio*100);
    end
    if nargout > 1
        keptidx = true(size(X, 2), 1);
    end
end
end