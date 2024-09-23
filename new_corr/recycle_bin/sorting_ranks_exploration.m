x = [1 4 7 8 2]';
y = [2 4 9 8 10]';
iid = 1:length(y);
iid = iid';
tbl = table(iid, x, y)

[~, x_rank] = sortrows(tbl, 2, 'ascend')
[~, y_rank] = sortrows(tbl, 3, 'ascend')

tbl_ranked = table(iid, x, x_rank, y, y_rank);

tbl(idx,:)