% Read the table
t = readtable("AllMarker.xlsx");

% Define the primary and secondary sorting columns and their orders
% Primary sort: avg_log2FC (descending for highest fold changes)
% Secondary sort: p_val_adj (ascending for most significant p-values)
primary_sort_column = 'avg_log2FC';
primary_sort_order = 'descend'; % 'descend' for highest avg_log2FC

secondary_sort_column = 'p_val_adj';
secondary_sort_order = 'ascend'; % 'ascend' for lowest p_val_adj (most significant)

% Convert group identifiers to a numerical array for easier comparison
all_grps = double(string(t.grp));
unique_grps = unique(all_grps);

% Initialize an empty table to store the top N markers from all groups
t_top_markers_combined = table();

% Define the number of top markers to extract
num_top_markers = 20;

% Loop through each unique group
for igrp = unique_grps'
    % Find rows belonging to the current group
    idx_grp = igrp == all_grps;
    t_grp = t(idx_grp, :);

    % --- Sort the current group's data by the specified criteria ---
    % When sorting by multiple columns, provide them as a cell array.
    % The order of columns in the cell array determines the sort priority.
    % The order of sort directions in the cell array corresponds to the columns.
    t_grp_sorted = sortrows(t_grp, ...
                           {primary_sort_column, secondary_sort_column}, ...
                           {primary_sort_order, secondary_sort_order});

    % Extract the top N markers
    num_rows_in_grp = size(t_grp_sorted, 1);
    if num_rows_in_grp > 0
        % Determine how many markers to actually take (min of num_top_markers or available rows)
        rows_to_take = min(num_top_markers, num_rows_in_grp);
        t_top_n_current_grp = t_grp_sorted(1:rows_to_take, :);

        % Add a column to identify the original group for these markers
        t_top_n_current_grp.OriginalGroup = repmat(string(igrp), rows_to_take, 1);

        % Vertically concatenate the top N markers for the current group to the combined table
        t_top_markers_combined = [t_top_markers_combined; t_top_n_current_grp];
    else
        fprintf('Warning: Group %s has no data.\n', num2str(igrp));
    end
end

% Write the combined table of top markers to a single Excel file
writetable(t_top_markers_combined, "extracted_top10_markers_all_groups_sorted.xlsx");

disp(['Extraction complete. Top ', num2str(num_top_markers), ' markers per group, sorted by avg_log2FC and p_val_adj, saved to extracted_top20_markers_all_groups_sorted.xlsx']);