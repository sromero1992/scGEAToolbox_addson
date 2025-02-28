fname = "AllSamplesMerged_HMOs6Removed_OldGenesRemoved_1500-20-15-500Filtered_DoubletsRemoved_AmbientRNARemoved_1500-20-15-500Filtered_AllGenesUMAP_Col_LeidenRes1_29Clusters_CellTypes_UseThisOne.mat";
% Load data
data = load(fname);
sce = data.sce;
clear data;

% Original batch IDs
batches_all = sce.c_batch_id; % Original batch IDs
batches = unique(batches_all); % Get unique batch IDs

% Group labels
gbatch = ["B. infantis #", "Control #", "HMOs #", "B. infantis + HMOs #"];
gbatch_cleaned = ["B_infantis", "Control", "HMOs", "B_infantis_+_HMOs"]; % Corresponding cleaned labels
batches_treat = batches_all;

% Rename batch IDs for clarity
for j = 1:length(gbatch)
    for i = 1:length(batches)
        % Check if the current batch contains the group label
        if contains(batches(i), gbatch(j))
            % Find indices where the batch matches
            batch_idx = batches_all == batches(i);

            % Assign cleaned group label to `batches_treat`
            batches_treat(batch_idx) = gbatch_cleaned(j);

        end
    end
end

batches_all = strrep(batches_all,' ','_');

% Convert to categorical variables
batches_treat = categorical(batches_treat, gbatch_cleaned, 'Ordinal', true);
batches_all = categorical(batches_all);

sce.c_batch_id = batches_treat;

fname = "simplified_batches_binf_hmo.mat";
save(fname,'sce',"-v7.3")