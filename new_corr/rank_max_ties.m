% function [ranks_forward, ranks_backward] = rank_min_ties_version2(Y)
%     % rank_min_ties_version2 computes ranks for a vector Y with ties handled using the 'max' method for forward and 'min' for backward.
%     % INPUT:
%     %   Y - Numeric vector to rank.
%     % OUTPUTS:
%     %   ranks_forward - Ranks with ties assigned the maximum rank in forward order.
%     %   ranks_backward - Ranks with ties assigned the minimum rank in backward order.
% 
%     % Compute initial ranks in ascending order
%     [~, sorted_indices] = sort(Y, 'ascend');  % Sort Y in ascending order to get indices
% 
%     % Initialize rank matrices
%     n = length(Y);
%     ranks_forward = zeros(n, 1);  % Allocate memory for forward ranks
% 
%     % Compute forward ranks using sorted indices
%     ranks_temp = 1:n;  % Temporary rank values in order of sorted Y
%     for i = 1:n
%         ranks_forward(sorted_indices(i)) = ranks_temp(i);  % Assign ranks back to original positions
%     end
% 
%     % Initialize backward ranks
%     ranks_backward = zeros(n, 1);
% 
%     % Handle ties: Assign the maximum rank for forward and minimum rank for backward
%     unique_values = unique(Y);
%     for i = 1:length(unique_values)
%         % Get logical indices of tied elements
%         tied_indices = (Y == unique_values(i));
%         if sum(tied_indices) > 1  % Only adjust if there are ties
%             % Get the maximum rank among the tied values for forward ranks
%             max_rank = max(ranks_forward(tied_indices));
%             ranks_forward(tied_indices) = max_rank;
% 
%             % Get the minimum rank among the tied values for backward ranks
%             min_rank = min(ranks_forward(tied_indices));
%             ranks_backward(tied_indices) = min_rank;
%         else
%             % Assign the original rank if no ties
%             ranks_forward(tied_indices) = ranks_temp(sorted_indices(tied_indices));
%             ranks_backward(tied_indices) = ranks_temp(sorted_indices(tied_indices));
%         end
%     end
% 
%     % Convert forward ranks to backward ranks
%     ranks_backward = max(ranks_backward) + 1 - ranks_backward;
% end

function [ranks_forward, ranks_backward] = rank_max_ties(Y)
    % rank_max_ties computes ranks for a vector Y with ties handled using the 'max' method.
    % INPUT:
    %   Y - Numeric vector to rank.
    % OUTPUTS:
    %   ranks_forward - Ranks with ties assigned the maximum rank in forward order.
    %   ranks_backward - Ranks with ties assigned the maximum rank in backward order.

    % Compute initial ranks in ascending order
    [~, sorted_indices] = sort(Y, 'ascend');  % Sort Y in ascending order to get indices

    % Initialize rank matrices
    n = length(Y);
    ranks_forward = zeros(n, 1);  % Allocate memory for forward ranks

    % Compute forward ranks using sorted indices
    ranks_temp = 1:n;  % Temporary rank values in order of sorted Y
    for i = 1:n
        ranks_forward(sorted_indices(i)) = ranks_temp(i);  % Assign ranks back to original positions
    end
    
    % Handle ties: Assign the maximum rank for all tied elements
    unique_values = unique(Y);
    for i = 1:length(unique_values)
        % Get logical indices of tied elements
        tied_indices = (Y == unique_values(i));
        if sum(tied_indices) > 1  % Only adjust if there are ties
            % Get the maximum rank among the tied values
            max_rank = max( ranks_forward(tied_indices) );
            % Assign the maximum rank to all tied values
            ranks_forward(tied_indices) = max_rank;
        end
    end

    % Compute backward ranks
    ranks_backward = n + 1 - ranks_forward;  % Convert forward ranks to backward ranks
end