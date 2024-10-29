% Read the table with column headers preserved
T = readtable("stability_analysisDV_final.csv", 'VariableNamingRule', 'preserve');

% Ensure the 'Type' column is treated as a categorical variable
T.Type = categorical(T.Type);

% Handle NaN values based on user selection
method = 1;  % You can modify this to test different NaN handling methods

% Handle NaN values based on selected method
switch method
    case 1
        T{:, 2:end} = fillmissing(T{:, 2:end}, 'constant', 0);       
    case 2
        T{:, 2:end} = fillmissing(T{:, 2:end}, 'mean');       
    case 3
        T{:, 2:end} = fillmissing(T{:, 2:end}, 'median');        
    case 4
        T{:, 2:end} = fillmissing(T{:, 2:end}, 'linear');        
    case 5
        T{:, 2:end} = fillmissing(T{:, 2:end}, 'previous');        
    case 6
        T{:, 2:end} = fillmissing(T{:, 2:end}, 'next');       
    otherwise
        disp('Invalid selection for NaN handling method.');
end

% Separate good and bad types
Tgood = T(T.Type == 'good', :);
Tgood_mat = Tgood{:, 2:end};
groups_good = categorical(Tgood.Properties.VariableNames(2:end));

Tbad = T(T.Type == 'bad', :);
Tbad_mat = Tbad{:, 2:end};
groups_bad = categorical(Tbad.Properties.VariableNames(2:end));

% Create a figure for the violin plots
figure;
% Plot the "Good" group with a specific color (e.g., green)
v1 = violinplot(groups_good, Tgood_mat);
for i = 1:length(v1)
    v1(i).FaceColor  = [0.49, 0.18, 0.56];  % Green color
    v1(i).EdgeColor  = [0.49, 0.18, 0.56];  % Green color
end
hold on;

% Plot the "Bad" group with a different color (e.g., red)
v2 = violinplot(groups_bad, Tbad_mat);
for i = 1:length(v2)
    v2(i).FaceColor  = [0.30, 0.75, 0.93];  % Red color
    v2(i).EdgeColor  = [0.30, 0.75, 0.93];  % Red color
end

% Add a legend indicating which color corresponds to "Good" and "Bad"
legend([v1(1), v2(1)], {'Good', 'Bad'}, 'Location', 'northeast');

% Title and labels
title('Violin Plots of Gene FoldChange(DiffDist) - Good vs Bad');
ylabel('FoldChange Values');
xtickangle(45);  % Rotate x-axis labels for better readability
hold off;