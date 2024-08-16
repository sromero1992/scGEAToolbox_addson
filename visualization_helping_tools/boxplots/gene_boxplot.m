function gene_boxplot(sce, my_gene)
    X = sce.X;
    X = sc_norm(X,'type','libsize');
    %X = log( X + 1);
    X = full(X);
    g = sce.g;
    nc = size(X,2);
    % Define your custom labels for x-axis
    my_labels = unique(sce.c_batch_id);
    nlabels = length(my_labels);

    ngenes = length(my_gene);
    
    % Statistics variables allocation
    min_val = zeros(nlabels,1);
    max_val = zeros(nlabels,1);
    median_val = zeros(nlabels,1);
    q1 = zeros(nlabels,1);
    q3 = zeros(nlabels,1);
    iqr = zeros(nlabels,1);
    for k = 1:ngenes
        idx = find( g == my_gene(k));
        %Xg = log( X(idx,:) + 1);
        fprintf("Processing gene %s \n",my_gene(k));
        if isempty(idx)
            continue
        end
        Xg = cell(1, nlabels);  % Pre-allocate cell array
        h = figure('Position',[0,0,1280,1024],"Visible","off");
        for i = 1:nlabels
            % Calculate statistics
            jdx = my_labels(i) == sce.c_batch_id;
            Xg{i} = X( idx, jdx);
            %fprintf("Xg size : %d %d \n", size(Xg));
            min_val(i) = min(Xg{i});
            max_val(i) = max(Xg{i});
            median_val(i) = median(Xg{i});
            q1(i) = prctile(Xg{i}, 25); % 25th percentile
            q3(i) = prctile(Xg{i}, 75); % 75th percentile
            iqr(i) = q3(i) - q1(i); % Interquartile range
        end
        for i = 1:nlabels
            % Initialize plot features
            subplot(1, nlabels, i);  % Divide figure into two subplots
            boxplot(Xg{i}, my_labels(i));
            ylim([min(min_val) max(max_val)])
            hold on;
            my_title = "Gene statistics ";
            my_title = {my_title, my_gene(k), my_labels(i)};
            title(my_title);
            text(1+0.02, max_val(i) - 0.2, ['Max: ', num2str(max_val(i))]);
            text(1+0.02, min_val(i) + 0.05, ['Min: ', num2str(min_val(i))]);
            text(1+0.02, median_val(i) + 0.3, ['Median: ', num2str(median_val(i))]);
            text(1-0.15, q3(i) - 0.05, ['Q3: ', num2str(q3(i))]);
            text(1-0.15, q1(i) + 0.05, ['Q1: ', num2str(q1(i))]);
        end
        % Create boxplot
        hold off;  % To add text on top of the plot
        saveas(h, sprintf('boxplot_%s.jpg', my_gene(k)))
    end
end