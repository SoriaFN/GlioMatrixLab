%DIRECTIONALITY_INDEX
%Federico N. Soria and Mario Fernandez Ballester (2024) 

% Prompt user to select an Excel file
[filename, filepath] = uigetfile({'*.xlsx', 'Excel Files (*.xlsx)'}, 'Select the Excel file');
if isequal(filename, 0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', fullfile(filepath, filename)]);
end

% Read data from the selected Excel file
data = readtable(fullfile(filepath, filename));

% Extract unique track numbers
trackNumbers = unique(data.Track);

% Create a figure for plotting tracks
track_plot_fig = figure;
hold on;
axis equal; % Ensure equal scaling for X and Y axes

% Initialize arrays to store major and minor axes lengths, aspect ratios, track lengths, areas, directionality indexes, and territory indexes
major_axes = zeros(length(trackNumbers), 1);
minor_axes = zeros(length(trackNumbers), 1);
aspect_ratios = zeros(length(trackNumbers), 1);
track_lengths = zeros(length(trackNumbers), 1);
areas_ellipse = zeros(length(trackNumbers), 1);
directionality_indexes= zeros(length(trackNumbers), 1);

% Loop through each track number
for i = 1:length(trackNumbers)
    % Get the current track number
    currentTrack = trackNumbers(i);
    
    % Extract x and y positions for the current track
    x = data.X(data.Track == currentTrack);
    y = data.Y(data.Track == currentTrack);
    
    % Plot the track
    plot(x, y, '-o', 'LineWidth', 1);
    
    % Fit a convex hull to the track points
    k = convhull(x, y);
    poly_x = x(k);
    poly_y = y(k);
    
    % Plot the convex hull
    plot(poly_x, poly_y, '-.', 'LineWidth', 1);
    
    % Calculate major and minor axes of the convex hull
    centroid = mean([poly_x, poly_y]);
    poly_x_centered = poly_x - centroid(1);
    poly_y_centered = poly_y - centroid(2);
    
    % Compute the covariance matrix
    covariance_matrix = [poly_x_centered, poly_y_centered]' * [poly_x_centered, poly_y_centered] / length(poly_x_centered);
    
    % Eigenvectors and eigenvalues of the covariance matrix
    [V, D] = eig(covariance_matrix);
    
    % Sort eigenvalues in descending order
    [~, idx] = sort(diag(D), 'descend');
    
    % Major and minor axes lengths
    major_axis = 2 * sqrt(D(idx(1), idx(1)));
    minor_axis = 2 * sqrt(D(idx(2), idx(2)));
    
    % Store major and minor axis lengths
    major_axes(i) = major_axis;
    minor_axes(i) = minor_axis;
    
    % Compute aspect ratio
    aspect_ratio = major_axis / minor_axis;
    aspect_ratios(i) = aspect_ratio;
    
    % Calculate track length
    track_length = sum(sqrt(diff(x).^2 + diff(y).^2));
    track_lengths(i) = track_length;
       
    % Compute area of the ellipse
    area_ellipse = pi * (major_axis/2) * (minor_axis/2);
    areas_ellipse(i) = area_ellipse;
    
    % Compute Directionality Index
    directionality_index = aspect_ratio * (area_ellipse/track_length);
    directionality_indexes(i) = directionality_index;

    % Display track number next to the track and polygon
    text(mean(x), mean(y), num2str(currentTrack), 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% Add titles, labels, and grid for the track plot
title('Track Plot with Convex Hulls');
xlabel('X Position');
ylabel('Y Position');
grid off;
hold off;

%%
% Plot Territory covered vs Directionality Index in red
plotRegressionScatter(areas_ellipse, directionality_indexes, ...
                      'Territory covered (area)', 'Directionality Index', ...
                      'Correlation between Territory and Directionality', 'r');

% Plot Track Length vs Directionality Index in blue
plotRegressionScatter(track_lengths, directionality_indexes, ...
                      'Track Length', 'Directionality Index', ...
                      'Correlation between Track Length and Directionality Index', 'b');

%% 
% Construct the output filename based on the input filename
[~, name, ~] = fileparts(filename);
output_filename = fullfile(filepath, ['track_properties_', name, '.xlsx']);

% Export data to Excel
output_table = table(trackNumbers, major_axes, minor_axes, aspect_ratios, areas_ellipse, track_lengths, directionality_indexes);
writetable(output_table, output_filename);

% Save the track plot as .fig and .tif
savefig(track_plot_fig, fullfile(filepath, ['track_plot_', name, '.fig']));
print(track_plot_fig, fullfile(filepath, ['track_plot_', name, '.tif']), '-dtiff', '-r300');

disp(['Data exported to ', output_filename]);

% === Display Directionality Index Table as a Figure ===

% Sort by directionality index in descending order
[directionality_indexes_sorted, sort_idx] = sort(directionality_indexes, 'descend');
trackNumbers_sorted = trackNumbers(sort_idx);

% Create a sorted table to display
directionality_table = table(trackNumbers_sorted, directionality_indexes_sorted, ...
    'VariableNames', {'TrackNumber', 'DirectionalityIndex'});

% Create a UI figure window
dir_table_fig = uifigure('Name', 'Directionality Index Table');

% Create a UI table inside the figure
uit = uitable(dir_table_fig, ...
    'Data', directionality_table, ...
    'Position', [20 20 360 500]);

% Auto-resize columns
uit.ColumnWidth = 'auto';



function plotRegressionScatter(xData, yData, xLabel, yLabel, plotTitle, lineColor)
    figure;
    scatter(xData, yData, 'filled');
    title(plotTitle);
    xlabel(xLabel);
    ylabel(yLabel);
    grid on;
    hold on;

    % Fit linear regression
    p = polyfit(xData, yData, 1);
    x_fit = linspace(min(xData), max(xData), 100);
    y_fit = polyval(p, x_fit);

    % Plot regression line
    plot(x_fit, y_fit, [lineColor '--'], 'LineWidth', 2);

    % Compute Pearson correlation coefficient
    R = corrcoef(xData, yData);
    r_value = R(1, 2);
    r_squared = r_value^2;

    % Display regression equation and R value
    eq_str = sprintf('y = %.2f x + %.2f\nR = %.2f, R^2 = %.2f', p(1), p(2), r_value, r_squared);
    text(mean(x_fit), max(y_fit)*0.95, eq_str, 'FontSize', 12, 'Color', lineColor, ...
         'BackgroundColor', 'w', 'EdgeColor', lineColor, 'Margin', 2);
    hold off;
end