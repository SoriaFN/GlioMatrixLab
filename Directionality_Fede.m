%DIRECTIONALITY_INDEX
%Federico Soria (2024) 

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
figure;
hold on;

% Initialize arrays to store major and minor axes lengths, aspect ratios, track lengths, areas, directionality indexes, and territory indexes
major_axes = zeros(length(trackNumbers), 1);
minor_axes = zeros(length(trackNumbers), 1);
aspect_ratios = zeros(length(trackNumbers), 1);
track_lengths = zeros(length(trackNumbers), 1);
areas_ellipse = zeros(length(trackNumbers), 1);
directionality_indexes = zeros(length(trackNumbers), 1);
territory_indexes = zeros(length(trackNumbers), 1);

% Loop through each track number
for i = 1:length(trackNumbers)
    % Get the current track number
    currentTrack = trackNumbers(i);
    
    % Extract x and y positions for the current track
    x = data.X(data.Track == currentTrack);
    y = data.Y(data.Track == currentTrack);
    
    % Plot the track
    plot(x, y, '-o', 'LineWidth', 1.5);
    
    % Fit a convex hull to the track points
    k = convhull(x, y);
    poly_x = x(k);
    poly_y = y(k);
    
    % Plot the convex hull
    plot(poly_x, poly_y, '-.', 'LineWidth', 1.5);
    
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
    
    % Compute Directionality index
    directionality_index = aspect_ratio * track_length;
    directionality_indexes(i) = directionality_index;

    % Compute territory index
    territory_index = area_ellipse * track_length;
    territory_indexes(i) = territory_index;

    % Display track number next to the track and polygon
    text(mean(x), mean(y), num2str(currentTrack), 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% Add titles, labels, and grid for the track plot
title('Track Plot with Fitted Ellipses and Convex Hulls');
xlabel('X Position');
ylabel('Y Position');
grid on;
hold off;

% Plot the correlation between Territory Index and Directionality Index
figure;
scatter(territory_indexes, directionality_indexes, 'filled');
title('Correlation between TI and DI');
xlabel('Territory Index');
ylabel('Directionality Index');
grid on;

% Construct the output filename based on the input filename
[~, name, ~] = fileparts(filename);
output_filename = fullfile(filepath, ['track_properties_', name, '.xlsx']);

% Export data to Excel
output_table = table(trackNumbers, major_axes, minor_axes, aspect_ratios, areas_ellipse, track_lengths, directionality_indexes, territory_indexes);
writetable(output_table, output_filename);

disp(['Data exported to ', output_filename]);
