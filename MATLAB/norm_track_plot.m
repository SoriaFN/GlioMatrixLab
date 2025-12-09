% === PLOT TRACKS WITH RELATIVE COORDINATES AND REFERENCE CIRCLES ===

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
axis equal;
grid on;

% Loop through each track number
for i = 1:length(trackNumbers)
    currentTrack = trackNumbers(i);
    
    % Extract original x and y positions
    x_original = data.X(data.Track == currentTrack);
    y_original = data.Y(data.Track == currentTrack);
    
    % Recalculate positions relative to the first point
    x = x_original - x_original(1);
    y = y_original - y_original(1);
    
    % Plot the track (line only, linewidth 2)
    plot(x, y, '-', 'LineWidth', 2);
    
    % Mark origin of each track with a dot
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
end

% Fixed axis limits
axis_limit = 100;
xlim([-axis_limit, axis_limit]);
ylim([-axis_limit, axis_limit]);

% Draw concentric dotted circles every 50 units
radii = 10:10:axis_limit;
theta = linspace(0, 2*pi, 300);
for r = radii
    x_circ = r * cos(theta);
    y_circ = r * sin(theta);
    plot(x_circ, y_circ, 'k:', 'LineWidth', 1); % dotted circle
end

% Labels, grid, title
title('Track Plot (Centered with Reference Circles)');
xlabel('X Displacement');
ylabel('Y Displacement');

hold off;

% Save figure
[~, name, ~] = fileparts(filename);
savefig(track_plot_fig, fullfile(filepath, ['track_plot_relative_', name, '.fig']));
print(track_plot_fig, fullfile(filepath, ['track_plot_relative_', name, '.tif']), '-dtiff', '-r300');
