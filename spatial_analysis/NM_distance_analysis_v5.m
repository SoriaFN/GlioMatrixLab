%% NM_distance_analysis.m  (v4)
%  Neuromelanin (NM) particle detection and distance-to-SN analysis
%  Federico Soria / Gliomatrix Lab - Achucarro Basque Center for Neuroscience
%
%  WHAT IT DOES:
%  1. Loads an RGB image of a midbrain section with NM pigment
%  2. User draws the SN boundary (polygon ROI)
%  3. User draws the outer analysis boundary (total NM-containing region)
%  4. Detects brown NM particles using HSV color thresholding (rejects
%     gray/black artifacts automatically)
%  5. Computes shortest distance from each particle centroid to the SN
%     boundary (negative = inside SN)
%  6. Exports results to an Excel (.xlsx) file and generates overlay figures
%
%  REQUIRES: Image Processing Toolbox
%
%  CHANGES vs v3:
%  - 3-frame TIFFs now ask for confirmation before being treated as RGB
%    (prevents silently merging a z-stack / time series into a fake RGB).
%  - Guarded XResolution access so missing-tag images fall back cleanly.
%  - Area filtering uses bwareafilt (one call instead of bwareaopen + loop).
%  - Distance-to-SN now uses bwdist's nearest-pixel map instead of an
%    O(particles x boundary) brute-force loop. Result is sub-pixel and
%    mathematically equivalent to v3, just much faster on large images.
%  - Inside/outside classification vectorised (no per-particle loop).
%  - Graceful handling of the zero-particles case.
%  - Colormap indexing in the distance figure guarded against /0.
%  - Removed an unused mask variable. No detection thresholds changed.
%
%  -----------------------------------------------------------------------
%  USER SETTINGS — adjust these before running
%  -----------------------------------------------------------------------

% Default pixel calibration (microns per pixel).
% The script will try to read this from the TIFF metadata (ImageJ format).
% If found, it will show you the value and ask you to confirm.
% If not found, it will use this default and ask you to confirm/change.
default_um_per_pixel = 0.409;

% HSV thresholds for brown NM detection
% Hue: 0-1 scale in MATLAB (0.028-0.139 corresponds to ~10-50 degrees)
% These defaults work well for unstained NM on pale background.
% If you're missing particles, try lowering hue_min or sat_min.
% If you're picking up too much background, raise sat_min or lower val_max.
hue_min = 0.022;   % ~10 degrees — lower bound of brown
hue_max = 0.180;   % ~50 degrees — upper bound of brown
sat_min = 0.04;    % minimum saturation (rejects gray/black artifacts)
val_max = 0.85;    % maximum brightness (rejects pale background)

% Median filter radius for cleaning (set to 0 to skip)
median_radius = 3;

% Minimum particle area in pixels (removes tiny specks)
min_area_px = 10;

% Maximum particle area in pixels (removes large blobs that aren't single
% NM granules — set to Inf to disable)
max_area_px = 500;

%  -----------------------------------------------------------------------
%  END OF USER SETTINGS — no need to edit below this line
%  -----------------------------------------------------------------------

%% 1. Load image
[filename, filepath] = uigetfile({'*.tif;*.tiff', 'TIFF images'; ...
    '*.jpg;*.jpeg;*.png', 'Other images'}, 'Select RGB image');
if isequal(filename, 0)
    error('No file selected.');
end
fullpath = fullfile(filepath, filename);

% Read image — handle both true RGB and ImageJ 3-frame composite TIFFs
info = imfinfo(fullpath);
if numel(info) == 3 && all([info.SamplesPerPixel] == 1) && info(1).BitDepth <= 16
    % Looks like an ImageJ composite (3 single-channel frames). This is
    % ambiguous — a z-stack or time series would look identical — so confirm.
    choice = questdlg(sprintf(['This TIFF has 3 single-channel frames.\n' ...
        'Treat them as the R, G and B channels of one RGB image?\n\n' ...
        'Choose "No" if this is actually a z-stack or time series.']), ...
        '3-frame TIFF detected', 'Yes, merge to RGB', 'No, cancel', ...
        'Yes, merge to RGB');
    if ~strcmp(choice, 'Yes, merge to RGB')
        error('Aborted: 3-frame TIFF was not confirmed as RGB channels.');
    end
    fprintf('Merging 3-frame composite TIFF to RGB...\n');
    R = imread(fullpath, 1);
    G = imread(fullpath, 2);
    B = imread(fullpath, 3);
    img = cat(3, R, G, B);
elseif numel(info) == 1
    img = imread(fullpath);
else
    error('Unexpected TIFF structure: %d frames. Expected 1 (RGB) or 3 (composite).', numel(info));
end

% Verify it's RGB
if size(img, 3) ~= 3
    error('Image must be RGB. Got %d channels.', size(img, 3));
end

% Convert to double [0,1] for processing
img_d = im2double(img);

%% 1b. Read calibration from TIFF metadata (ImageJ format)
um_per_pixel = default_um_per_pixel;
cal_source = 'default';

try
    tif_info = info(1);
    % ImageJ stores calibration as XResolution = pixels per unit,
    % with the unit specified in ImageDescription.
    if isfield(tif_info, 'ImageDescription') && ...
       isfield(tif_info, 'XResolution') && ...
       contains(tif_info.ImageDescription, 'unit=micron', 'IgnoreCase', true)
        xres = tif_info.XResolution;        % pixels per micron
        if ~isempty(xres) && xres > 0
            um_per_pixel = 1.0 / xres;
            cal_source = 'TIFF metadata (ImageJ)';
        end
    end
catch
    % If metadata reading fails, we'll just use the default
end

% Show calibration and ask user to confirm or change
if strcmp(cal_source, 'default')
    prompt_msg = sprintf('No calibration found in TIFF metadata.\nDefault: %.4f um/pixel.\n\nEnter calibration (um/pixel), or press Enter to accept default:', ...
        um_per_pixel);
else
    prompt_msg = sprintf('Calibration read from %s: %.4f um/pixel.\n\nPress Enter to accept, or type a new value:', ...
        cal_source, um_per_pixel);
end

user_cal = inputdlg(prompt_msg, 'Pixel calibration', 1, {num2str(um_per_pixel, '%.4f')});
if ~isempty(user_cal) && ~isempty(user_cal{1})
    new_val = str2double(user_cal{1});
    if ~isnan(new_val) && new_val > 0
        um_per_pixel = new_val;
    end
end
fprintf('Using calibration: %.4f um/pixel\n', um_per_pixel);

fprintf('Loaded: %s (%d x %d pixels, %.1f x %.1f um)\n', ...
    filename, size(img,2), size(img,1), ...
    size(img,2)*um_per_pixel, size(img,1)*um_per_pixel);

%% 2. Draw ROIs interactively
figure('Name', 'Draw ROIs', 'NumberTitle', 'off');
imshow(img);
title('STEP 1: Draw SN boundary (polygon). Double-click to close.', ...
    'FontSize', 14);

roi_sn = drawpolygon('Color', 'cyan', 'LineWidth', 2);
sn_vertices = roi_sn.Position;  % Nx2 [x, y]
sn_mask = createMask(roi_sn);
fprintf('SN ROI: %d vertices, %d pixels\n', size(sn_vertices,1), sum(sn_mask(:)));

title('STEP 2: Draw OUTER analysis boundary. Double-click to close.', ...
    'FontSize', 14);
roi_outer = drawpolygon('Color', 'yellow', 'LineWidth', 2);
outer_mask = createMask(roi_outer);
fprintf('Outer ROI: %d pixels\n', sum(outer_mask(:)));

%% 3. Color thresholding in HSV
hsv = rgb2hsv(img_d);
H = hsv(:,:,1);
S = hsv(:,:,2);
V = hsv(:,:,3);

% Brown detection mask
brown_mask = (H >= hue_min) & (H <= hue_max) & ...
             (S >= sat_min) & ...
             (V <= val_max);

% Apply median filter to remove salt-and-pepper noise
if median_radius > 0
    brown_mask = medfilt2(uint8(brown_mask), [median_radius median_radius]) > 0;
end

% Restrict detection to the outer analysis ROI only
brown_mask = brown_mask & outer_mask;

% Size filtering: keep components within [min_area_px, max_area_px].
% bwareafilt does min and max in one pass (replaces bwareaopen + manual loop).
if isfinite(max_area_px)
    if any(brown_mask(:))
        brown_mask = bwareafilt(brown_mask, [min_area_px, max_area_px]);
    end
else
    brown_mask = bwareaopen(brown_mask, min_area_px);
end

%% 4. Measure particles
cc = bwconncomp(brown_mask);
props = regionprops(cc, 'Centroid', 'Area', 'EquivDiameter', ...
    'MajorAxisLength', 'MinorAxisLength', 'Perimeter');

n_particles = cc.NumObjects;
fprintf('\nDetected %d NM particles\n', n_particles);

[~, basename, ~] = fileparts(filename);
outdir = filepath;  % save in same folder as input image
results_file = fullfile(outdir, [basename '_NM_results.xlsx']);

% Bail out gracefully if nothing was detected
if n_particles == 0
    warning('No NM particles detected. Check thresholds or ROIs. Saving empty results table.');
    T = cell2table(cell(0,11), 'VariableNames', {'ParticleID', ...
        'Centroid_X_px', 'Centroid_Y_px', 'Area_px', 'Area_um2', ...
        'EquivDiameter_um', 'Perimeter_um', 'Location', ...
        'Dist_to_SN_border_px', 'Dist_to_SN_border_um', 'Signed_dist_um'});
    writetable(T, results_file);
    fprintf('Empty results saved to: %s\nDone (no particles).\n', results_file);
    return;
end

% Extract per-particle measurements
centroids = vertcat(props.Centroid);  % Nx2 [x, y]
areas = [props.Area]';
equiv_diameters = [props.EquivDiameter]';
perimeters = [props.Perimeter]';

%% 5. Classify inside/outside SN and compute distances
% Centroid pixel coordinates (rounded, clamped to image bounds)
cx_px = min(max(round(centroids(:,1)), 1), size(img, 2));
cy_px = min(max(round(centroids(:,2)), 1), size(img, 1));
lin   = sub2ind(size(sn_mask), cy_px, cx_px);

% Inside-SN test: is the centroid pixel within the SN mask?
is_inside_sn = sn_mask(lin);

n_inside  = sum(is_inside_sn);
n_outside = sum(~is_inside_sn);
fprintf('  Inside SN:  %d particles\n', n_inside);
fprintf('  Outside SN: %d particles\n', n_outside);

% Shortest distance from each centroid to the SN boundary.
% bwdist returns, for every pixel, the distance to and index of the nearest
% boundary pixel. We look up the nearest boundary pixel for each centroid's
% pixel, then take the exact (sub-pixel) Euclidean distance to it. This
% reproduces the v3 brute-force result but in O(image) rather than
% O(particles x boundary).
sn_boundary = bwperim(sn_mask);
[~, nearestIdx] = bwdist(sn_boundary);
[nearest_boundary_y, nearest_boundary_x] = ...
    ind2sub(size(sn_mask), double(nearestIdx(lin)));

dist_to_sn_px = hypot(centroids(:,1) - nearest_boundary_x, ...
                      centroids(:,2) - nearest_boundary_y);

% Convert to microns
dist_to_sn_um      = dist_to_sn_px * um_per_pixel;
areas_um2          = areas * um_per_pixel^2;
equiv_diameters_um = equiv_diameters * um_per_pixel;
perimeters_um      = perimeters * um_per_pixel;

% Signed distance: negative for particles inside the SN (convention)
signed_dist_um = dist_to_sn_um;
signed_dist_um(is_inside_sn) = -signed_dist_um(is_inside_sn);

%% 6. Build results table
particle_id = (1:n_particles)';
location = cell(n_particles, 1);
location(is_inside_sn)  = {'inside_SN'};
location(~is_inside_sn) = {'outside_SN'};

T = table(particle_id, ...
    centroids(:,1), centroids(:,2), ...
    areas, areas_um2, ...
    equiv_diameters_um, perimeters_um, ...
    location, ...
    dist_to_sn_px, dist_to_sn_um, signed_dist_um, ...
    'VariableNames', {'ParticleID', 'Centroid_X_px', 'Centroid_Y_px', ...
    'Area_px', 'Area_um2', 'EquivDiameter_um', 'Perimeter_um', ...
    'Location', 'Dist_to_SN_border_px', 'Dist_to_SN_border_um', ...
    'Signed_dist_um'});

% Sort by distance (furthest outside first)
T = sortrows(T, 'Signed_dist_um', 'descend');

% Display summary
fprintf('\n--- SUMMARY ---\n');
fprintf('Total particles:     %d\n', n_particles);
fprintf('Inside SN:           %d\n', n_inside);
fprintf('Outside SN:          %d\n', n_outside);
if n_outside > 0
    outside_dists = dist_to_sn_um(~is_inside_sn);
    fprintf('Distance outside SN (um):\n');
    fprintf('  Mean:   %.1f\n', mean(outside_dists));
    fprintf('  Median: %.1f\n', median(outside_dists));
    fprintf('  Max:    %.1f\n', max(outside_dists));
    fprintf('  Min:    %.1f\n', min(outside_dists));
end
if n_inside > 0
    inside_areas = areas_um2(is_inside_sn);
    fprintf('Mean area inside SN:  %.1f um2\n', mean(inside_areas));
    if n_outside > 0
        outside_areas = areas_um2(~is_inside_sn);
        fprintf('Mean area outside SN: %.1f um2\n', mean(outside_areas));
    end
end

%% 7. Save results
writetable(T, results_file);
fprintf('\nResults saved to: %s\n', results_file);

%% 7b. Save summary to separate Excel file
summary_file = fullfile(outdir, [basename '_NM_summary.xlsx']);

summary_names  = {};
summary_values = [];

summary_names{end+1}  = 'Total_particles';
summary_values(end+1) = n_particles;
summary_names{end+1}  = 'Inside_SN';
summary_values(end+1) = n_inside;
summary_names{end+1}  = 'Outside_SN';
summary_values(end+1) = n_outside;

if n_outside > 0
    outside_dists = dist_to_sn_um(~is_inside_sn);
    summary_names{end+1}  = 'Dist_outside_mean_um';
    summary_values(end+1) = mean(outside_dists);
    summary_names{end+1}  = 'Dist_outside_median_um';
    summary_values(end+1) = median(outside_dists);
    summary_names{end+1}  = 'Dist_outside_max_um';
    summary_values(end+1) = max(outside_dists);
    summary_names{end+1}  = 'Dist_outside_min_um';
    summary_values(end+1) = min(outside_dists);
end

if n_inside > 0
    summary_names{end+1}  = 'Area_inside_mean_um2';
    summary_values(end+1) = mean(areas_um2(is_inside_sn));
end
if n_inside > 0 && n_outside > 0
    summary_names{end+1}  = 'Area_outside_mean_um2';
    summary_values(end+1) = mean(areas_um2(~is_inside_sn));
end

% Add calibration and thresholds for reproducibility
summary_names{end+1}  = 'um_per_pixel';
summary_values(end+1) = um_per_pixel;
summary_names{end+1}  = 'hue_min';
summary_values(end+1) = hue_min;
summary_names{end+1}  = 'hue_max';
summary_values(end+1) = hue_max;
summary_names{end+1}  = 'sat_min';
summary_values(end+1) = sat_min;
summary_names{end+1}  = 'val_max';
summary_values(end+1) = val_max;

T_summary = table(summary_names', summary_values', ...
    'VariableNames', {'Parameter', 'Value'});
writetable(T_summary, summary_file);
fprintf('Summary saved to: %s\n', summary_file);

%% 8. Generate figures

% --- Figure 1: Detection overlay ---
fig1 = figure('Name', 'NM Detection', 'NumberTitle', 'off', ...
    'Position', [50 50 900 900]);
imshow(img); hold on;

% Draw SN boundary
h_sn = plot(sn_vertices([1:end, 1], 1), sn_vertices([1:end, 1], 2), ...
    'c-', 'LineWidth', 2);

% Plot particles color-coded by location
h_inside = []; h_outside = [];
if n_inside > 0
    h_inside = scatter(centroids(is_inside_sn, 1), centroids(is_inside_sn, 2), ...
        30, 'g', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end
if n_outside > 0
    h_outside = scatter(centroids(~is_inside_sn, 1), centroids(~is_inside_sn, 2), ...
        30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

legend_handles = h_sn;
legend_labels = {'SN boundary'};
if n_inside > 0
    legend_handles(end+1) = h_inside;
    legend_labels{end+1} = 'Inside SN';
end
if n_outside > 0
    legend_handles(end+1) = h_outside;
    legend_labels{end+1} = 'Outside SN';
end
legend(legend_handles, legend_labels, 'FontSize', 12, 'Location', 'southwest', ...
    'TextColor', 'w', 'Color', [0.2 0.2 0.2]);
title(sprintf('NM particles: %d inside, %d outside SN', n_inside, n_outside), ...
    'FontSize', 14);
hold off;

saveas(fig1, fullfile(outdir, [basename '_NM_overlay.png']));

% --- Figure 2: Distance lines for outside particles ---
fig2 = figure('Name', 'Distances to SN', 'NumberTitle', 'off', ...
    'Position', [100 100 900 900]);
imshow(img); hold on;

% Draw SN boundary
plot(sn_vertices([1:end, 1], 1), sn_vertices([1:end, 1], 2), ...
    'c-', 'LineWidth', 2);

% Draw distance lines from each outside particle to nearest SN border
outside_idx = find(~is_inside_sn);
if ~isempty(outside_idx)
    % Color-code lines by distance
    dists_out = dist_to_sn_um(outside_idx);
    denom = max(max(dists_out), eps);   % guard against all-zero distances
    cmap = hot(256);

    for k = 1:length(outside_idx)
        idx = outside_idx(k);
        % Map distance to colormap
        ci = max(1, min(256, round(dists_out(k) / denom * 255) + 1));
        line_color = cmap(ci, :);

        plot([centroids(idx, 1), nearest_boundary_x(idx)], ...
             [centroids(idx, 2), nearest_boundary_y(idx)], ...
             '-', 'Color', [line_color 0.6], 'LineWidth', 1);
    end
    scatter(centroids(outside_idx, 1), centroids(outside_idx, 2), ...
        20, dists_out, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    colormap(hot);
    cb = colorbar;
    cb.Label.String = 'Distance to SN border (\mum)';
    cb.Label.FontSize = 12;
end
title('Distance of each NM particle to SN border', 'FontSize', 14);
hold off;

saveas(fig2, fullfile(outdir, [basename '_NM_distances.png']));

% --- Figure 3: Histogram of distances ---
fig3 = figure('Name', 'Distance Distribution', 'NumberTitle', 'off', ...
    'Position', [150 150 700 500]);

if n_outside > 0
    histogram(dist_to_sn_um(~is_inside_sn), 20, ...
        'FaceColor', [0.85 0.33 0.1], 'EdgeColor', 'w');
    xlabel('Distance to SN border (\mum)', 'FontSize', 13);
    ylabel('Number of NM particles', 'FontSize', 13);
    title('Distribution of NM particle distances from SN', 'FontSize', 14);

    % Add summary stats as text
    text(0.97, 0.93, sprintf('n = %d\nMedian = %.0f \\mum\nMean = %.0f \\mum', ...
        n_outside, median(outside_dists), mean(outside_dists)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'FontSize', 11, ...
        'BackgroundColor', 'w', 'EdgeColor', 'k');
end
set(gca, 'FontSize', 11);

saveas(fig3, fullfile(outdir, [basename '_NM_histogram.png']));

% --- Figure 4: Area vs distance scatter ---
fig4 = figure('Name', 'Area vs Distance', 'NumberTitle', 'off', ...
    'Position', [200 200 700 500]);

scatter(signed_dist_um, areas_um2, 40, ...
    'MarkerFaceColor', [0.85 0.33 0.1], 'MarkerEdgeColor', 'k', ...
    'MarkerFaceAlpha', 0.6);
xlabel('Signed distance to SN border (\mum)  [negative = inside]', 'FontSize', 13);
ylabel('Particle area (\mum^2)', 'FontSize', 13);
title('NM particle area vs. distance to SN', 'FontSize', 14);
xline(0, '--k', 'SN border', 'FontSize', 11, 'LabelVerticalAlignment', 'top');
set(gca, 'FontSize', 11);
grid on;

saveas(fig4, fullfile(outdir, [basename '_NM_area_vs_dist.png']));

fprintf('\nAll figures saved to: %s\n', outdir);
fprintf('Done!\n');
