%% DETECT_MULTIEVENTS
% Scans a binary calcium event matrix (ROIs x time) using an expanding
% window to detect "multievents" — timepoints where 2 or more ROIs fire
% within a short time window. Starts with a 3s window; if a multievent is
% found, the window expands in 3s steps as long as new ROIs keep joining.
% The multievent level is stamped at the row of the first event in the
% window. Scanning resumes after the end of the expanded window.
%
% INPUT:  .xlsx file with ROIs as columns, rows as time (1 row = 1s),
%         cell values are 0 or 1.
% OUTPUT: Same .xlsx file with an added column "Multievent_level".
%
% Multievent levels:
%   0 = no multievent
%   2 = 2 ROIs active within the window
%   3 = 3 ROIs active, etc.
%
% Federico Soria lab — calcium imaging analysis

clear; clc;

%% --- 1. Load the Excel file ---

[fileName, filePath] = uigetfile({'*.xlsx;*.xls', 'Excel files'}, ...
    'Select the calcium event table');
if isequal(fileName, 0)
    error('No file selected. Aborting.');
end

fullPath = fullfile(filePath, fileName);
T = readtable(fullPath);

% Convert table to numeric matrix (rows = time, cols = ROIs)
data = table2array(T);
[nTimepoints, nROIs] = size(data);

fprintf('Loaded %s\n', fileName);
fprintf('  %d timepoints, %d ROIs\n', nTimepoints, nROIs);

%% --- 2. Parameters ---

baseWindow = 3;    % initial rolling window in seconds (= rows)
expandStep = 3;    % expansion step if new ROIs are found (seconds)
minROIs    = 2;    % minimum number of ROIs to qualify as multievent

%% --- 3. Scan with expanding window ---

% Preallocate output: one value per timepoint
multieventLevel = zeros(nTimepoints, 1);   % for the Excel column (stamped at first event)
heatmapLevel    = zeros(nTimepoints, 1);   % for the plot (painted across full span)

t = 1;  % current timepoint index

while t <= nTimepoints

    % --- Step A: Check the initial 3s window ---
    winEnd = min(t + baseWindow - 1, nTimepoints);
    window = data(t:winEnd, :);
    activeROIs = any(window, 1);  % logical: which ROIs fired in this window
    nActive = sum(activeROIs);

    if nActive >= minROIs
        % --- Step B: Expand window while new ROIs keep joining ---
        % Each expansion adds another 3s chunk. If any NEW ROI appears
        % in the new chunk, keep expanding. Stop when no new ROIs are found.
        expanding = true;
        while expanding
            nextEnd = min(winEnd + expandStep, nTimepoints);
            if nextEnd == winEnd
                break;  % reached end of data
            end
            % Check only the NEW chunk for additional ROIs
            newChunk = data((winEnd + 1):nextEnd, :);
            newActiveROIs = any(newChunk, 1);

            % Are there ROIs in the new chunk that weren't already active?
            addedROIs = newActiveROIs & ~activeROIs;
            if any(addedROIs)
                % Expand: merge the new ROIs and extend the window
                activeROIs = activeROIs | newActiveROIs;
                winEnd = nextEnd;
            else
                expanding = false;
            end
        end

        % Final count after all expansions
        nActive = sum(activeROIs);

        % --- Step C: Find first and last rows with events in this window ---
        % Only consider rows with activity in the ROIs that are part of
        % this multievent (i.e. the activeROIs).
        firstEventRow = t;
        lastEventRow  = t;
        for row = t:winEnd
            if any(data(row, activeROIs))
                firstEventRow = row;
                break;
            end
        end
        for row = winEnd:-1:t
            if any(data(row, activeROIs))
                lastEventRow = row;
                break;
            end
        end

        % Stamp the multievent level at the first event's row (for Excel)
        multieventLevel(firstEventRow) = nActive;

        % Paint the full span in the heatmap vector
        heatmapLevel(firstEventRow:lastEventRow) = nActive;

        % Resume scanning from the end of the expanded window
        t = winEnd + 1;
    else
        % No multievent — advance one second
        t = t + 1;
    end
end

%% --- 4. Report summary ---

nMultievents = sum(multieventLevel > 0);
fprintf('\nResults:\n');
fprintf('  Total multievents detected: %d\n', nMultievents);

% Break down by level
levels = unique(multieventLevel(multieventLevel > 0));
for i = 1:length(levels)
    count = sum(multieventLevel == levels(i));
    fprintf('  Level %d: %d events\n', levels(i), count);
end

%% --- 4b. Raster plot with multievent heatmap ---

figure('Name', 'Multievent Detection', 'Color', 'w', ...
    'Position', [100 100 800 800]);

% --- Top subplot: raster of individual ROI events (bigger dots) ---
ax1 = subplot(4, 1, 1:3);  % raster gets most of the vertical space
hold on;
tickHeight = 0.4;  % half-height of each vertical tick
for roi = 1:nROIs
    eventTimes = find(data(:, roi));
    if ~isempty(eventTimes)
        xCoords = [eventTimes'; eventTimes'; nan(1, length(eventTimes))];
        yCoords = [repmat(roi - tickHeight, 1, length(eventTimes)); ...
                   repmat(roi + tickHeight, 1, length(eventTimes)); ...
                   nan(1, length(eventTimes))];
        plot(xCoords(:), yCoords(:), '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 3);
    end
end
ylabel('ROI');
title('Calcium events raster');
set(ax1, 'YLim', [0 nROIs + 1], 'XLim', [1 nTimepoints]);
set(ax1, 'XTickLabel', []);  % hide x labels, shared with bottom
box on;
hold off;

% --- Bottom subplot: heatmap strip of multievent levels ---
ax2 = subplot(5, 1, 5);

% imagesc displays the multievent level as a 1-row color strip
% heatmapLevel spans the full duration of each multievent (first to last event)
imagesc(heatmapLevel');

% Viridis colormap
viridis = [0.267004 0.004874 0.329415; 0.282327 0.140926 0.457517; ...
           0.253935 0.265254 0.529983; 0.206756 0.371758 0.553117; ...
           0.163625 0.471133 0.558148; 0.127568 0.566949 0.550556; ...
           0.134692 0.658636 0.517649; 0.208030 0.745110 0.451882; ...
           0.360741 0.822786 0.358206; 0.553392 0.889434 0.226055; ...
           0.778816 0.942459 0.104616; 0.993248 0.906157 0.143936];
viridisMap = interp1(linspace(1,12,12), viridis, linspace(1,12,256));

colormap(ax2, viridisMap);
clim([0 8]);
cb = colorbar;
cb.Label.String = 'Multievent level';

% Format: no y-axis (it's just one row), x-axis = time
set(ax2, 'YTick', [], 'XLim', [0.5 nTimepoints + 0.5]);
set(ax2, 'TickLength', [0 0]);
xlabel('Time (s)');
title('Multievent level');
box on;

% --- Force both axes to the same horizontal extent ---
% The colorbar shrinks ax2, so we read ax1's position and apply it to ax2,
% then place the colorbar just outside the figure area.
drawnow;  % let MATLAB compute layout first
set([ax1, ax2], 'FontName', 'Arial', 'FontSize', 14)
set(ax1, 'LineWidth', 2);
set(ax2, 'LineWidth', 0.1);
pos1 = get(ax1, 'Position');  % [left bottom width height]
pos2 = get(ax2, 'Position');
set(ax2, 'Position', [pos1(1) pos2(2) pos1(3) pos2(4)]);

% Move colorbar to the right, outside the plot area
cbPos = get(cb, 'Position');
set(cb, 'Position', [pos1(1) + pos1(3) + 0.01, pos2(2), cbPos(3), pos2(4)]);

% Link x-axes so zooming/panning is synchronized
linkaxes([ax1, ax2], 'x');

% Append the multievent level as a new column to the original table
T.Multievent_level = multieventLevel;

% Overwrite the original Excel file with the added column
writetable(T, fullPath);
fprintf('\nColumn "Multievent_level" added to: %s\n', fileName);
