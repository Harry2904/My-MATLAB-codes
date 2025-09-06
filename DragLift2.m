% ---- Step 1: Pick the file ----
[filename, pathname] = uigetfile({'*.dat;*.txt;*.csv','Data files (*.dat,*.txt,*.csv)'}, ...
                                 'Select the data file');
if isequal(filename,0), error('No file selected.'); end
filepath = fullfile(pathname, filename);

% ---- Step 2: Read header line and clean it ----
fid = fopen(filepath,'r');
assert(fid>0,'Cannot open file: %s', filepath);

firstLine = fgetl(fid);
if ~ischar(firstLine)
    fclose(fid); error('File is empty.');
end

% Remove UTF-8 BOM if present and the "variables =" prefix
firstLine = regexprep(firstLine, '^\xEF\xBB\xBF', '');                % BOM
firstLine = regexprep(firstLine, '^\s*variables\s*=\s*', '', 'ignorecase');

% Split header into variable names (commas and/or whitespace)
rawNames = regexp(firstLine, '[,\s]+', 'split');
rawNames = rawNames(~cellfun('isempty',rawNames));                     % drop empties

% ---- Step 3: Read numeric data (handles spaces/tabs/commas) ----
fmt = repmat('%f', 1, numel(rawNames));                                % one %f per column
C = textscan(fid, fmt, 'Delimiter', {' ','\t',','}, ...
                    'MultipleDelimsAsOne', true, 'CollectOutput', true);
fclose(fid);

if isempty(C) || isempty(C{1})
    error('No numeric data found. Check delimiters and number of columns.');
end

M = C{1};
if size(M,2) ~= numel(rawNames)
    % If the file has a different number of numeric columns than names, explain why
    error(['Column count mismatch: header has %d names (%s) but data rows have %d columns.\n' ...
           'This usually happens when the header uses commas but rows are space-delimited, ' ...
           'or there are extra tokens in the header.'],numel(rawNames), strjoin(rawNames, ', '), size(M,2));
end

% Make safe/unique MATLAB variable names but keep cx/cy intact
varNames = matlab.lang.makeValidName(rawNames, 'ReplacementStyle', 'delete');
varNames = matlab.lang.makeUniqueStrings(varNames);
data = array2table(M, 'VariableNames', varNames);

% ---- Step 4: Verify cx & cy exist and proceed ----
if ~ismember('cx', data.Properties.VariableNames) || ~ismember('cy', data.Properties.VariableNames)
    error('Columns "cx" and/or "cy" not found. Found: %s', strjoin(data.Properties.VariableNames, ', '));
end

% Drag and Lift coefficient values for plots
cx = data.cx;
cy = data.cy;

% Example: select range and compute stats
% min_no = input('Enter the starting index (min_no): ');
% max_no = input('Enter the ending index (max_no): ');
% if min_no < 1 || max_no > height(data) || min_no > max_no
%     error('Invalid range specified.');
% end
% cx_values = data.cx;
% cy_values = data.cy;
% 
% cx_mean = mean(cx_values);
% cy_rms  = rms(cy_values);
% 
% fprintf('Mean of cx: %.6f\n', cx_mean);
% fprintf('RMS of cy: %.6f\n',  cy_rms);

% ---- User inputs stable time range ----
t_start = input('Enter start time of stable region: ');
t_end   = input('Enter end time of stable region: ');

% Extract time column
time = data.time;

% Logical mask for time range
mask = (time >= t_start & time <= t_end);

% Filter cx and cy within that time range
cx_values = data.cx(mask);
cy_values = data.cy(mask);

% Calculate statistics
cx_mean = mean(cx_values);
cy_rms  = rms(cy_values);

fprintf('Stable region (%.2f - %.2f): Mean(Cx)=%.6f, RMS(Cy)=%.6f\n', ...
        t_start, t_end, cx_mean, cy_rms);

%Plotting
figure;
plot(time, cx, 'b-', 'LineWidth', 1.5); hold on;  % Plot Cd in blue
plot(time, cy, 'r-', 'LineWidth', 2.0);           % Plot Cl in red
hold off;
xlabel('Time');
ylabel('C_d & C_l');
title('Drag and Lift Coefficients vs Time');
legend('C_d (Drag Coefficient)', 'C_l (Lift Coefficient)');
grid on;