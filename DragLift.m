clc; clear;
% ---- Step 1: Let user select the .xlsx file ----
[filename, pathname] = uigetfile('*.xlsx', 'Select the .xlsx file');

if isequal(filename, 0)
    error('No file selected. Exiting...');
end

% Full path to file
filepath = fullfile(pathname, filename);
data = readtable(filepath);

% ---- Step 3: Extract columns ----
cx_values = data.cx;
cy_values = data.cy;

% ---- Step 4: Selecting the range for cx and cy values ----
min_no = input('Enter the starting index (min_no): ');
max_no = input('Enter the ending index (max_no): ');

% Ensure valid indices
if min_no < 1 || max_no > height(data) || min_no > max_no
    error('Invalid range specified.');
end

% Slice the data based on the given range
cx_values = cx_values(min_no:max_no);
cy_values = cy_values(min_no:max_no);

% ---- Step 5: Calculate mean for cx (indexed loop) ----
cx_mean = mean(cx_values);

% ---- Step 6: Calculate RMS for cy (direct iteration loop) ----
cy_rms = rms(cy_values);

% ---- Step 7: Output results ----
fprintf('Mean of cx: %.6f\n', cx_mean);
fprintf('RMS of cy: %.6f\n', cy_rms);