clear; clc;

% Step 1: Upload .dat file
[filename, pathname] = uigetfile('*.dat', 'Select a .dat file');
if isequal(filename, 0)
    disp('User cancelled file selection.');
    return;
end

fullFilePath = fullfile(pathname, filename);
disp(['Selected file: ', fullFilePath]);

% Step 2: Try reading the file using readmatrix, then fallback
try
    data = readmatrix(fullFilePath);
catch
    data = [];
end

% If readmatrix failed or data is empty, try textscan
if isempty(data)
    fid = fopen(fullFilePath, 'r');
    rawText = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = rawText{1};

    % Split by space or tab
    data = cellfun(@(x) strsplit(x), lines, 'UniformOutput', false);
    data = vertcat(data{:});
end

% Show preview
disp('Preview of data:');
disp(data(1:min(5,end), :));  % Show first 5 rows

% Step 3: Ask where to save the Excel file
[saveFileName, savePath] = uiputfile('*.xlsx', 'Save Excel file as');
if isequal(saveFileName, 0)
    disp('User cancelled saving.');
    return;
end

excelFilePath = fullfile(savePath, saveFileName);

% Step 4: Write to Excel
try
    if iscell(data)
        writecell(data, excelFilePath);
    else
        writematrix(data, excelFilePath);
    end
    disp(['Data successfully saved to ', excelFilePath]);
catch ME
    disp('Failed to write Excel file.');
    disp(ME.message);
end