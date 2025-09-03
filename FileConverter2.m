% Clear previous data
clear; clc;

% --- Step 1: Select Excel (.xlsx) file ---
[filename, pathname] = uigetfile('*.xlsx', 'Select an Excel file to convert');
if isequal(filename, 0)
    disp('User cancelled file selection.');
    return;
end

excelFilePath = fullfile(pathname, filename);
disp(['Selected file: ', excelFilePath]);

% --- Step 2: Read Excel file ---
try
    data = readcell(excelFilePath);  % Works with both text and numeric data
catch ME
    error(['Error reading Excel file: ', ME.message]);
end

% --- Step 3: Prompt user to choose save location for .dat file ---
[saveFileName, savePath] = uiputfile('*.dat', 'Save as .dat file');
if isequal(saveFileName, 0)
    disp('User cancelled save.');
    return;
end

datFilePath = fullfile(savePath, saveFileName);

% --- Step 4: Write data to .dat file (tab-delimited) ---
try
    fid = fopen(datFilePath, 'w');
    
    for i = 1:size(data, 1)
        for j = 1:size(data, 2)
            entry = data{i, j};
            % Convert to string
            if isnumeric(entry)
                fprintf(fid, '%g', entry);  % For numbers
            elseif ischar(entry) || isstring(entry)
                fprintf(fid, '%s', entry);  % For text
            elseif isempty(entry)
                fprintf(fid, '');  % Leave blank for empty cells
            else
                fprintf(fid, '%s', string(entry));  % Fallback
            end
            
            % Add tab between columns except the last one
            if j < size(data, 2)
                fprintf(fid, '\t');
            end
        end
        fprintf(fid, '\n');  % Newline after each row
    end
    
    fclose(fid);
    disp(['Data successfully saved to ', datFilePath]);

catch ME
    if fid > 0
        fclose(fid);
    end
    error(['Failed to write .dat file: ', ME.message]);
end
