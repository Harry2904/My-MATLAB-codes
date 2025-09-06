% ---- Step 1: Ask user to select 3 files ----
numFiles = 3;
filenames = cell(numFiles,1);
pathnames = cell(numFiles,1);

for k = 1:numFiles
    [fname, pname] = uigetfile({'*.dat;*.txt;*.csv','Data files (*.dat,*.txt,*.csv)'}, ...
                               sprintf('Select data file %d of %d', k, numFiles));
    if isequal(fname,0)
        error('File %d not selected. Aborting.', k);
    end
    filenames{k} = fname;
    pathnames{k} = pname;
end

% ---- Step 2: Loop through files and load data ----
allData = cell(numFiles,1);

for k = 1:numFiles
    filepath = fullfile(pathnames{k}, filenames{k});
    fprintf('Reading file: %s\n', filepath);

    fid = fopen(filepath,'r');
    assert(fid>0,'Cannot open file: %s', filepath);

    % Read header line
    firstLine = fgetl(fid);
    if ~ischar(firstLine)
        fclose(fid); error('File is empty: %s', filepath);
    end

    % Clean header
    firstLine = regexprep(firstLine, '^\xEF\xBB\xBF', ''); % BOM
    firstLine = regexprep(firstLine, '^\s*variables\s*=\s*', '', 'ignorecase');
    rawNames = regexp(firstLine, '[,\s]+', 'split');
    rawNames = rawNames(~cellfun('isempty',rawNames));

    % Read numeric data
    fmt = repmat('%f', 1, numel(rawNames));
    C = textscan(fid, fmt, 'Delimiter', {' ','\t',','}, ...
                        'MultipleDelimsAsOne', true, 'CollectOutput', true);
    fclose(fid);

    if isempty(C) || isempty(C{1})
        error('No numeric data in file: %s', filepath);
    end

    M = C{1};
    if size(M,2) ~= numel(rawNames)
        error('Column mismatch in file: %s', filepath);
    end

    % Variable names
    varNames = matlab.lang.makeValidName(rawNames, 'ReplacementStyle', 'delete');
    varNames = matlab.lang.makeUniqueStrings(varNames);

    % Store as table
    data = array2table(M, 'VariableNames', varNames);

    % Verify required columns
    if ~ismember('cx', data.Properties.VariableNames) || ~ismember('cy', data.Properties.VariableNames)
        error('Columns "cx" and/or "cy" not found in file: %s', filepath);
    end

    allData{k} = data;
end

% Drag and Lift coefficient values for plots
cx = data.cx;
cy = data.cy;

% % ---- Step 3: Plot Cd (cx) and Cl (cy) for all 3 files ----
% figure; hold on;
% colors = lines(numFiles); % different colors for clarity
% 
% for k = 1:numFiles
%     data = allData{k};
% 
%     % If you have a 'time' column, use it instead of row index
%     if ismember('time', data.Properties.VariableNames)
%         time = data.time;
%     else
%         time = (1:height(data))'; % row index as time
%     end
% 
%     plot(time, data.cx, '-', 'LineWidth', 1.5, 'Color', colors(k,:));
%     plot(time, data.cy, '--', 'LineWidth', 1.5, 'Color', colors(k,:));
% end
% 
% xlabel('Time (s)');
% ylabel('Coefficient');
% % legendEntries = cell(1,2*numFiles);
% % for k = 1:numFiles
% %     [~,name,ext] = fileparts(filenames{k});
% %     legendEntries{2*k-1} = sprintf('Cd (%s)', [name ext]);
% %     legendEntries{2*k}   = sprintf('Cl (%s)', [name ext]);
% % end
% % legend(legendEntries, 'Interpreter', 'none');
% title('Drag (C_d) and Lift (C_l) vs Time for 3 Files');
% grid on;
% hold off;

% ---- Step 3: Plot only Drag Coefficient (Cd = cx) for all 3 files ----
figure; hold on;
colors = lines(numFiles); % different colors for clarity

for k = 1:numFiles
    data = allData{k};
    
    % If you have a 'time' column, use it instead of row index
    if ismember('time', data.Properties.VariableNames)
        time = data.time;
    else
        time = (1:height(data))'; % row index as time
    end
    
    % Plot only cx (Cd)
    plot(time, data.cx, '-', 'LineWidth', 1.5, 'Color', colors(k,:));
end

xlabel('Time');
ylabel('Drag Coefficient C_D');
legendEntries = cell(1,numFiles);
for k = 1:numFiles
    [~,name,ext] = fileparts(filenames{k});
    legendEntries{k} = sprintf('Cd (%s)', [name ext]);
end
legend(legendEntries, 'Interpreter', 'none');
title('Drag Coefficient (C_D) vs Time for Different Cases');
grid on;
hold off;
