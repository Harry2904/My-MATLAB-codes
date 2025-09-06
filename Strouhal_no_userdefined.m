% Step 1: Pick the file
[filename, pathname] = uigetfile({'*.dat;*.txt;*.csv','Data files (*.dat,*.txt,*.csv)'}, ...
                                 'Select the data file');
if isequal(filename,0), error('No file selected.'); end
filepath = fullfile(pathname, filename);

% Step 2: Read header line and clean it
fid = fopen(filepath,'r');
assert(fid>0,'Cannot open file: %s', filepath);

firstLine = fgetl(fid);
if ~ischar(firstLine)
    fclose(fid); error('File is empty.');
end

% Remove UTF-8 BOM if present and the "variables =" prefix
firstLine = regexprep(firstLine, '^\xEF\xBB\xBF', ''); %BOM
firstLine = regexprep(firstLine, '^\s*variables\s*=\s*', '', 'ignorecase');

% Split header into variable names (commas and/or whitespace)
rawNames = regexp(firstLine, '[,\s]+', 'split');
rawNames = rawNames(~cellfun('isempty',rawNames)); % drop empties

% Step 3: Read numeric data (handles spaces/tabs/commas)
fmt = repmat('%f', 1, numel(rawNames)); % one %f per column
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

% Make safe/unique MATLAB variable names but keep time/cy intact
varNames = matlab.lang.makeValidName(rawNames, 'ReplacementStyle', 'delete');
varNames = matlab.lang.makeUniqueStrings(varNames);
data = array2table(M, 'VariableNames', varNames);

% Step 4: Verify time & cy exist and proceed
if ~ismember('time', data.Properties.VariableNames) || ~ismember('cy', data.Properties.VariableNames)
    error('Columns "time" and/or "cy" not found. Found: %s', strjoin(data.Properties.VariableNames, ', '));
end

% User inputs stable time range
t_start = input('Enter start time of stable region: ');
t_end   = input('Enter end time of stable region: ');

% Extract time column
time = data.time;

% Logical mask for time range
mask = (time >= t_start & time <= t_end);

% Filter time and cy within that time range
time_values = data.time(mask);
cy_values = data.cy(mask);
cy_rms = rms(cy_values);

% Step 5: Define Flow Parameters
U_inf = 1; % Freestream velocity (dimensionless)
D = 1; % Characteristic length, (dimensionless)

% Step 6: Preprocess the signal
% Ensure the signal is suitable for FFT
dt = time(2) - time(1); % Calculate time step (s)
Fs = 1 / dt; % Sampling frequency (Hz)
N = length(cy_values); % Number of data points

% Remove the mean to eliminate the DC (0 Hz) component
cl_detrended = cy_values - mean(cy_values);

% Step 7: Perform the FFT
Y = fft(cl_detrended); % Compute the FFT
P2 = abs(Y / N); % Two-sided spectrum magnitude
P1 = P2(1:floor(N/2)+1); % Take first half of the spectrum
P1(2:end-1) = 2 * P1(2:end-1); % Scale the magnitudes correctly (except DC)

% Create the frequency axis (f) for the one-sided spectrum
f = Fs * (0:(N/2)) / N;

% Step 8: Find the Dominant Vortex Shedding Frequency
% Find the peak frequency in the spectrum, ignoring very low frequencies
min_freq_index = find(f > 0.16, 1); % Adjust f as needed
[~, peak_index] = max(P1(min_freq_index:end));
peak_index = peak_index + min_freq_index - 1; % Adjust index

f_shed = f(peak_index); % The dominant shedding frequency (Hz)

% Step 9: Calculate the Strouhal Number
Strouhal_number = f_shed * D / U_inf;

% Step 10: Display and Plot Results
fprintf('Vortex Shedding Frequency (f): %.4f Hz\n', f_shed);
fprintf('Strouhal Number (St): %.4f\n', Strouhal_number);

figure;
subplot(2, 1, 1)
plot(time_values, cy_values);
xlabel('Time (s)');
ylabel('Lift Coefficient (C_L)');
title('Time Domain Signal');
grid on;

subplot(2, 1, 2)
plot(f, P1);
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
title('Single-Sided Amplitude Spectrum of C_L');
grid on;
xlim([0, 5]); % Adjust x-axis limit to zoom in on relevant frequencies
hold on;
plot(f_shed, P1(peak_index), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Mark the peak
legend('Spectrum', 'Dominant Frequency');