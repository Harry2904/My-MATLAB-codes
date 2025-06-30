% Incubation Coefficients Data
N = [1,2,3,5,10,100];
Fthg = [6.033e3,5.957e3,5.856e3,6.311e3,4.866e3,4.823e3];
Fths = [5890e3,2.324e3,15320e3,0.036e3,2.556e3,4.16e3];

F_th = Fths;  % Select the data for fitting

% Compute logarithms
log_F = log(F_th);
log_N = log(N);

% Scatter plot
figure;
scatter(log_N, log_F, 'b', 'DisplayName', 'Data Points');
hold on;

% Linear regression
p = polyfit(log_N, log_F, 1);
slope = p(1);
intercept = p(2);

% Fit line
log_N_fit = linspace(min(log_N), max(log_N), 100);
log_F_fit = polyval(p, log_N_fit);

% Plot linear fit
plot(log_N_fit, log_F_fit, 'r', 'DisplayName', ...
    sprintf('Linear Fit: y = %.2e x + %.2e', slope, intercept));

% Labeling
xlabel('log(N)');
ylabel('log(F_{th})');
title('log(F_{th}) vs log(N) with Linear Fit');
legend('Location', 'best');
grid on;

% Display results
fprintf('Slope: %.2e\n', slope);
fprintf('Intercept: %.2e\n', intercept);
