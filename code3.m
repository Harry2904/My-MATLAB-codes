% Given data
w = 2e-5;
F = [7.958e3,15.916e3,23.874e3,31.832e3,39.79e3];

% Diameter of craters for Stainless Steel (using SD100 for this analysis)
SD1 = [16.35e-6,22.85e-6,11e-6,12.25e-6,18.85e-6];
SD2 = [18.3e-6,21.45e-6,15.05e-6,20.45e-6,20.15e-6];
SD3 = [20.65e-6,34.2e-6,20.25e-6,22.95e-6,26.95e-6];
SD5 = [19.05e-6,34.2e-6,26.25e-6,27.95e-6,25.05e-6];
SD10 = [27.5e-6,35.85e-6,39.1e-6,35.35e-6,46.1e-6];
SD100 = [31.2e-6,40.85e-6,47.1e-6,51.45e-6,57.4e-6];

% Square the crater diameters
SD_squared = SD100.^2;

% Natural logarithm of force
log_F = log(F);

% Scatter plot
figure;
scatter(log_F, SD_squared, 'b', 'DisplayName', 'Data Points');
hold on;

% Linear regression
p = polyfit(log_F, SD_squared, 1);
slope = p(1);
intercept = p(2);

% Generate fit line
log_F_fit = linspace(min(log_F), max(log_F), 100);
SD_squared_fit = polyval(p, log_F_fit);

% Plot fit line
plot(log_F_fit, SD_squared_fit, 'r', 'DisplayName', ...
    sprintf('Linear Fit: y = %.2e x + %.2e', slope, intercept));

% Axis labels and title
xlabel('Log(F/F_{th})');
ylabel('Squared D (m²)');
title('Squared D vs Log(F/F_{th}) with Linear Fit');
legend('Location', 'best');
grid on;

% Output slope and intercept
fprintf('Slope: %.2e\n', slope);
fprintf('Intercept: %.2e\n', intercept);
