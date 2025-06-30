% Given data
w = 2e-5;
F = [7.958e3,15.916e3,23.874e3,31.832e3,39.79e3];

% Depths of craters for Soda Glass
h5 = [1.115e-6,1.439e-6,6.622e-6,1.51e-6,2.398e-6];
h100 = [18.248e-6,22.062e-6,31.417e-6,35.314e-6,35.786e-6];

% Choose the h data to use
h = h100;

% Compute the natural logarithm of F
log_F = log(F);

% Scatter plot
figure;
scatter(log_F, h, 'b', 'DisplayName', 'Data Points');
hold on;

% Linear fit
p = polyfit(log_F, h, 1);
slope = p(1);
intercept = p(2);

% Generate line data for fit
log_F_th = linspace(min(log_F), max(log_F), 100);
h_fit = polyval(p, log_F_th);

% Plot the linear fit
plot(log_F_th, h_fit, 'r', 'DisplayName', ...
    sprintf('Linear Fit: y = %.2e x + %.2e', slope, intercept));

% Labels and title
xlabel('Log(F/F_{th})');
ylabel('Depth h (m)');
title('Depth h vs Log(F/F_{th}) with Linear Fit');
legend('Location', 'best');
grid on;

% Display slope and intercept
fprintf('Slope: %.2e\n', slope);
fprintf('Intercept: %.2e\n', intercept);
