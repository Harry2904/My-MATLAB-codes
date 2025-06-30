% Given data
w = 2e-5;
F = [7.958e3,15.916e3,23.874e3,31.832e3,39.79e3];

% Diameter of craters for Soda Glass (selecting GD100 for processing)
Gl_D100 = [23.6e-6,28.9e-6,32.7e-6,31e-6,47.1e-6];

% Square the GD values
Gl_D_squared = Gl_D100.^2;

% Natural logarithm of force
log_F = log(F);

% Scatter plot
figure;
scatter(log_F, Gl_D_squared, 'b', 'DisplayName', 'Data Points');
hold on;

% Linear fit
p = polyfit(log_F, Gl_D_squared, 1);
slope = p(1);
intercept = p(2);

% Generate fit line data
log_F_th = linspace(min(log_F), max(log_F), 100);
Gl_D_squared_fit = polyval(p, log_F_th);

% Plot the linear fit
plot(log_F_th, Gl_D_squared_fit, 'r', 'DisplayName', ...
    sprintf('Linear Fit: y = %.2e x + %.2e', slope, intercept));

% Label the axes
xlabel('Log(F/F_{th})');
ylabel('Squared D (m²)');
title('Squared D vs Log(F/F_{th}) with Linear Fit');
legend('Location', 'best');
grid on;

% Display slope and intercept
fprintf('Slope: %.2e\n', slope);
fprintf('Intercept: %.2e\n', intercept);
