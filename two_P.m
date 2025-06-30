% Known and given values
q0 = 600e3;       % W/m^2
L_total = 3.6;    % m
L_tp = 1.8;      % two-phase length from previous result (approx)
dz = 0.05;        % discretization step
D = 0.025;        % m (pipe diameter)
z_vals = 0:dz:(L_tp + dz);  % positions from start of 2φ region

% Thermophysical properties
h_fg = 1513.6e3;  % J/kg
m_dot = 0.125;      % kg/s
P_tp = pi * D * dz;  % heat transfer area per segment

% Heat flux profile in two-phase region
q_vals = q0 * sin(pi * (L_total - L_tp + z_vals) / L_total);  % W/m^2

% Heat added in each dz
delta_Q = q_vals * P_tp;  % W (J/s)
delta_x = delta_Q / (m_dot * h_fg);  % Quality increment

% Cumulative quality
x_vals = cumsum(delta_x);
x_vals = min(x_vals, 1);  % Clamp x to 1 (saturated vapor)

% Create table (equivalent to pandas DataFrame)
df_quality = table(z_vals', q_vals', delta_Q', delta_x', x_vals', ...
    'VariableNames', {'z_m', 'q_W_m2', 'delta_Q_W', 'delta_x', 'x'});

% Display first 10 rows
disp('First 10 rows:');
disp(df_quality(1:min(10, height(df_quality)), :));

% Plotting
figure;
subplot(2,1,1);
plot(z_vals, q_vals, 'b-', 'LineWidth', 2);
xlabel('Axial Distance z (m)');
ylabel('Heat Flux q'''' (W/m^2)');
title('Heat Flux Profile in Two-Phase Region');
grid on;

subplot(2,1,2);
plot(z_vals, x_vals, 'r-', 'LineWidth', 2);
xlabel('Axial Distance z (m)');
ylabel('Quality x');
title('Quality Variation Along Channel');
grid on;