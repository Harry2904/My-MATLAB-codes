% Constants
L = 3.6;  % m
D = 0.025; %m
z = linspace(0, L, 500);
q0 = 600e3;  % W/m^2
h_f = 3622.22;  % W/m^2·K
T_inlet = 80;  % °C
G = 254.648; %kg/m^2.s
Cp_f = 5400; %J/kg.K

% Heat flux profiles
q_sine = q0 * sin(pi * z / L);  % W/m^2
q_const = q0 * ones(size(z));  % W/m^2

% Wall temperature rise for sine profile
T_wall_sine = T_inlet + q0*((1/h_f) + (4*L/(pi*G*D*Cp_f))*(1 - cos(pi*z/L)));

% Wall temperature rise for constant profile
T_wall_const = T_inlet + q0*((1/h_f) + (4*z/(G*D*Cp_f)));

% Plotting
figure('Position', [100, 100, 1000, 500]);
plot(z, T_wall_sine, 'b', 'DisplayName', 'Sine Heat Flux Profile');
hold on;
plot(z, T_wall_const, 'r--', 'DisplayName', 'Constant Heat Flux Profile');
xlabel('Axial Distance z (m)');
ylabel('Wall Temperature (°C)');
title('Wall Temperature Variation in Single-Phase Region');
legend('show');
grid on;
box on;
hold off;