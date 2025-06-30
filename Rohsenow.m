% Given properties of n-butyl alcohol at 1 atm
h_fg = 552000;         % J/kg
cp_f = 2550;           % J/kg.K
Pr_f = 13.5;
sigma = 0.024;         % N/m
rho_f = 810;           % kg/m^3
rho_g = 2.2;           % kg/m^3
mu_f = 0.0024;         % Pa.s
g = 9.81;              % m/s^2

% Rohsenow constants for copper-n-butyl alcohol
C_sf = 0.013;
n = 1.7;

% Wall superheat range (°C)
dT = linspace(5, 30, 100);

% Rohsenow correlation
term = (cp_f .* dT) ./ (h_fg .* C_sf .* Pr_f.^n);
q_roh = mu_f .* h_fg .* sqrt(g * (rho_f - rho_g) / sigma) .* (term).^n;

% Stephan–Abdelsalam correlation
q_sa = 0.00122 .* h_fg .* rho_f.^0.79 .* sigma.^0.45 .* cp_f.^(-0.49) .* mu_f.^(-0.29) .* (dT).^3.06;

% Plotting
figure;
plot(dT, q_roh/1000, 'LineWidth', 2); hold on;
plot(dT, q_sa/1000, '--', 'LineWidth', 2);
xlabel('Wall Superheat \DeltaT_{sat} (°C)');
ylabel('Heat Flux q'''' (kW/m^2)');
title('Boiling Curves for n-Butyl Alcohol on Copper');
legend('Rohsenow', 'Stephan–Abdelsalam', 'Location', 'NorthWest');
grid on;
