T = 20.000; % °C
P = 572250; % Pa
L = 3;
x = 0.1:0.1:0.8;

D = input('Enter pipe diameter (m): ');
mass_flow = input('Enter mass flow rate (kg/s): ');
dT_sat = input('Enter temperature difference (K): ');
dp_sat = input('Enter pressure difference Psat(t_w) - Psat(T_sat) (Pa): ');

% Liquid phase
liq.rho = 1224.4;        % Density (kg/m^3)
liq.v = 0.00081610;      % Specific volume (m^3/kg)
liq.h = 227500;          % Enthalpy (kJ/kg)
liq.cp = 1405;         % Cp (J/g·K)
liq.mu = 0.000207;     % Viscosity (Pa·s)
liq.k = 0.083284;        % Thermal conductivity (W/m·K)
sigma = 0.00869;       % Surface tension (N/m)

% Vapor phase
vap.rho = 27.791;
vap.v = 0.035997;
vap.h = 410000;
vap.cp = 1000.7;
vap.mu = 1.1488e-05;
vap.k = 0.013335;

hfg = vap.h - liq.h; 


hc = zeros(size(x));
Re_f = zeros(size(x));
Re_2phi = zeros(size(x));
Xtt = zeros(size(x));
inv_Xtt = zeros(size(x));
F = zeros(size(x));
h_NB = zeros(size(x));
h_total = zeros(size(x));


for i = 1:length(x)
    % Martinelli parameter
    Xtt(i) = ((1 - x(i))/x(i))^0.9 * (vap.rho / liq.rho)^0.5 * (liq.mu / vap.mu)^0.1;
    inv_Xtt(i) = 1 / Xtt(i);
    
    % Correction factor F
    if inv_Xtt(i) < 0.1
        F(i) = 1;
    else
        F(i) = 2.35 * (0.213 + inv_Xtt(i))^0.736;
    end
    
    A = (pi / 4) * D^2;
    G = mass_flow / A;
    
    % Convective heat transfer coefficient
    hc(i) = 0.023 * ((G * (1 - x(i)) * D) / liq.mu)^0.8 * ...
            ((liq.mu * liq.cp / liq.k)^0.4) * (liq.k / D) * F(i);
    
    % Two-phase Reynolds number
    Re_f(i) = (G * (1 - x(i)) * D) / liq.mu;
    Re_2phi(i) = F(i)^(1/0.8) * Re_f(i);
    
    % Suppression factor
    S = 1 / (1 + 2.53e-6 * Re_2phi(i)^1.17);
    
    % Nucleate boiling heat transfer coefficient
    h_NB(i) = 0.00122 * S * (dT_sat)^0.24 * (dp_sat)^0.75 * ...
              (liq.k^0.79 * liq.rho^0.49 * liq.cp^0.45) / ...
              (sigma^0.5 * liq.mu^0.29 * hfg^0.24 * vap.rho^0.24);
    
    h_total(i) = hc(i) + h_NB(i);
end


figure;
plot(x, hc, '-o', 'DisplayName', 'Convective h_c');
hold on;
plot(x, h_NB, '-s', 'DisplayName', 'Nucleate Boiling h_{NB}');
plot(x, h_total, '-d', 'DisplayName', 'Total h_{total}');
hold off;
% Customize plot with formatted title
title(sprintf('Heat Transfer Coefficients vs Quality\nD = %.4f m, m = %.4f kg/s, dT_sat = %.4f K', ...
              D, mass_flow, dT_sat), 'Interpreter', 'none');

xlabel('Quality (x)');
ylabel('Heat Transfer Coefficient (W/m^2·K)');
legend('Location', 'best');
grid on;


disp('Results for different x values:');
disp('x        hc          h_NB        h_total');
for i = 1:length(x)
    fprintf('%.1f    %.6f    %.6f    %.6f\n', x(i), hc(i), h_NB(i), h_total(i));
end