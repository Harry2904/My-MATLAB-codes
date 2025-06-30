clear all;
clc;

%% Properties of n-butyl alcohol (C4H9OH) at 1 atm
T_sat = 390.85;
h_fg = 5.91e5;
rho_f = 712;
rho_v = 2.30;
sigma = 0.0171;
Cp_f = 3200;
mu_f = 4.4e-4;
k_f = 0.127;
Pr_f = 10.3;
g = 9.81;
C2 = 2.5;

%% Rohsenow correlation parameters
C_sf = 0.003;
m = 0.7;
n = 0.33;

%% Superheat range (5°C to 30°C)
delta_T = 5:1:30;
delta_T_K = delta_T;
T_w = T_sat + delta_T_K;

%% Preallocate heat flux arrays
q_Rohsenow = zeros(size(delta_T));
q_SA = zeros(size(delta_T));

%% Calculate heat flux using both correlations
for i = 1:length(delta_T_K)
    lhs = (Cp_f * delta_T_K(i)) / h_fg;
    E = (sigma / (g * (rho_f - rho_v)))^0.5;

    % Rohsenow correlation
    q_Rohsenow(i) = ((lhs / (C_sf * Pr_f^(m + 1)))^(1/n)) * mu_f * h_fg / E;

    % Stephan-Abdelsalam correlation
    q_SA(i) = (C2 * delta_T(i))^(1/0.330);
end

%% Convert to kW/m^2 for plotting
q_Rohsenow_kW = q_Rohsenow / 1e3;
q_SA_kW = q_SA / 1e3;

disp(q_SA_kW);
disp(q_Rohsenow_kW)

%% Plot the boiling curves
figure;
plot(delta_T, q_Rohsenow_kW, 'b-', 'LineWidth', 2, 'DisplayName', 'Rohsenow');
hold on;
plot(delta_T, q_SA_kW, 'r--', 'LineWidth', 2, 'DisplayName', 'Stephan-Abdelsalam');
grid on;
xlabel('Wall Superheat (\DeltaT) [°C]');
ylabel('Heat Flux q'''' [kW/m^2]');
title('Nucleate Pool Boiling Curve of n-Butyl Alcohol on Copper at 1 atm');
legend('Location', 'northwest');
xlim([5, 30]);
ylim([0, max([q_Rohsenow_kW q_SA_kW]) * 1.1]);
