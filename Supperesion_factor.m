% Given properties of R-22 at 218 kPa
rho_g = 13.8;      % kg/m^3
rho_f = 1057;      % kg/m^3
mu_g = 1.5e-5;     % Pa.s
mu_f = 2.5e-4;     % Pa.s
kf = 0.075;        % W/m.K
cp_f = 1.28e3;     % J/kg.K
hfg = 176e3;       % J/kg
sigma = 0.0136;    % N/m
Tsat = 233.9;      % K
q = 6e3;           % W/m^2
D = 0.007;         % m
Nu = 4.36;         % constant heat flux
% Preallocate arrays
G_list = 100:100:1000;              % kg/m^2.s
gamma_list = zeros(size(G_list));
xsupp_list = zeros(size(G_list));
hfo_list = zeros(size(G_list));

% Loop over G values
for i = 1:length(G_list)
    G = G_list(i);                 % Mass flux for current iteration

    % Reynolds and Prandtl numbers
    Re = G * D / mu_f;
    Pr = (mu_f * cp_f) / kf;

    % Dittus-Boelter Nusselt number
    Nu = 0.023 * Re^0.8 * Pr^0.4;

    % Heat transfer coefficient hfo
    hfo = Nu * kf / D;
    hfo_list(i) = hfo;

    % Gamma (with hfo varying)
    gamma = (rho_g / rho_f)^0.56 * (mu_f / mu_g)^0.11 * ...
        ((q * hfg * kf * rho_g) / (98 * hfo^2 * sigma * Tsat))^1.11;

    % Suppression quality
    xsupp = gamma / (1 + gamma);

    % Store results
    gamma_list(i) = gamma;
    xsupp_list(i) = xsupp;
end

% Display results in a table
results = table(G_list', hfo_list', gamma_list', xsupp_list', ...
    'VariableNames', {'G_kg_m2s', 'hfo_W_m2K', 'Gamma', 'x_supp'});
disp(results);

% Plot suppression quality vs gamma
figure;
plot(gamma_list, xsupp_list, 'bo-', 'LineWidth', 1.5);
xlabel('\gamma (dimensionless)');
ylabel('Suppression Quality x_{supp}');
title('Suppression Quality vs \gamma for R-22');
grid on;
writetable(results, 'suppression_quality_results.xlsx');