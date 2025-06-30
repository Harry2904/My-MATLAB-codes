% Input parameters
R_i = 0.05; %in(m)
R_o = 0.2;  %in(m)
L = 2;      %in(m)
imax = 42;  
jmax = 22;  
rho = 7850; %(kg/m^3)
cp = 500;   %(J/kg-K)
k = 45;     %(W/m-K)
h = 25;     %(W/m^2-K)
beta = 1.2; % Non-uniform grid control parameter
q_flux = 500e4; %(W/m^2)
epsilon_st = 1e-4;
T_wb = 20;  % Left boundary temperature (°C)
T_ini = 50; % Initial temperature (°C)
T_inf_avg = 80; % Average ambient temperature (°C)
T_inf_amp = 30; % Ambient temperature amplitude (°C)

% Derived properties
alpha = k / (rho * cp);
Dz = L / (imax - 2); % Uniform axial step size
zi = linspace(0, L, imax-1);

dt = 360;
total_time = 24 * 3600;
t_max = 3 * total_time;

% Grid generation for non-uniform radial direction
dr = (R_o - R_i);
yi = linspace(0, 1, jmax-1);

% Generate radial coordinates (non-uniform spacing)
y = zeros(1, jmax-1);
beta_1 = beta + 1;
beta_2 = beta - 1;
for j = 1:jmax-1
   beta_1_div_beta_2 = (beta_1 / beta_2) ^ (2 * yi(j) - 1);
   num = (beta_1 * beta_1_div_beta_2) - beta_2;
   den = 2 * (1 + beta_1_div_beta_2);
   y(j) = dr * (num / den);
end

% Calculate cell centers
yc = zeros(1, jmax);
yc(2:jmax-1) = (y(2:jmax-1) + y(1:jmax-2)) / 2;
yc(1) = y(1);
yc(jmax) = y(jmax-1);

% Compute differences Dy and dy
Dy_vals = zeros(1, jmax);
dy = zeros(1, jmax-1);
Dy_vals(1) = y(1);
Dy_vals(jmax) = y(jmax-1);
Dy_vals(2:jmax-1) = y(2:jmax-1) - y(1:jmax-2);
dy(1:jmax-1) = yc(2:jmax) - yc(1:jmax-1);

% Initialize coefficient matrices
aE = zeros(imax, jmax);
aN = zeros(imax, jmax);
aP = zeros(imax, jmax);
aP_0 = zeros(imax, jmax);

% Compute coefficients in axial and radial directions
for j = 2:jmax-1
    for i = 2:imax-1
        aE(i,j) = k * (Dy_vals(j-1)) / Dz;
        aN(i,j) = k * (Dz) / Dy_vals(j);
        aP_0(i,j) = rho * cp * (Dz) * Dy_vals(j-1) / dt;
        aP(i,j) = aP_0(i,j) + aE(i,j) + aE(i-1,j) + aN(i,j) + aN(i,j-1);
    end
end

% Initialize temperature field
T = ones(imax, jmax) * T_ini;
T(1, :) = T_wb;

% Time-stepping loop
unsteadiness_nd = 1; 
n = 0;

while unsteadiness_nd >= epsilon_st
    n = n + 1;
    
    % Update boundary conditions
    T(imax, :) = T(imax-1, :);
    T(1, :) = T_wb;
    
    % Apply ambient sinusoidal variation on outer surface
    T_inf_t = T_inf_avg + T_inf_amp * sin(2 * pi * dt / (24 * 3600));
    aS_conv = k / dy(jmax-1);
    aP_conv = aS_conv + h;

    T(:, jmax) = (aS_conv * T(:, jmax-1) + h * T_inf_t) / aP_conv;


    T_old = T;
    
    %source term
    b = zeros(imax, jmax);
    for j = 2:jmax-1
        for i = 2:imax-1
            b(i,j) = aP_0(i,j) * T_old(i,j) + q_flux *dr* Dy_vals(j);
        end
    end

    %implicit scheme
    Error = 1;
    while Error >= epsilon_st
        T_old_iter = T;
        
        for j = 2:jmax-1
            for i = 2:imax-1
                T(i,j) = (aE(i,j) * T(i+1,j) + aE(i-1,j) * T(i-1,j) + ...
          aN(i,j) * T(i,j+1) + aN(i,j-1) * T(i,j-1) + b(i,j)) ...
         / aP(i,j);

            end
        end
        
        Error = max(max(abs(T - T_old_iter)));
    end
    
    unsteadiness_nd = max(max(abs(T - T_old))) / dt;

end

disp('Simulation complete!');

%% Plotting
time_steps = [3600, 2*3600, 72*3600];
colors = {'r', 'b', 'g'};

figure;
hold on;
for t_idx = 1:length(time_steps)
    
    T_radial = T(round(imax/2), :);
    
    plot(yc, T_radial, 'Color', colors{t_idx}, 'LineWidth', 2, ...
         'DisplayName', sprintf('t = %.1f hrs', time_steps(t_idx) / 3600));
end
xlabel('Radial Position (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution in Radial Direction');
legend;
grid on;
hold off