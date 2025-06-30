% Parameters
L = 1.0;  % Length of the bar (m)
dx = 0.1; % Spatial step (m)
dt = 0.1; % Time step (s)
alpha = 400 / (8000 * 385); % Thermal diffusivity (m^2/s)
x = 0:dx:L;  % Spatial grid
Nx = length(x); % Number of spatial points
T = 900;  % Total time (s)
Nt = T/dt; % Number of time steps

% Initial and boundary conditions
T_initial = 25;  % Initial temperature (°C)
T_boundary = 400; % Boundary temperature (°C)

% FTCS Scheme
T_ftcs = T_initial * ones(Nx, 1); 
T_ftcs(1) = T_boundary;  % Left boundary
T_ftcs(end) = T_boundary; % Right boundary

% Coefficient for FTCS scheme
r = alpha * dt / dx^2;

% Preallocate for temperature at middle over time
T_mid_ftcs = zeros(1, Nt);

% FTCS loop
for n = 1:Nt
    T_new = T_ftcs;  % Temporary storage for updated values
    for i = 2:Nx-1
        T_new(i) = T_ftcs(i) + r * (T_ftcs(i+1) - 2 * T_ftcs(i) + T_ftcs(i-1));
    end
    T_ftcs = T_new;
    
    % Store temperature at the middle length
    T_mid_ftcs(n) = T_ftcs(floor(Nx/2));
end

% Plot FTCS
figure;
plot(0:dt:T-dt, T_mid_ftcs, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('FTCS: Temperature at Middle Length Over Time');
grid on;

