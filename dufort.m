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

% Du-Fort Frenkel Scheme
T_dufort = T_initial * ones(Nx, 1); 
T_old = T_dufort; % For holding the previous time step
T_dufort(1) = T_boundary;
T_dufort(end) = T_boundary;

% Coefficient for Du-Fort Frenkel scheme
r_dufort = 2 * alpha * dt / dx^2;

% Preallocate for temperature at middle over time
T_mid_dufort = zeros(1, Nt);

% Du-Fort Frenkel loop
for n = 1:Nt
    T_new = T_dufort; % Temporary storage
    for i = 2:Nx-1
        T_new(i) = (1 / (1 + r_dufort)) * (r_dufort * (T_dufort(i+1) + T_dufort(i-1)) + (1 - r_dufort) * T_old(i));
    end
    T_old = T_dufort;
    T_dufort = T_new;
    
    % Store temperature at the middle length
    T_mid_dufort(n) = T_dufort(floor(Nx/2));
end

% Plot Du-Fort Frenkel
figure;
plot(0:dt:T-dt, T_mid_dufort, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Du-Fort Frenkel: Temperature at Middle Length Over Time');
grid on;
