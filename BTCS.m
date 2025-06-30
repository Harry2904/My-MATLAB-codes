% Parameters
L = 1.0;  % Length of the bar (m)
dx = 0.1; % Spatial step (m)
dt = 1; % Time step (s)
alpha = 400 / (8000 * 385); % Thermal diffusivity (m^2/s)
x = 0:dx:L;  % Spatial grid
Nx = length(x); % Number of spatial points
T = 900;  % Total time (s)
Nt = T/dt; % Number of time steps

% Initial and boundary conditions
T_initial = 25;  % Initial temperature (°C)
T_boundary = 400; % Boundary temperature (°C)

% BTCS Scheme
T_btcs = T_initial * ones(Nx, 1); 
T_btcs(1) = T_boundary;  
T_btcs(end) = T_boundary;

% Coefficient for BTCS scheme
r_btcs = alpha * dt / dx^2;

% Coefficient matrix for the implicit scheme (tridiagonal matrix A)
A = (1 + 2*r_btcs) * eye(Nx-2) + (-r_btcs) * diag(ones(Nx-3,1),1) + (-r_btcs) * diag(ones(Nx-3,1),-1);

% Preallocate for temperature at middle over time
T_mid_btcs = zeros(1, Nt);

% BTCS loop
for n = 1:Nt
    % Prepare the right-hand side of the equation
    b = T_btcs(2:Nx-1);  % Internal points only
    
    % Adjust the boundary conditions
    b(1) = b(1) + r_btcs * T_btcs(1);  % Left boundary condition
    b(end) = b(end) + r_btcs * T_btcs(end);  % Right boundary condition
    
    % Solve the linear system A * T_new_internal = b
    T_new_internal = A \ b;
    
    % Update the temperature profile
    T_btcs(2:Nx-1) = T_new_internal;
    
    % Store temperature at the middle length
    T_mid_btcs(n) = T_btcs(floor(Nx/2));
end

% Plot BTCS temperature at the middle over time
figure;
plot(0:dt:T-dt, T_mid_btcs, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('BTCS: Temperature at Middle Length Over Time');
grid on;

