%% Staggered Grid Navier-Stokes Solver for Lid-Driven Cavity Flow
% Based on the formulation from "Computational Fluid Dynamics on a Staggered Grid"

clear; clc;

%% Grid and Flow Parameters
Lx = 1;  % Domain length (x-direction)
Ly = 1;  % Domain height (y-direction)
Nx = 42;  % Number of pressure grid points in x-direction
Ny = 42;  % Number of pressure grid points in y-direction
Re = 100; % Reynolds number
nu = 1 / Re;  % Kinematic viscosity (nondimensional)
U_lid = 1; % Lid velocity

%% Staggered Grid Generation
dx = Lx / Nx; % Grid spacing in x
dy = Ly / Ny; % Grid spacing in y

% Velocity grid points (u, v) are staggered w.r.t pressure (p)
x_u = linspace(0, Lx, Nx+1); % u-velocity at cell faces
y_v = linspace(0, Ly, Ny+1); % v-velocity at cell faces
x_p = linspace(dx/2, Lx-dx/2, Nx); % Pressure at cell centers
y_p = linspace(dy/2, Ly-dy/2, Ny);

%% Initialization
u_star = zeros(Nx+1, Ny);  % Intermediate u-velocity
v_star = zeros(Nx, Ny+1);  % Intermediate v-velocity
p = zeros(Nx, Ny);         % Pressure field

%% Time-stepping parameters
dt = min([dx, dy])^2 / (4*nu); % CFL condition for stability
max_iter = 1000;

%% Main Time-Stepping Loop
for iter = 1:max_iter
    % Step 1: Solve Intermediate Momentum Equations
    u_star = compute_nu_star(u_star, v_star, p, Nx, Ny, dx, dy, dt, nu, U_lid);
    
    % Step 2: Solve Pressure-Poisson Equation
    div_vel = (u_star(2:end, :) - u_star(1:end-1, :)) / dx + ...
              (v_star(:, 2:end) - v_star(:, 1:end-1)) / dy;
    p_rhs = div_vel / dt;
    p = poisson_solver(p_rhs, Nx, Ny, dx, dy);
    
    % Step 3: Correct Velocities
    for i = 2:Nx
        for j = 2:Ny-1
            u_star(i,j) = u_star(i,j) - dt * (p(i,j) - p(i-1,j)) / dx;
            v_star(i,j) = v_star(i,j) - dt * (p(i,j) - p(i,j-1)) / dy;
        end
    end
end

%% Visualization
[X, Y] = meshgrid(x_p, y_p);
U_plot = 0.5 * (u_star(1:end-1, :) + u_star(2:end, :)); % Interpolate u to cell centers
V_plot = 0.5 * (v_star(:, 1:end-1) + v_star(:, 2:end)); % Interpolate v to cell centers
quiver(X, Y, U_plot', V_plot', 'k');
title('Velocity Vector Plot (Staggered Grid)'); xlabel('x'); ylabel('y');

%% Function to Compute nu_star (Intermediate Velocity)
function nu_star = compute_nu_star(nu_star, v_star, p, Nx, Ny, dx, dy, dt, nu, U_lid)
    for i = 2:Nx
        for j = 2:Ny-1
            u_conv = -(nu_star(i+1,j) - nu_star(i-1,j)) / (2*dx);
            v_conv = -(v_star(i,j+1) - v_star(i,j-1)) / (2*dy);
            laplacian_u = (nu_star(i+1,j) - 2*nu_star(i,j) + nu_star(i-1,j)) / dx^2 + ...
                          (nu_star(i,j+1) - 2*nu_star(i,j) + nu_star(i,j-1)) / dy^2;
            nu_star(i,j) = nu_star(i,j) + dt * (-u_conv - v_conv + nu * laplacian_u);
        end
    end
    
    % Apply Lid-Driven Cavity Boundary Conditions
    nu_star(1, :) = 0;  nu_star(end, :) = 0;
    nu_star(:, 1) = 0;  nu_star(:, end) = U_lid;
end

%% Poisson Solver Function (Using Gauss-Seidel)
function p = poisson_solver(rhs, Nx, Ny, dx, dy)
    p = zeros(Nx, Ny);
    tol = 1e-6;
    error = 1;
    while error > tol
        p_old = p;
        for i = 2:Nx-1
            for j = 2:Ny-1
                p(i,j) = 0.25 * (p_old(i+1,j) + p_old(i-1,j) + p_old(i,j+1) + p_old(i,j-1) - rhs(i,j) * dx * dy);
            end
        end
        error = max(max(abs(p - p_old)));
    end
end
