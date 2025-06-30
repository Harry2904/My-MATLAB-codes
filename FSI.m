% 2D Immersed Boundary Method - Flow over Cylinder at Re=100
clear; clc;

% Grid setup
Lx = 1.0; Ly = 1.0;
Nx = 200; Ny = 200;
dx = Lx / (Nx - 1); dy = Ly / (Ny - 1);
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[xx, yy] = meshgrid(x, y);

% Cylinder parameters
R = 0.05; xc = 0.3; yc = 0.5;
Re = 100; U_inf = 1.0;
nu = U_inf * 2 * R / Re;

% Time setup
dt = 0.001; t_final = 2.0;
nt = floor(t_final / dt);

% Field initialization
u = U_inf * ones(Ny, Nx);
v = zeros(Ny, Nx);
psi = zeros(Ny, Nx);
omega = zeros(Ny, Nx);

% Boundary points on cylinder
theta = linspace(0, 2*pi, 200);
Xb = xc + R * cos(theta);
Yb = yc + R * sin(theta);

% Time loop
for t = 1:nt
    % Step 1: Vorticity from velocity
    omega(2:end-1, 2:end-1) = (v(2:end-1,3:end) - v(2:end-1,1:end-2)) / (2*dx) ...
                            - (u(3:end,2:end-1) - u(1:end-2,2:end-1)) / (2*dy);

    % Step 2: Solve Poisson’s equation for streamfunction
    for k = 1:200
        psi(2:end-1,2:end-1) = 0.25 * ( ...
            psi(2:end-1,3:end) + psi(2:end-1,1:end-2) + ...
            psi(3:end,2:end-1) + psi(1:end-2,2:end-1) + ...
            dx^2 * omega(2:end-1,2:end-1));
    end

    % Step 3: Velocity from streamfunction
    u(2:end-1,2:end-1) = (psi(3:end,2:end-1) - psi(1:end-2,2:end-1)) / (2*dy);
    v(2:end-1,2:end-1) = -(psi(2:end-1,3:end) - psi(2:end-1,1:end-2)) / (2*dx);

    % Step 4: Direct forcing IBM to impose no-slip on boundary points
    for k = 1:length(Xb)
        i = round((Yb(k) - y(1)) / dy) + 1;
        j = round((Xb(k) - x(1)) / dx) + 1;
        if i >= 1 && i <= Ny && j >= 1 && j <= Nx
            u(i,j) = 0;
            v(i,j) = 0;
        end
    end

    % Step 5: Visualization
    if mod(t,100) == 0
        clf;
        contourf(xx, yy, sqrt(u.^2 + v.^2), 30, 'LineColor', 'none');
        hold on;
        fill(Xb, Yb, 'k'); % Cylinder
        quiver(xx(1:5:end,1:5:end), yy(1:5:end,1:5:end), ...
               u(1:5:end,1:5:end), v(1:5:end,1:5:end), 'r');
        title(['Time = ', num2str(t*dt), ' s']);
        colorbar; axis equal; drawnow;
    end
end
