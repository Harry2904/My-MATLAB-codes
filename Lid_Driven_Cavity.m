% Grid parameters
L = 1;
Nx = 42;
Ny = 42;

dx = L / Nx;
dy = L / Ny;

% Staggered grid definitions
xP = linspace(dx/2, L-dx/2, Nx);
yP = linspace(dy/2, L-dy/2, Ny);

xU = linspace(0, L, Nx+1);
yU = linspace(dy/2, L-dy/2, Ny);

xV = linspace(dx/2, L-dx/2, Nx);
yV = linspace(0, L, Ny+1);
[Xp, Yp] = meshgrid(xP, yP);
[Xu, Yu] = meshgrid(xU, yU);
[Xv, Yv] = meshgrid(xV, yV);

figure;
hold on; grid on;
plot(Xp, Yp, 'ro', 'MarkerSize', 6, 'DisplayName', 'Pressure Nodes (P)');
plot(Xu, Yu, 'bs', 'MarkerSize', 6, 'DisplayName', 'U-Velocity Nodes (U)');
plot(Xv, Yv, 'g^', 'MarkerSize', 6, 'DisplayName', 'V-Velocity Nodes (V)');
xlabel('X-Direction'); ylabel('Y-Direction');
title('Staggered Grid for Lid-Driven Cavity Flow');
axis equal;
hold off;
