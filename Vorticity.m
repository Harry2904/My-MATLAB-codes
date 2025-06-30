L = 1;
Re = 100;
U = 1;
dx = 0.02;
dy = dx;
nx = (L / dx) + 1;
ny = (L / dx) + 1;
error = 10000;

psi = zeros(nx, ny);
omega = zeros(nx, ny);
u = zeros(nx, ny);
u(1, :) = U;
v = zeros(nx, ny);

k = 0;
while k < 1000
    psi_old = psi;
    omega_old = omega;
    u_old = u;
    v_old = v;

    if error == 10000
        A(:,:,1) = psi_old;
        A(:,:,2) = omega_old;
        A(:,:,3) = u_old;
        A(:,:,4) = v_old;
        M_old = max(abs(A),[],'all');
    else
        M_old = M_current;
    end

    % Update u and v
    for j = 2:nx-1
        for i = 2:ny-1
            u(i,j) = (psi(i-1, j) - psi(i+1, j)) / (2 * dy);
            v(i,j) = (psi(i, j-1) - psi(i, j+1)) / (2 * dx);
        end
    end


    % Update streamfunction si
    for j = 2:nx-1
        for i = 2:ny-1
            psi(i,j) = (psi(i+1, j) + psi(i-1, j) + psi(i, j+1) + psi(i, j-1) + (omega(i, j)* dx^2)) / 4;
        end
    end
        
    % Boundary conditions for omega
    omega(1, :) = (-2 * psi(2, :) - 2*U*dy) / dy^2;  %bottom
    omega(nx, :) = -((2 * psi(nx-1, :))) / dy^2;  %top
    omega(:, 1) = (-2 * psi(:, 2)) / dx^2;  %left
    omega(:, ny) = (-2 * psi(:, ny-1)) / dx^2;  %right

    % Update vorticity omega
    for j = 2:nx-1
        for i = 2:ny-1
            w = (1-(u(i,j)*dx*Re/2))*omega(i,j+1) + (1+(u(i,j)*dx*Re/2))*omega(i,j-1) + (1-(v(i,j)*dx*Re/2))*omega(i-1,j) + (1+(v(i,j)*dx*Re/2))*omega(i+1,j);
            omega(i,j) = w / 4;
        end
    end
% Compute error and update the maximums
    B(:,:,1) = psi;
    B(:,:,2) = omega;
    B(:,:,3) = u;
    B(:,:,4) = v;
    M_current = max(abs(B),[],'all');
    error = abs(M_current - M_old);
    k = k + 1;
end

% Plotting
figure;
contourf(linspace(0, L, nx), ...
    -linspace(0, L, ny), psi, 40);
colorbar;
title('Streamfunction \psi');
xlabel('x'); ylabel('y');

figure;
contourf(linspace(0, L, nx), -linspace(0, L, ny), omega, 40);
colorbar;
title('Vorticity \omega');
xlabel('x'); ylabel('y');

figure;
contourf(linspace(0, L, nx), -linspace(0, L, ny), u, 40);
colorbar;
title('Velocity u');
xlabel('x'); ylabel('y');

figure;
contourf(linspace(0, L, nx), -linspace(0, L, ny), v, 40);
colorbar;
title('Velocity v');
xlabel('x'); ylabel('y');

figure;
plot(u(:,26), 1-linspace(0, L, ny), 'r', 'LineWidth', 2);
title('Velocity Profile at Centreline');
xlabel('u');
ylabel('y');