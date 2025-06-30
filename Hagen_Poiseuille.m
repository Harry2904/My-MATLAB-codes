clc; clear; close all;

%% 1. Input Parameters
L = 15;
R = 0.5;
Re = 100;
U_in = 1;
dt = 0.001;
nt = 10000;
convergence_criterion = 1e-3;
imax = 64;
jmax = 152;

%% 2. Non-Uniform Grid (Algebraic Stretching)
eta = linspace(0, 1, jmax);
r = R * eta.^2;
dr = diff(r); dr = [dr, dr(end)];
dz = L / (imax - 1);

%% 3. Corner Cells Arrays (Storage for Future Use)
Uz = zeros(jmax, imax);
Ur = zeros(jmax, imax);
P = zeros(jmax, imax);
P_prime = zeros(jmax, imax);

%% 4. Initialize U-Velocity and Its Components
Uz(:,1) = U_in;
Uz(jmax,:) = 0;
Uz(1,:) = Uz(2,:);

%% 5. Initialize V-Velocity and Its Components
Ur(:,1) = 0;
Ur(:,imax) = 0;
Ur(jmax,:) = 0;
Ur(1,:) = 0;

%% 6. Initialize Pressure and Its Components
P(:,:) = 0;
P_prime(:,:) = 0;

%% 7. Initialize Mass Source Term and Mass Flux
mass_residual = 1;

%% 8. Convergence Criterion
iter = 0;
while iter < nt && mass_residual > convergence_criterion
    iter = iter + 1;
    
    %% 9. Predictor Step (Solving Momentum Equations)
    % Update Uz using predictor formula
    j = 2;
    while j < jmax
        i = 2;
        while i < imax
            dUz_dt = -(Uz(j,i)*(Uz(j,i+1) - Uz(j,i-1))/(2*dz)) ...
                     -(Ur(j,i)*(Uz(j+1,i) - Uz(j-1,i))/(2*dr(j))) ...
                     +(1/Re)*((Uz(j,i+1) - 2*Uz(j,i) + Uz(j,i-1))/(dz^2) + ...
                              (1/r(j))*(Uz(j+1,i) - Uz(j-1,i))/(2*dr(j)) + ...
                              (Uz(j+1,i) - 2*Uz(j,i) + Uz(j-1,i))/(dr(j)^2));
            Uz(j,i) = Uz(j,i) + dt*dUz_dt;
            i = i + 1;
        end
        j = j + 1;
    end
    
    % Update Ur using predictor formula
    j = 2;
    while j < jmax
        i = 2;
        while i < imax
            dUr_dt = -(Uz(j,i)*(Ur(j,i+1) - Ur(j,i-1))/(2*dz)) ...
                     -(Ur(j,i)*(Ur(j+1,i) - Ur(j-1,i))/(2*dr(j))) ...
                     +(1/Re)*((Ur(j,i+1) - 2*Ur(j,i) + Ur(j,i-1))/(dz^2) + ...
                              (1/r(j))*(Ur(j+1,i) - Ur(j-1,i))/(2*dr(j)) + ...
                              (Ur(j+1,i) - 2*Ur(j,i) + Ur(j-1,i))/(dr(j)^2));
            Ur(j,i) = Ur(j,i) + dt*dUr_dt;
            i = i + 1;
        end
        j = j + 1;
    end
 end
