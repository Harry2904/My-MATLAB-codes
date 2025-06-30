% Input parameters
rho = 1000; % in kg/m^3
cp = 4187;  % in J/kgK
k = 1.5;    % in W/m-K
L = 6;      % in m
H = 1;      % in m
u = 1;      % in m/s
v = 0;      % in m/s
alpha = 3.58e-7; % in m^2/s
imax = 42;
jmax = 22;
epsilon_st = 1e-4; % Convergence tolerance
beta = 1.2;
Re = 10;
Pr = 1;

% Non-uniform grid generation in X-direction
zi = linspace(0, 1, imax-1);
beta_1 = beta + 1;
beta_2 = beta - 1;

z = L * ((beta_1 - (beta_2 * ((beta_1 / beta_2) .^ (zi)))) ./ ...
        (1 + ((beta_1 / beta_2) .^ (zi))));  

% Corner points in X-direction
zc = zeros(1, imax);
zc(2:imax-1) = (z(2:imax-1) + z(1:imax-2)) / 2; 
zc(1) = 0;
zc(imax) = L;
Dz(2:imax -1) = z(2:imax-1) - z(1:imax-2);
Dz(1) =0;Dz(imax)=0; 
Dy=Dz;
Dz=Dz';

% Non-uniform grid generation in Y-direction
yi = linspace(0, 1, jmax-1);

y = H * ((beta_1 * ((beta_1 / beta_2) .^ (2 * yi - 1)) - beta_2) ./ ...
         (2 * (1 + ((beta_1 / beta_2) .^ (2 * yi - 1)))));

% Corner points in Y-direction
yc = zeros(1, jmax);
yc(2:jmax-1) = (y(2:jmax-1) + y(1:jmax-2)) / 2;
yc(1) = 0;
yc(jmax) = H;

% Meshgrid
[Zc, Yc] = meshgrid(zc, yc);
[Z, Y] = meshgrid(z, y);

% Plot the grid
figure;
plot(Z', Y', 'k'); 
hold on;
plot(Z, Y, 'k');
xlabel('X');
ylabel('Y');
title('Non-Uniform Grid');
grid on;

% Initializing scheme
disp("SELECT THE ADVECTION SCHEME (1/2)");
fprintf("1. FOU \n");
fprintf("2. QUICK \n");
scheme = input("ENTER THE ADVECTION SCHEME : ");

%Mass flux initiation
mz = rho * u;
my = rho * v;
mz_plus = max(mz,0);
mz_minus = min(mz,0);
my_plus = max(my,0);
my_minus = min(my,0);

%weights in non-uniform grid
w_ne_p = zeros(imax-1,3);
w_ne_m = zeros(imax-1,3);
w_nn_p = zeros(jmax-1,3);
w_nn_m = zeros(jmax-1,3);

for i = 2:imax-2
    w_ne = weights(scheme, Dz(i+1), Dz(i), Dz(i-1));
    w_ne_p(i,3) = w_ne';

    w_ne = weights(scheme, Dz(i), Dz(i+1), Dz(i+2));
    w_ne_m(i,3) = w_ne';
end
for j = 2:jmax-2
    w_nn = weights(scheme, Dy(j+1), Dy(j), Dy(j-1));
    w_nn_p(j,3) = w_nn';

    w_nn = weights(scheme, Dy(j), Dy(j+1), Dy(j+2));
    w_nn_m(j,3) = w_nn';
end

%Non-dimensionalizing temperature(IC and BCs)

theta = zeros(imax, jmax); %initial condition
theta(1, :) = 1; %left wall
theta(:, 1) = 0; %bottom wall
theta(1:imax, jmax) = 0; %top wall

% Convergence criteria
unsteadiness_nd = 1;
n = 0;
max_iter = 5000;

while unsteadiness_nd > epsilon_st && n < max_iter
    n = n + 1;
    theta_old = theta;
 end

% Plot results
figure;
contourf(Zc, Yc, theta', 20, 'LineColor', 'none');
colorbar;
xlabel('X'); ylabel('Y');
title('Steady-State Temperature Contour');

plot(z ,y ,'r-', 'LineWidth', 2);
hold on;
xlabel('T');
ylabel('y');
grid on;