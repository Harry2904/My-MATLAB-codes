imax = 42;
jmax = 22;

R_inner = 0.05; 
R_outer = 0.2;  
L = 2.0;
h = 0.2 - 0.05;
beta = 1.2;
rho = 7850; 
Cp = 500;  
k = 45;    
dt = 360; 
t_end = 72 * 3600; 
time_steps = t_end / dt;
error_steady = 100;


% Hypothetical Commputational space
delta_zeta = 1 / (jmax - 2);
delta_xi = 1 / (imax - 2);
zeta = delta_zeta * linspace(0, jmax-2, jmax-1);
xi = delta_xi * linspace(0, imax-2, imax-1);


% Calculate vertices
term = ((beta + 1) / (beta - 1)).^(1 - xi);
x = L * ((beta + 1) - (beta - 1) * term) / (1 + term);
r = h * ((1 + beta) * ((beta + 1) / (beta - 1)).^(2*zeta - 1) + (1 - beta)) ./ ...
    (2 * (1 + ((beta + 1) / (beta - 1)).^(2*zeta - 1)));

figure;
plot(x, r, '-o', 'MarkerSize', 3, 'LineWidth', 1.5)
xlabel('x (Non-uniform grid)')
ylabel('r (Non-uniform grid)')
title('Non-uniform Grid Distribution')
grid on