a = 250;  
x_min = 0; x_max = 400;  
dx = 0.1;  
x = x_min:dx:x_max;  
Nx = length(x);
c = 0.1;
dt_lw = c * dx / a; 
dt_btcs = 0.1; 
t_final = 0.5;
Nt_lw = round(t_final / dt_lw);
Nt_btcs = round(t_final / dt_btcs);

u_initial = zeros(1, Nx);
for i = 1:Nx
    if x(i) >= 50 && x(i) <= 110
        u_initial(i) = 100 * sin(pi * (x(i) - 50) / 60);
    end
end
u_initial(1) = 0;
u_initial(end) = 0;


%% Lax-Wendroff (Explicit)
u_lw = u_initial;  
for n = 1:Nt_lw
    u_new = u_lw; 
    for i = 2:Nx-1
        u_new(i) = u_lw(i) - 0.5 * c * (u_lw(i+1) - u_lw(i-1)) ...
            + 0.5 * c^2 * (u_lw(i+1) - 2*u_lw(i) + u_lw(i-1));
    end
    u_lw = u_new; 
end


%% BTCS (Implicit)
u_btcs = u_initial;
c = (a * dt_btcs) / (2 * dx);
for i = 1:Nt_btcs
    u_old = u_btcs;
    A = zeros(Nx-2, Nx-2);
    B = zeros(1, Nx-2);
    A(1, 1:2) = [1, c];
    A(Nx-2, Nx-3:Nx-2) = [-c, 1];
    for j = 2:Nx-3
        A(j, j-1:j+1) = [-c, 1, c];
        B(j) = u_old(j);
    end
    u_btcs = [0, (A \ B')', 0];
end


%% Plotting the results
plot(x, u_lw, 'r', 'DisplayName', 'Lax-Wendroff');
hold on;
plot(x, u_btcs, 'b--', 'DisplayName', 'BTCS');
legend('show');
xlabel('x');
ylabel('u');
title('Comparison of Lax-Wendroff and BTCS at t = 0.5s');