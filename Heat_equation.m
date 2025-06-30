L_x = 1; 
L_y = 1; 
N_x = 40;
N_y = 40;
dx = L_x / (N_x-1);
dy = L_y / (N_y-1);
alpha = 0.0001;
dt = 0.1;
total_time = 600;
Nt = total_time / dt;
k = 100;
h = 10;
T_inf = 25;

T = zeros(N_x, N_y);

T(1,:) = 200;
T(:,end) = T_inf;
q_right = 2;

T_02 = zeros(1, Nt);
T_05 = zeros(1, Nt);
T_075 = zeros(1, Nt);

for n = 1:Nt
    T_new = T;
    for i = 2:N_x-1
        for j = 2:N_y-1
            T_new(i,j) = T(i,j) + alpha*dt * ( ...
                (T(i+1,j) - 2*T(i,j) + T(i-1,j)) / dx^2 + ...
                (T(i,j+1) - 2*T(i,j) + T(i,j-1)) / dy^2 );
        end
    end
    
    T_new(1,:) = 200;
    
    T_new(end,:) = T(end,:) + (q_right * dx) / k;
   
    T_new(:,1) = T_new(:,2);
    
    T_new(:,end) = (T_new(:,end-1) + (h*dy/k)*T_inf) / (1 + h*dy/k);
    
    T = T_new;
    
    T_02(n) = T(round(0.2/dx), round(0.2/dy));
    T_05(n) = T(round(0.5/dx), round(0.5/dy));
    T_075(n) = T(round(0.75/dx), round(0.75/dy));
    
    if mod(n,10) == 0
        contourf(T', 20, 'LineColor', 'none');
        colorbar;
        title(['Temperature distribution at t = ', num2str(n*dt), ' s']);
        xlabel('x');
        ylabel('y');
        drawnow;
    end
end
figure;
plot(0:dt:total_time-dt, T_02, 'r', 'DisplayName', '(0.2, 0.2)');
hold on;
plot(0:dt:total_time-dt, T_05, 'g', 'DisplayName', '(0.5, 0.5)');
plot(0:dt:total_time-dt, T_075, 'b', 'DisplayName', '(0.75, 0.75)');
xlabel('Time (s)');
ylabel('Temperature (C)');
legend('show');
title('Temperature history at selected points');
grid on;