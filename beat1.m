% Beta factor
beta = 0.5;  % Adjust as needed (0 < beta <= 1)

 % Apply Boundary Conditions
    T_new(1,:) = 200;  % Left wall: constant temperature (200°C)
    T_new(end,:) = T(end,:) + (q_right * dx) / k;  % Right wall: heat flux
    T_new(:,1) = T_new(:,2);  % Bottom wall: insulated
    T_new(:,end) = (T_new(:,end-1) + (h*dy/k)*T_inf) / (1 + h*dy/k);  % Top wall: convection

% Time-stepping loop using Beta formulation
for n = 1:Nt
    T_new = T;  % Initialize new temperature array
    
    % Update the interior points using the Beta formulation
    for i = 2:N_x-1
        for j = 2:N_y-1
            explicit_part = alpha * dt * ( ...
                (T(i+1,j) - 2*T(i,j) + T(i-1,j)) / dx^2 + ...
                (T(i,j+1) - 2*T(i,j) + T(i,j-1)) / dy^2 );
            T_new(i,j) = beta * (T(i,j) + explicit_part) + (1 - beta) * T(i,j);
        end
    end
    
    % Update the temperature array
    T = T_new;
    
    % Store temperatures at selected points
    T_02(n) = T(round(0.2/dx), round(0.2/dy));
    T_05(n) = T(round(0.5/dx), round(0.5/dy));
    T_075(n) = T(round(0.75/dx), round(0.75/dy));
    
    % Plot temperature distribution every 10 time steps
    if mod(n,10) == 0
        contourf(T', 20, 'LineColor', 'none');
        colorbar;
        title(['Temperature distribution at t = ', num2str(n*dt), ' s']);
        xlabel('x');
        ylabel('y');
        drawnow;
    end
end

% Plot temperature history at selected points
figure;
plot(0:dt:total_time-dt, T_02, 'r', 'DisplayName', '(0.2, 0.2)');
hold on;
plot(0:dt:total_time-dt, T_05, 'g', 'DisplayName', '(0.5, 0.5)');
plot(0:dt:total_time-dt, T_075, 'b', 'DisplayName', '(0.75, 0.75)');
xlabel('Time (s)');
ylabel('Temperature (°C)');
legend('show');
title('Temperature history at selected points');
grid on;
