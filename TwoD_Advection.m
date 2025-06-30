% Input parameters 
rho = 1000;
cp = 4187;
L1 = 1;
L2 = 1;
u = 1;
v = 1;
imax = 52;
jmax = 52;
epsilon_st = 1e-6;

% Temperature (given condition)
T_ini = 50;
T_wb = 100;
T_sb = 0;
dTc = 100;

% Grid initialization
dz = L1 / (imax - 2);
dy = L2 / (jmax - 2);
z = linspace(0, L1, imax);
y = linspace(0, L2, jmax);

% Initializing the corner points correctly
z_corner = [z(1), z(imax-1)];
y_corner = [y(1), y(jmax-1)];

% Initializing scheme
disp("SELECT THE ADVECTION SCHEME (1/2/3)");
fprintf("1. FOU \n");
fprintf("2. SOU \n");
fprintf("3. QUICK\n");
scheme = input("ENTER THE ADVECTION SCHEME : ");

% Mass flux initialization
mz = rho * u;
my = rho * v;
mz_plus = max(mz, 0);
mz_minus = min(mz, 0);
my_plus = max(my, 0);
my_minus = min(my, 0);

% Initialize temperature field
T = T_ini * ones(imax, jmax); 

% Applying boundary conditions correctly
T(2:imax -1 ,2:jmax -1)=T_ini;
T(:, 1) = T_sb; 
T(1, :) = T_wb; 

% Time-step criteria
dt = 1 / ((abs(u) / dz) + (abs(v) / dy)); % FOU scheme
   if (scheme == 2)
       dt = (0.4)*dt;
      elseif (scheme == 3)
          dt = (4/9)*dt;
   else
   end
unsteadiness_nd = 1;
n = 0;

while unsteadiness_nd >= epsilon_st
    n = n + 1;
    
    % Apply non dirichlet boundary conditions at each step
    T(imax, 1:jmax) = T(imax-1, 1:jmax);
    T(1:imax, jmax) = T(1:imax, jmax-1);
    
    % Store old temperature values
    T_old = T;

%Initializing enthaply-flux temperature
hz_old = zeros(imax, jmax);
hy_old = zeros(imax, jmax);
Q_adv_old = zeros(imax, jmax);
for j = 2:jmax-1
    for i = 1:imax-1
        if (i == 1) || (i == imax-1)
            Te_p = T_old(i, j); 
            Te_m = T_old(i+1, j);
        else
            if i+2 <= imax  
                w = weights(scheme, i);
                [Te_p, Te_m] = T_f(w, T_old(i-1, j), T_old(i, j), T_old(i+1, j), T_old(i+2, j));
            else  
                Te_p = T_old(i, j);
                Te_m = T_old(i+1, j);
            end
        end
        hz_old(i, j) = cp * (mz_plus * Te_p + mz_minus * Te_m);
    end
end
for j = 1:jmax-1
    for i = 2:imax-1
        if (j == 1) || (j == jmax-1)
            Tn_p = T_old(i, j); 
            Tn_m = T_old(i, j+1);
        else
            if j+2 <= jmax  
                w = weights(scheme, i);
                [Tn_p, Tn_m] = T_f(w, T_old(i, j-1), T_old(i, j), T_old(i, j+1), T_old(i, j+2));
            else  
                Tn_p = T_old(i, j);
                Tn_m = T_old(i, j+1);
            end
        end
        hy_old(i, j) = cp * (my_plus * Tn_p + my_minus * Tn_m);
    end
end
for j=2:jmax -1
     for i=2:imax -1
         Q_adv_old(i,j)=((hz_old(i,j)-hz_old(i-1,j))*dy)+(( hy_old(i,j)-hy_old(i,j-1))*dz);
            T(i,j)=T_old(i,j) -(dt/(rho*cp*dz*dy))*Q_adv_old(i,j);
     end
end

 unsteadiness_nd = max(max(abs(T - T_old)));
%fprintf("Time step no. %5d , unsteadiness_nd) = %8.4e\n", n , unsteadiness_nd);
end
T_FOU = zeros(imax, jmax);
T_SOU = zeros(imax, jmax);
T_QUICK = zeros(imax, jmax);

figure;
contourf(z, y, T', 20, 'LineColor', 'none');
colorbar;
xlabel('Length (m)');
ylabel('Height (m)');
title('Temperature Contour Plot');
colormap(jet);

%Plotting for temperature at centerline
figure;
 plot(T(imax/2,:), y, 'r-', 'LineWidth', 2); hold on;
 xlabel('Y-axis (m)');
 ylabel('Temperature (°C)');
 title('Temperature Profile at Vertical Centerline (z=0.5)');
 grid on;