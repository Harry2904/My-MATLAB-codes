clc;
clear all;
matfiles = dir('*.mat') ;
N = length(matfiles) ;
u_cell = cell(N,1) ;
v_cell = cell(N,1) ;
x_cell = cell(N,1) ;
y_cell = cell(N,1) ;
vort_cell = cell(N,1) ;
for i = 1:N
M = load(matfiles(i).name);
u_cell{i}=M.("u");
v_cell{i}=M.("v");
x_cell{i}=M.("x");
y_cell{i}=M.("y");
vort_cell{i}=M.("vort");
end
u_sum = 0;
v_sum = 0;
vort_sum = 0;
for i = 1:N
u_sum = u_sum + u_cell{i};
v_sum = v_sum + v_cell{i};
vort_sum = vort_sum + vort_cell{i};
end
u = u_sum./N;
v = v_sum./N;
vort = vort_sum./N;
x = x_cell{N};
y = y_cell{N};
x = x - 0.019;
y = y - 0.036;
%% Calculation of Free Stream Velocity %%
U = 0;
for i = 1:63
U = U + u(i,1) + u(i,2) + u(i,3) + u(i,4) + u(i,5);
end
U = U/63/5;
D = 30e-3;
X = x/D;
Y = y/D;
r=D/4;
% Create a vectortheta.
theta=linspace(0,2*pi,200);
% Generate x-coordinates.
x1=r*cos(theta)/D;
% Generate y-coordinate.
y1=r*sin(theta)/D;
%% Velocity vector plot and Streamlines
figure(1)
hold on
q = quiver(X,Y,u/U,v/U);
% Plotting streamlines
n = 63;
z = 1:7:n;
xstart = ones(length(z),1) * min(min(X));
ystart = Y(z,1);
hold on
h=streamline(X,Y,u/U,v/U,xstart,ystart);
set(h,'color','red')
xlabel("x/D")
ylabel("y/D")
legend([q,h(1)],"Non-Dim Velocity Vector","Streamlines")
title('Time Averaged Velocity Vectors and Streamlines')
saveas(gcf,'vector_stream.jpg')
% Set equal scale on axes.
axis('equal');
%% Vorticity Contour %%
figure(2)
hold on
contourf(X,Y,vort)
xlabel("x/D")
ylabel("y/D")
title('Time Averaged Vorticity Contours')
colorbar
% Set equal scale on axes.
axis('equal');
saveas(gcf,'Av_vorticity.jpg')
%% Streamlines %%
figure (3)
n = 63;
xstart = ones(n,1) * min(min(X));
ystart = Y(1:n,1);
hold on
h=streamline(X,Y,u/U,v/U,xstart,ystart);
set(h,'color','Cyan')
xlabel("x/D")
ylabel("y/D")
title('Streamlines')
% Set equal scale on axes.
axis('equal');
%% Calculation of Drag Coefficient %%
Cd = 0;
imax = size(y,1);
jmax = size(y,2);
j = jmax-2;
for i = 1:imax-1
Dy = (y(i+1,j) - y(i,j));
u_star = (u(i+1,j) + u(i,j))/(2*U);
Cd = Cd + (1/D) * u_star * (1 - u_star) * Dy
end
u_star = u(:,j)./U;
xx = 1/D .* u_star .*(1-u_star);
cd1 = trapz(y(:,j),xx)
%% Finding Coefficient of Drag (Using Simpson's 3/8 Rule) Based on Non-Dimensionalized Time-Averaged U-Velocities
% finding step size
dy = diff(y/D);
stepsize = mean(dy,'all');
% setting first grid number and last grid number for numerical integration
imax = size(y,1); % no of rows
jmax = size(y,2); % no of columns
istart = 1;
iend = floor(imax/3)*3 + 1;
excluded_gridpoints = imax-iend;
if excluded_gridpoints == -1
istart = istart + 1;
iend = iend - 2;
end
% numerical integration by Simpson's 3/8 Rule
for integration_xlocation = 1:2
if integration_xlocation == 1
% taking first third column values
vertical_gridline = 3;
end
if integration_xlocation == 2
% taking last third column values
vertical_gridline = jmax-2;
end
integrand(1,:) = u(:,vertical_gridline)/U;
% taking |U*|.U* where |U*| corresponds to Non-Dimensional Mass Flow Rate, and
% U* corresponds to Non-Dimensional U-Velocity vectors with its specified directions
integrand = (abs(integrand)).*integrand;
sum = 0;
for i = istart:3:iend-3
sum = sum + integrand(1,i);
sum = sum + 3*integrand(1,i+1);
sum = sum + 3*integrand(1,i+2);
sum = sum + integrand(1,i+3);
end
integration = (3/8)*stepsize*sum;
if integration_xlocation == 1
ndinitial_momentumrate = integration;
end
if integration_xlocation == 2
ndfinal_momentumrate = integration;
end
end
% Using Cd* = 2(Fnet*)
ndnetforce_xdir = ndfinal_momentumrate-ndinitial_momentumrate;
dragcoefficient = 2*abs(ndnetforce_xdir);
fprintf('Coefficient of Drag (Cd): %.3f',dragcoefficient);
%% Plotting Instantaneous Velocity Vectors %%
Frame = [1 2 3 4 5 10 20 30 40 50];
for f = 1:length(Frame)
figure(f + 3);
q = quiver(X,Y,u_cell{Frame(f)},v_cell{Frame(f)});
n = 63;
z = 1:7:n;
xstart = ones(length(z),1) * min(min(X));
ystart = Y(z,1);
h=streamline(X,Y,u_cell{Frame(f)}/U,v_cell{Frame(f)}/U,xstart,ystart);
set(h,'color','red')
xlabel('x/D')
ylabel('y/D')
legend([q,h(1)],"Non-Dim Velocity Vector","Streamlines")
title("Inst. Non-Dimensional Velocity Vector Plot & Streamlines [Frame "+ Frame(f) + "]")
%disp(Frame(f))
% Set equal scale on axes.
axis('equal');
filename = ['Inst_Vel ' num2str(Frame(f)) '.jpg'];
saveas(gcf,filename)
%imwrite(img,filename,'jpg')
end
%% Plotting Instantaneous Vorticity Contours %%
Frame = [1 2 3 4 5 10 20 30 40 50];
for f = 1:length(Frame)
figure(f + 3 + 5);
contourf(X,Y, vort_cell{Frame(f)})
xlabel("x/D")
ylabel("y/D")
title("Inst. Vorticity Contour [Frame "+ Frame(f) + "]")
colorbar
%legend('Non-Dim Velocity Vector')
% Set equal scale on axes.
axis('equal');
filename = ['Inst_Vort ' num2str(Frame(f)) '.jpg'];
saveas(gcf,filename)
end
close all
figure(50)
hold on
j = [10 50 60 70 80];
for ele = 1:length(j)
%name = 'x/D = ' + num2str(X(1,j(ele)));
plot(u(:,j(ele))./U, Y(:,j(ele)))
Legend{ele} = strcat('x/D = ', num2str(X(1,j(ele))));
end
xlabel("u/U")
ylabel("y/D")
legend(Legend)
title("Variation of Non-Dim u-Velocity wrt Non-Dim Vertical Distance y/D")
saveas(gcf,'nondim_U_Y.jpg')