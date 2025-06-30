clc;
clear all;
close all;

%% Robust file loading section
matfiles = dir('*.mat');
if isempty(matfiles)
    error('No .mat files found in the directory');
end

N = length(matfiles);
u_cell = cell(N,1);
v_cell = cell(N,1);
x_cell = cell(N,1);
y_cell = cell(N,1);
vort_cell = cell(N,1);

% Load first file to get reference dimensions
first_file = load(matfiles(1).name);
ref_size = size(first_file.u);

% Load all files with error checking
for i = 1:N
    try
        M = load(matfiles(i).name);
        
        % Check for required fields
        if ~isfield(M, 'u') || ~isfield(M, 'v') || ~isfield(M, 'x') || ~isfield(M, 'y') || ~isfield(M, 'vort')
            error('Missing required fields in file %s', matfiles(i).name);
        end
        
        % Verify array sizes match first file
        if ~isequal(size(M.u), ref_size)
            error('Array size mismatch in file %s', matfiles(i).name);
        end
        
        % Store data (preserving your original variable names)
        u_cell{i} = M.u;
        v_cell{i} = M.v;
        x_cell{i} = M.x;
        y_cell{i} = M.y;
        vort_cell{i} = M.vort;
        
    catch ME
        warning('Skipping file %s: %s', matfiles(i).name, ME.message);
        continue;
    end
end

% Remove any empty cells from failed loads
valid_files = ~cellfun(@isempty, u_cell);
u_cell = u_cell(valid_files);
v_cell = v_cell(valid_files);
x_cell = x_cell(valid_files);
y_cell = y_cell(valid_files);
vort_cell = vort_cell(valid_files);
N = sum(valid_files);

if N == 0
    error('No valid files were loaded');
end

%% Averaging (original logic preserved exactly)
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
x = x_cell{N}; % Using last valid file's coordinates
y = y_cell{N};

% Your original coordinate adjustment
x = x - 0.019;
y = y - 0.036;

%% Free Stream Velocity Calculation (made more robust)
[rows, cols] = size(u);
if rows < 63 || cols < 5
    warning('Using adaptive free stream calculation due to small array size');
    U = mean(u(:,1:min(5,end)), 'all'); % Safe alternative
else
    % Your original calculation
    U = 0;
    for i = 1:63
        U = U + u(i,1) + u(i,2) + u(i,3) + u(i,4) + u(i,5);
    end
    U = U/(63*5);
end

D = 25.4e-3;
X = x/D;
Y = y/D;
r = D/4;

% Create vector theta
theta = linspace(0,2*pi,200);
x1 = r*cos(theta)/D;
y1 = r*sin(theta)/D;

%% Plot 1: Velocity vectors and streamlines (original logic preserved)
figure(1)
clf; hold on;
q = quiver(X,Y,u/U,v/U);

n = 63;
z = 1:7:n;
xstart = ones(length(z),1) * min(X(:));
ystart = Y(z,1);
h = streamline(X,Y,u/U,v/U,xstart,ystart);
set(h,'color','red')

xlabel("x/D")
ylabel("y/D")
legend([q,h(1)],"Non-Dim Velocity Vector","Streamlines")
title('Time Averaged Velocity Vectors and Streamlines')
axis equal;
saveas(gcf,'vector_stream.jpg')

%% Plot 2: Vorticity contour (original logic preserved)
figure(2)
clf; hold on;
contourf(X,Y,vort)
xlabel("x/D")
ylabel("y/D")
title('Time Averaged Vorticity Contours')
colorbar
axis equal;
saveas(gcf,'Av_vorticity.jpg')

%% Plot 3: Streamlines only (original logic preserved)
figure(3)
clf; hold on;
n = 63;
xstart = ones(n,1) * min(X(:));
ystart = Y(1:n,1);
h = streamline(X,Y,u/U,v/U,xstart,ystart);
set(h,'color','Cyan')
xlabel("x/D")
ylabel("y/D")
title('Streamlines')
axis equal;

%% Drag Coefficient Calculation (original logic with safety checks)
% Your original trapezoidal method
Cd = 0;
imax = size(y,1);
jmax = size(y,2);
j = min(jmax-2, jmax); % Ensure we don't exceed array bounds

for i = 1:imax-1
    Dy = (y(i+1,j) - y(i,j));
    u_star = (u(i+1,j) + u(i,j))/(2*U);
    Cd = Cd + (1/D) * u_star * (1 - u_star) * Dy;
end

u_star = u(:,j)./U;
xx = 1/D .* u_star .*(1-u_star);
cd1 = trapz(y(:,j),xx);

%% Simpson's 3/8 Rule (original logic preserved exactly)
dy = diff(y/D);
stepsize = mean(dy,'all');

imax = size(y,1);
jmax = size(y,2);
istart = 1;
iend = floor(imax/3)*3 + 1;
excluded_gridpoints = imax-iend;

if excluded_gridpoints == -1
    istart = istart + 1;
    iend = iend - 2;
end

for integration_xlocation = 1:2
    if integration_xlocation == 1
        vertical_gridline = 3;
    end
    if integration_xlocation == 2
        vertical_gridline = min(jmax-2, jmax); % Safety check
    end
    
    integrand(1,:) = u(:,vertical_gridline)/U;
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

ndnetforce_xdir = ndfinal_momentumrate-ndinitial_momentumrate;
dragcoefficient = 2*abs(ndnetforce_xdir);
fprintf('Coefficient of Drag (Cd): %.3f\n',dragcoefficient);

%% Instantaneous Velocity Plots (with bounds checking)
Frame = [1 2 3 4 5 10 20 30 40 50];
Frame = Frame(Frame <= N); % Only use frames that exist

for f = 1:length(Frame)
    fig_num = f + 3;
    figure(fig_num);
    clf; hold on;
    
    q = quiver(X,Y,u_cell{Frame(f)},v_cell{Frame(f)});
    n = 63;
    z = 1:7:n;
    xstart = ones(length(z),1) * min(X(:));
    ystart = Y(z,1);
    h = streamline(X,Y,u_cell{Frame(f)}/U,v_cell{Frame(f)}/U,xstart,ystart);
    set(h,'color','red')
    
    xlabel('x/D')
    ylabel('y/D')
    legend([q,h(1)],"Non-Dim Velocity Vector","Streamlines")
    title(sprintf('Inst. Non-Dim Velocity Vector Plot & Streamlines [Frame %d]', Frame(f)))
    axis equal;
    
    filename = sprintf('Inst_Vel %d.jpg', Frame(f));
    saveas(gcf,filename)
end

%% Instantaneous Vorticity Plots (with bounds checking)
for f = 1:length(Frame)
    fig_num = f + 3 + 5;
    figure(fig_num);
    clf; hold on;
    
    contourf(X,Y, vort_cell{Frame(f)})
    xlabel("x/D")
    ylabel("y/D")
    title(sprintf('Inst. Vorticity Contour [Frame %d]', Frame(f)))
    colorbar
    axis equal;
    
    filename = sprintf('Inst_Vort %d.jpg', Frame(f));
    saveas(gcf,filename)
end

%% Final Plot: Velocity Profiles (original logic preserved)
figure(50)
clf; hold on;
j = [10 50 60 70 80];
j = j(j <= size(X,2)); % Ensure columns exist

Legend = cell(1,length(j));
for ele = 1:length(j)
    plot(u(:,j(ele))./U, Y(:,j(ele)))
    Legend{ele} = sprintf('x/D = %.2f', X(1,j(ele)));
end

xlabel("u/U")
ylabel("y/D")
legend(Legend)
title("Variation of Non-Dim u-Velocity wrt Non-Dim Vertical Distance y/D")
saveas(gcf,'nondim_U_Y.jpg')

close all;