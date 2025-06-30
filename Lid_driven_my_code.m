%Input Parameters
rho_o=1;
L2=1;
ratio=input('Enter the ratio = ');
L1=ratio*L2;
imax=input('Enter the imax value = ');
jmax=input('Enter the jmax value = ');
U_o=1; %Lid Velocity
Re=100;
mu=1/Re;
nu=mu/rho_o;

%Unifrom grid generation
delta_x = L1 / (imax - 2);
delta_y = L2 / (jmax - 2);
%cell center arrays(zc,yc)
xc=zeros(1,imax);
xc(1,2)=delta_x/2;
xc(1,imax-1)=L1-(delta_x/2);
xc(1,imax)=L1;
for i=3:imax-2
    xc(1,i)=xc(1,2)+(i-2)*delta_x;
end

yc=zeros(1,jmax);
yc(1,2)=delta_y/2;
yc(1,jmax-1)=L2-(delta_y/2);
yc(1,jmax)=L2;
for i=3:jmax-2
    yc(1,i)=yc(1,2)+(i-2)*delta_y;
end

%Arrays to store distance between cell centres
del_x=zeros(1,imax-1); 
del_y=zeros(1,jmax-1);
for i=1:imax-1
    del_x(1,i)=xc(1,i+1)-xc(1,i);
end
for i=1:jmax-1
    del_y(1,i)=yc(1,i+1)-yc(1,i);
end

%Initialising u-Velocity
[u,u_old,u_star,u_dash]=deal(zeros(imax-1,jmax),zeros(imax-1,jmax),zeros(imax-1,jmax),zeros(imax-1,jmax));

%Initialising v-Velocity
[v,v_old,v_star,v_dash]=deal(zeros(imax,jmax-1),zeros(imax,jmax-1),zeros(imax,jmax-1),zeros(imax,jmax-1));

%Initialising Pressure
[p,p_old,p_dash,po_dash,p_star,po_star,ap]=deal(zeros(imax,jmax),zeros(imax,jmax),zeros(imax,jmax),zeros(imax,jmax),zeros(imax,jmax),zeros(imax,jmax),zeros(imax,jmax));

%Initialising fluxes for mass conservation
[mx_star,mx_staro,my_star,my_staro]=deal(zeros(imax-1,jmax),zeros(imax-1,jmax),zeros(imax,jmax-1),zeros(imax-1,jmax));
[mx_dash,my_dash]=deal(zeros(imax-1,jmax),zeros(imax,jmax-1));

%Initialising source term for mass conservation
[Sm_star,Sm_staro,Sm_dash]=deal(zeros(imax,jmax),zeros(imax,jmax),zeros(imax,jmax));

%Initialising fluxes for x-momentum eqn (muxp=Mux+,muxm=Mux-,muyp=Muy+,muym=Muy-,a=advection,d=diffusion)
[muxp,muxm,muyp,muym,aux,auy,dux,duy]=deal(zeros(imax-1,jmax),zeros(imax-1,jmax),zeros(imax-1,jmax),zeros(imax-1,jmax),zeros(imax-1,jmax),zeros(imax-1,jmax),zeros(imax-1,jmax),zeros(imax-1,jmax));

%Initialising fluxes for y-momentum eqn(mvxp=Mvx+,mvxm=Mvx-,mvyp=Mvy+,mvym=Mvy-,a=advection,d=diffusion)
[mvxp,mvxm,mvyp,mvym,avx,avy,dvx,dvy]=deal(zeros(imax,jmax-1),zeros(imax,jmax-1),zeros(imax,jmax-1),zeros(imax,jmax-1),zeros(imax,jmax-1),zeros(imax,jmax-1),zeros(imax,jmax-1),zeros(imax,jmax-1));

%Initial Conditons for velocities and pressure
u(:,:)=0;
v(:,:)=0;
p(:,:)=0;
p_star(:,:)=0; %predicted pressure
u_dash=u;
v_dash=v;

%Dirichlet BC for velocities and pressure [p_dash=Non-Dirichlet]
u(1,:)=0;v(1,:)=0;
u(imax-1,:)=0;v(imax,:)=0;
u(:,1)=0;v(:,1)=0;
u(:,jmax)=U_o;v(:,jmax-1)=0;

%Grid Fourier Number Criteria
delta_Tdiff=(0.5/nu)*(1/delta_x^2+1/delta_y^2)^-1;

%Loop Variables
unsteadiness_nd = 1;
n_unsteadiness=0;

while unsteadiness_nd > 1e-3
    u_old=u;v_old=v;p_old=p;
    
    % CFL Criteria
    delta_Tadv = (max(max(abs(u_old)))/delta_x+max(max(abs(v_old)))/delta_y)^-1;
    delta_T = 0.5*min(delta_Tdiff, delta_Tadv);
    
    % MASS FLUXES FOR X & Y-MOMENTUM EQUATION
    % Computation of x-mass flux for x-Momentum Conservation Law
    i = 1;
    while i <= imax - 2
        j = 2;
        while j <= jmax - 1
            if (u_old(i+1,j) + u_old(i,j)) > 0
                muxp(i,j) = 0.5 * rho_o * (u_old(i+1,j) + u_old(i,j));
                muxm(i,j) = 0;
            else
                muxp(i,j) = 0;
                muxm(i,j) = 0.5 * rho_o * (u_old(i+1,j) + u_old(i,j));
            end
            j = j + 1;
        end
        i = i + 1;
    end
    
     %Computation of y-mass flux for x-Momentum Conservation Law
    i = 2;
    while i <= imax - 2
        j = 1;
        while j <= jmax - 1
            if (v_old(i,j) + v_old(i+1,j)) > 0
                muyp(i,j) = 0.5 * rho_o * (v_old(i,j) + v_old(i+1,j));
                muym(i,j) = 0;
            else
                muyp(i,j) = 0;
                muym(i,j) = 0.5 * rho_o * (v_old(i,j) + v_old(i+1,j));
            end
            j = j + 1;
        end
        i = i + 1;
    end
    
     %Computation of x-mass flux for y-Momentum Conservation Law
    i = 1;
    while i <= imax - 1
        j = 2;
        while j <= jmax - 2
            if (u_old(i,j) + u_old(i,j+1)) > 0
                mvxp(i,j) = 0.5 * rho_o * (u_old(i,j) + u_old(i,j+1));
                mvxm(i,j) = 0;
            else
                mvxp(i,j) = 0;
                mvxm(i,j) = 0.5 * rho_o * (u_old(i,j) + u_old(i,j+1));
            end
            j = j + 1;
        end
        i = i + 1;
    end
        
    %Computation of y-mass flux for y-Momentum Conservation Law
    i = 2;
    while i <= imax - 1
        j = 1;
        while j <= jmax - 2
            if (v_old(i,j+1) + v_old(i,j)) > 0
                mvyp(i,j) = 0.5 * rho_o * (v_old(i,j+1) + v_old(i,j));
                mvym(i,j) = 0;
            else
                mvyp(i,j) = 0;
                mvym(i,j) = 0.5 * rho_o * (v_old(i,j+1) + v_old(i,j));
            end
            j = j + 1;
        end
        i = i + 1;
    end
      
    % ADVECTION & DIFFUSION FLUXES FOR X & Y-MOMENTUM EQUATION
    %Computation of x-advection & diffusion flux for x-momentum eqn
   i = 1;
while i <= imax - 2
    j = 2;
    while j <= jmax - 1
        aux(i,j) = muxp(i,j) * u_old(i,j) + muxm(i,j) * u_old(i+1,j);
        dux(i,j) = mu * (u_old(i+1,j) - u_old(i,j)) / del_x(1,i);
        j = j + 1;
    end
    i = i + 1;
end
    
    %Computation of y-advection & diffusion flux for x-momentum eqn
   i = 2;
while i <= imax - 2
    j = 1;
    while j <= jmax - 1
        auy(i,j) = muyp(i,j) * u_old(i,j) + muym(i,j) * u_old(i,j+1);
        duy(i,j) = mu * (u_old(i,j+1) - u_old(i,j)) / del_y(1,j);
        j = j + 1;
    end
    i = i + 1;
end
    
    %Computation of x-advection & diffusion flux for y-momentum eqn
    i = 1;
while i <= imax - 1
    j = 2;
    while j <= jmax - 2
        avx(i,j) = mvxp(i,j) * v_old(i,j) + mvxm(i,j) * v_old(i+1,j);
        dvx(i,j) = mu * (v_old(i+1,j) - v_old(i,j)) / del_x(1,i);
        j = j + 1;
    end
    i = i + 1;
end
    
     %Computation of y-advection & diffusion flux for x-momentum eqn
     i = 2;
while i <= imax - 1
    j = 1;
    while j <= jmax - 2
        avy(i,j) = mvyp(i,j) * v_old(i,j) + mvym(i,j) * v_old(i,j+1);
        dvy(i,j) = mu * (v_old(i,j+1) - v_old(i,j)) / del_y(1,j);
        j = j + 1;
    end
    i = i + 1;
end 
    
     %PREDICTOR STEP for U-velocity
      i = 2;
    while i <= imax - 2
        j = 2;
        while j <= jmax - 1
            u_star(i,j) = u_old(i,j) + (delta_T / (rho_o * delta_y * delta_x)) * (((dux(i,j) - dux(i-1,j)) * delta_y + (duy(i,j) - duy(i,j-1)) * delta_x) - ((aux(i,j) - aux(i-1,j)) * delta_y + (auy(i,j) - auy(i,j-1)) * delta_x) + (p_old(i,j) - p_old(i+1,j)) * delta_y);
            j = j + 1;
        end
        i = i + 1;
    end
    
    u_star(:, jmax) = U_o;
    
    
     %PREDICTOR STEP for V-velocity
     i = 2;
while i <= imax - 1
    j = 2;
    while j <= jmax - 2
        v_star(i,j) = v_old(i,j) + (delta_T / (rho_o * del_y(1,j) * delta_x)) * ...
                     (((dvx(i,j) - dvx(i-1,j)) * del_y(1,j) + (dvy(i,j) - dvy(i,j-1)) * delta_x) - ...
                      ((avx(i,j) - avx(i-1,j)) * del_y(1,j) + (avy(i,j) - avy(i,j-1)) * delta_x) + ...
                      (p_old(i,j) - p_old(i,j+1)) * delta_x);
        j = j + 1;
    end
    i = i + 1;
end
    % PREDICTION OF MASS SOURCE TERM
    %Predicting mass flux
    mx_star=rho_o*u_star;
    my_star=rho_o*v_star;
   
    %Predicting mass source term
    i = 2;
while i <= imax - 1
    j = 2;
    while j <= jmax - 1
        Sm_star(i,j) = (mx_star(i,j) - mx_star(i-1,j)) * delta_y + (my_star(i,j) - my_star(i,j-1)) * delta_x;
        j = j + 1;
    end
    i = i + 1;
end
    %%
    if max(max(abs(Sm_star))) < 1e-8
        u=u_star;
        v=v_star;
        p=p_old;
    else
        p_dash(:,:)=0;
        p_star=p_old;
                
        %% Computing Pressure Correction 
        Error_source = 1; 
        n_pressure=0;
        while Error_source > 1e-8
            
            po_dash=p_dash; %pod=P'(N):Pressure Correction from previous GS iteration
            po_star=p_star; 
            mx_staro=mx_star;
            my_staro=my_star;

            %% Pressure Correction Pp'(N+1)
            %Computing Sm*(N)
            for i=2:imax-1
                for j=2:jmax-1
                    Sm_staro(i,j)=(mx_staro(i,j)-mx_staro(i-1,j))*delta_y+(my_staro(i,j)-my_staro(i,j-1))*delta_x;
                end
            end

            %Computing ap
            for i=2:imax-1
                for j=2:jmax-1
                    ap(i,j)=delta_T*((1/del_x(1,i)+1/del_x(1,i-1))*delta_y+(1/del_y(1,j)+1/del_y(1,j-1))*delta_x);
                end
            end
            
            %Computing Pressure Correction Pp'(N+1)
            for i=2:imax-1
                for j=2:jmax-1
                    Sm_dash(i,j)=-delta_T*delta_y*((po_dash(i+1,j)-po_dash(i,j))/del_x(1,i)-(po_dash(i,j)-p_dash(i-1,j))/del_x(1,i-1))-delta_T*delta_x*((po_dash(i,j+1)-po_dash(i,j))/del_y(1,j)-(po_dash(i,j)-p_dash(i,j-1))/del_y(1,j-1));
                    p_dash(i,j)=po_dash(i,j)-(Sm_staro(i,j)+Sm_dash(i,j))/ap(i,j);
                end
            end
            %Non-Dirichlet BC for Pressure Correction
            p_dash(1,:)=p_dash(2,:);
            p_dash(imax,:)=p_dash(imax-1,:);
            p_dash(:,1)=p_dash(:,2);
            p_dash(:,jmax)=p_dash(:,jmax-1);
            
            %% Mass-Flux Correction
            for i=1:imax-1
                for j=2:jmax-1
                    mx_dash(i,j)=-delta_T*(p_dash(i+1,j)-p_dash(i,j))/del_x(1,i);
                end
            end
             for i=2:imax-1
                for j=1:jmax-1
                    my_dash(i,j)=-delta_T*(p_dash(i,j+1)-p_dash(i,j))/del_y(1,j);
                end
             end

             %% Updating Predicted Mass Flux & Velocity
             mx_star=mx_staro+mx_dash;
             my_star=my_staro+my_dash;
             u_star=mx_star/rho_o;
             v_star=my_star/rho_o;
             
             %% New Iterative Values of Mass-source & Pressure
             for i=2:imax-1
                 for j=2:jmax-1
                     Sm_star(i,j)=(mx_star(i,j)-mx_star(i-1,j))*delta_y+(my_star(i,j)-my_star(i,j-1))*delta_x;
                 end
             end
             Error_source=max(max(abs(Sm_star)));
             p_star=po_star+p_dash;          
             n_pressure=n_pressure+1;
             
        end
        
        p=p_star;
        u=u_star;
        v=v_star;
    end 
    
    Cn=U_o*delta_T/L1;
    error_U=max(max(abs(u-u_old)))/Cn; 
    error_V=max(max(abs(v-v_old)))/Cn;
    unsteadiness_nd = max(error_U,error_V);
    
    n_unsteadiness = n_unsteadiness + 1;
    storage(1, n_unsteadiness)= n_pressure;
   
end

%% Velocity at Grid Points
up=zeros(imax,jmax);
up(1,:)=u(1,:);up(imax,:)=u(imax-1,:);
for i=2:imax-1
    up(i,:)=0.5*(u(i,:)+u(i-1,:));
end

vp=zeros(imax,jmax);
vp(:,1)=v(:,1);vp(:,jmax)=v(:,jmax-1);
for i=2:jmax-1
    vp(:,i)=0.5*(v(:,i)+v(:,i-1));
end


%% SELECTING GHIA DATA FILES
% Open file selection dialog for multiple files
[files, path] = uigetfile('*.xlsx', 'Select Two Data Files', 'MultiSelect', 'on');

% Check if files are selected
if isequal(files, 0)
    disp('No files selected. Exiting...');
    return;
end

% Ensure files are stored as a cell array (if only one file is selected)
if ischar(files)
    files = {files};  % Convert single file string to cell array
end

% Check if exactly two files are selected
if length(files) ~= 2
    disp('Please select exactly two files.');
    return;
end

% Read and store both files in separate tables
filePath1 = fullfile(path, files{1});
filePath2 = fullfile(path, files{2});

disp(['Selected File 1: ', filePath1]);
disp(['Selected File 2: ', filePath2]);

% Read tables from both files
data1 = readtable(filePath1, 'PreserveVariableNames', true);
data2 = readtable(filePath2, 'PreserveVariableNames', true);

% Display first few rows of each file
disp('First few rows of File 1:');
disp(head(data1));

disp('First few rows of File 2:');
disp(head(data2));

 %% U-centreline
figure(1)
plot(0.5*(up(floor((imax+1)/2),:)+up(floor((imax+2)/2),:)),yc,'LineWidth', 1);
title('Variation of U-velocity along vertical centreline')
xlabel('U-velocity')
ylabel('Y')
hold on;

%% V-centreline
figure(2)
plot(xc/ratio,0.5*(vp(:,floor((jmax+1)/2))+vp(:,floor((jmax+2)/2))),'LineWidth', 1);
title('Variation of V-velocity along horizontal centreline')
xlabel('X')
ylabel('V-velocity')
hold on;

%% Streamline Function
streamline=zeros(imax,jmax);
streamline(:,1)=0;
for i=2:imax-1
    for j=2:jmax-1
        streamline(i,j)=streamline(i,j-1)+u(i,j)*delta_y;
    end
 end

%% Pressure contour
figure(3)
contourf(xc,yc,transpose(p),25);
colormap(jet)
title('Pressure Contour')
xlabel('X')
ylabel('Y')

figure(4)
contourf(xc,yc,transpose(streamline),(min(min(streamline)):0.01:max(max(streamline))));
colormap(jet)
title('Streamline Contour')
xlabel('X')
ylabel('Y')

figure(5)
contourf(xc,yc,transpose((up.^2+vp.^2).^0.5));
colormap(jet)
title('U-velocity Contour')
xlabel('X')
ylabel('Y')

%% Velocity Vector
figure(6)
quiver(xc,yc,transpose(up),transpose(vp)) 
title('Velocity Vector Plot')
xlabel('X')
ylabel('Y')
hold on;