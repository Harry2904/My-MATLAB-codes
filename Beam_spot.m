close all
wavelength=1030e-9; % wavelength of laser in m
rad_laser_exit=(3.65e-3)/2; %diameter at laser exit
focal_len= input('Enter the focal length = ');
M2=1.07; %M-square value of the laser
z1=174.5e-2;
W1 = (rad_laser_exit) * sqrt(1 + ((M2 * wavelength * z1) ./ (pi * (rad_laser_exit)^2)).^2);

%W1 is the beam radius at beam expander inlet
W2=2.5*W1;
z2=(8+9)*1e-2;
W3= (W2) * sqrt(1 + ((M2 * wavelength * z2) ./ (pi * (W2)^2)).^2);
D=W3*2;

spot_dia_foc=(4*focal_len*M2*wavelength)/(pi*D);
fprintf('Spot diameter at focus: %.3f microns\n', spot_dia_foc * 1e6);

DOF=(0.64*pi*(spot_dia_foc/2)^2)/(M2*wavelength);
fprintf('Depth of focus: %.3f microns\n', DOF * 1e6);
hold off;