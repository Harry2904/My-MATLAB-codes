wavelength1=1030e-9; % wavelength of laser in m
wavelength2=515e-9; % wavelength of laser in m
focal_len= linspace(0, 100e-3, 100); %focal length of lens used in m
M2=1.07; %M-square value of the laser
D=0.00926; % diameter of laser falling on lens in m
DOF_1=(2.56*focal_len.^2*M2*wavelength1)/(pi*D^2);
DOF_1=DOF_1*1e6; % converting to um
DOF_2=(2.56*focal_len.^2*M2*wavelength2)/(pi*D^2);
DOF_2=DOF_2*1e6; % converting to um
plot(focal_len*1000,DOF_1,'r-',"linewidth",2);
hold on
plot(focal_len*1000,DOF_2,'b-',"linewidth",2);
ax = gca;
ax.FontSize = 28;
set(gca, 'fontname', 'Times New Roman');
xlabel(' focal length (mm)');
ylabel('DOF (\mum)');
grid on;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickLength = [0.02, 0.02];
ax.LineWidth = 1;
box on
legend('\lambda = 1030 nm','\lambda = 515 nm');
legend box off