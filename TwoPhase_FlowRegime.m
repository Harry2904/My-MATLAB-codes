mass_flux = input("G: ");
rho_f = 1053.1;
rho_g = 87.271;
heat_flux = 100e3;
D = 0.01;
hf = 287.43e3;
hg = 426.62e3;
H_fg = hg - hf;

m = (4 * heat_flux) / (mass_flux * D * H_fg);

L_2phase = 1 / m;

z = linspace(0,L_2phase,10);

x = m * z;

SM_f = (mass_flux^2 / rho_f) * (1 - x).^2;
SM_g = (mass_flux^2 / rho_g) * (x).^2;

plot(SM_f, SM_g, 'bo-', 'LineWidth', 2, 'MarkerSize', 6)
xlabel('SM_f')
ylabel('SM_g')
title(sprintf('Mass Flux (G) = %.2f', mass_flux))
grid on

data = ones(10, 3);
data(:,1) = z;
data(:,2) = SM_f;
data(:,3) = SM_g;
