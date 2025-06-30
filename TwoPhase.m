clc; clear; close all;

% Given Data
D = 0.037;
L = 3.0;
G_f = 1500;
G_g = 130;         
G = G_f + G_g;
P = 10e5;
T = 25;
S = 1;

% Fluid Properties (Approximate at 25°C, 10 bar)
rho_f = 997;
rho_g = 0.0231;
mu_f = 8.91e-4;
mu_g = 9.08e-5;
sigma = 0.072;
v_f = 0.001003;
v_g = 43.38;
v_fg = 43.378;

% Quality (x)
x = 0.0797;

%Void Fraction (α)
alpha = 1 / (1 + ((1-x) / x) * (rho_g / rho_f) * S);

%Two-Phase Viscosity and Density(Chichitti Rule)
mu_tp = (x / mu_g + (1-x) / mu_f)^(-1);
rho_tp = (alpha * rho_g + (1-alpha) * rho_f);

%Homogeneous Model Pressure Drop
Re_tp = (G * D) / mu_tp;
f_tp = 0.079 * Re_tp^(-0.25); % Blasius equation for turbulent flow

dP_Hom = (2 * f_tp * (G^2) * (v_f + x * v_fg)) / D;

% ---- Lockhart-Martinelli Model ----
Xtt = ((1-x)/x)^0.9 * (rho_g/rho_f)^0.5 * (mu_f/mu_g)^0.1;
phi_2 = 1 + (20 / Xtt) + (1 / Xtt^2);
dP_LM = phi_2 * dP_Hom;

% ---- Chisholm Model ----
C = 21; % Empirical constant for Chisholm model
phi_2C = 1 + (C / Xtt) + (1 / Xtt^2);
dP_Chisholm = phi_2C * dP_Hom;

% ---- Friedel Correlation ----
E = (1 - x)^2 + x^2 * (rho_f / rho_g);
F = (x^0.78 * (1-x)^0.224);
H = (rho_f / rho_g)^0.91 * (mu_g / mu_f)^0.19 * (1 - mu_g / mu_f)^0.7;
phi_2F = E + 3.24 * F * H;
dP_Friedel = phi_2F * dP_Hom;

% ---- Thom Correlation ----
dP_Thom = 0.23 * (G^1.8) * (x^0.8) * (L/D)^1.2 / (rho_f^0.2);

% Display Results
fprintf('\nCorrected Pressure Drop Results (Pa):\n');
fprintf('Homogeneous Model: %.2f Pa\n', dP_Hom);
fprintf('Lockhart-Martinelli: %.2f Pa\n', dP_LM);
fprintf('Chisholm Model: %.2f Pa\n', dP_Chisholm);
fprintf('Friedel Correlation: %.2f Pa\n', dP_Friedel);
fprintf('Thom Correlation: %.2f Pa\n', dP_Thom);
