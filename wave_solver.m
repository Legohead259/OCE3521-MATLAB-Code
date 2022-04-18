% Wave Solver
% Braidan Duffy
% 3/29/2021

% Program to solve for wave conditions across a profile
clear all; close all;

% Constants
g = 9.81;       % m/s^2
rho_sw = 1025;  % kg/m^3

% Defined wave parameters
m = 1/50;
T = 15;               % s - DW wave period
H_0 = 2;              % m - DW wave height
theta_0 = 30;         % deg - DW wave angle to shore
h = 8;                % m - water depth of interest
wave_freq = 2*pi/T;   % Hz

%% 1) Find L_0 and L
[L_0, L, k_0, k, kh_0, kh, ~] = dispersion(T, h);

%% 2) Find C_0, Cg_0, C, Cg
n_0 = 0.5;
C_0 = L_0/T;
Cg_0 = n_0 * C_0;

n = 0.5*(1+(2*kh/sinh(2*kh)));
C = L/T;
Cg = n * C;

%% 3) Find theta
theta = asind(C/C_0 * sind(theta_0));

%% 4) Use conservation of energy to find H
H = H_0 * sqrt(Cg_0 / Cg) * sqrt(cosd(theta_0)/cosd(theta));