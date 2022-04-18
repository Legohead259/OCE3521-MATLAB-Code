% OCE3521 - Homework 2
% Braidan Duffy
% Due: 03/02/21

%% Problem 2 - DD 3.4
% The equation for the upper moving boundary of a fluid is
%           Zeta(x, t) = 30e^(-(0.02x-0.1t))
% The lower boundary is expressed by
%           Zeta(x, t) = 0
% a) Sketch the boundaries for t=0 s

x_vector = 0:0.1:100; % Generate x values
upper_bound = 30 * exp(-0.02 .* x_vector); % Generate upper boundary

figure(1)
hold on
plot(x_vector, upper_bound)
plot(x_vector, zeros([1 length(x_vector)]))
hold off
title("DD 3.4 Part A - Upper and Lower Boundaries at t=0s")

% b) Discuss the motional characteristics of the upper boundary

%% Problem 5 - DD 3.11
% Draw the streamlines for the following function
% psi = -H/2*g/sigma*sinh(k(h+z))/cosh(kh)*cos(kx-sigma*t)

T = 5; % s
h = 10; % m
H = 2; % m
g = 9.81; % m/s^2
k = calculate_wavenumber_dispersion(2*pi/T, h); % caluclates k iteratively using dispersion equation
psi = -1:1;
Z = 0:0.1:15;

% Plot streamlines
figure(2)
hold on
for N=1:length(psi)
    for z=1:length(Z)
        X(z) = k*acos(-4*pi*cosh(k*h)*psi(N)/(H*g*T*sinh(k*(h+z))));
    end
    plot(X, Z)
end
hold off
title("DD 3.11 - Streamline Plot (Psi = -1, 0, 1)")
legend("Psi = -1", "Psi = 0", "Psi = 1")
ylim([0 5])

%% Problem 8 - DD 3.15
% Develop and iterative technique to solve the dispersion relationship for
% k five sigma and h.
% (sigma)^2 = gktanh(kh) => k = (sigma)^2 / (gtanh(kh))

% See calculate_wavenumber_dispersion.m for the code