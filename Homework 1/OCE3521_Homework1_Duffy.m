% OCE3521 - Homework 1
% Braidan Duffy
% Due: 02/11/21

%% Problem 2 - DD2.3
% Sketch two streamlines for t=T/8

% psi(t=T/8) = -cos(pi/4) * (3z+5x)

x = -10:0.1:10;
z_0 = -5/3 .* x; % Z function when psi = 0
z_1 = -1/3 * (5.*x+sec(pi/4)); % z function when psi = 1

% Plotting
figure(1)
title("DD2.3 - Sketch of two streamlines (Psi=0 and Psi=1)")
hold on
plot(x, z_0)
plot(x, z_1)
xlabel("X-Values")
ylabel("Z-Values")
legend("Psi=0", "Psi=1")
hold off

%% Problem 5 - DD2.8
% Sketch the two streamlines through (1, 1) and (1, 2)
x = -10:0.1:10;

z1 = atan(1) ./ x;
z2 = atan(0.5) ./ x;

% Plotting
figure(2)
title("DD2.8 - Sketch the stream lines through (1,1) and (1,2)")
hold on
plot(x, z1)
plot(x, z2)
xlabel("X-Values")
ylabel("Z-Values")
legend("Psi=0", "Psi=1")
hold off

%% Problem 6 - DD2.10
% Sketch the streamlines for psi=0 and psi=6A

x = 0:0.1:10; % Domain x = [0, 10]
z_2 = 2 ./ x .^ 2; % Implicit z function for 

% Plotting
figure(3)
title("DD2.10 - Sketch the streamlines for Psi=0 and Psi=6A")
hold on
plot(x, zeros(1, length(x))) % first function of z is implicitly 0
plot(x, z_2)
xlabel("X-Values")
ylabel("Z-Values")
legend("Psi=0", "Psi=6A")
hold off