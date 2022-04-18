% OCE3521 - Wave Solver
% Braidan Duffy
% Due: 04/08/2021

% Inputs
H0 = 1; % m
T = 9; % s
theta0 = 20; % deg
h = 10; % m - design water depth

%% Determine Design Wave Parameters

[L, C, n, Cg, theta, Ks, Kr, H, E, F] = wave_params(H0, T, h, theta0);

% Display
disp("Design Wave Parameters")
fprintf("Design wavelength:       %0.3f m\n", L)
fprintf("Design celerity:         %0.3f m/s\n", C)
fprintf("Design group celerity:   %0.3f m/s\n", Cg)
fprintf("Design shore angle:      %0.3f deg\n", theta)
fprintf("Design wave height:      %0.3f m\n", H)
fprintf("Design wave energy:      %0.3f J/m^2\n", E)
fprintf("Design wave energy flux: %0.3f J/m/s\n\n", F)

%% Determine Design Wave Trajectories and Pressure

[Sx, Sz, u, w, Kp, P, z] = wave_particles(H, T, h);

% Display
disp("Design Wave Trajectories and Pressures")
fprintf("Max design x-displacement:               %0.3f m\n", max(Sx))
fprintf("Max design z-displacement:               %0.3f m\n", max(Sz))
fprintf("Max design horizontal velocity (u):     %0.3f m/s\n", max(u))
fprintf("Max design vertical velocity (w):       %0.3f m/s\n", max(w))
fprintf("Hydro pressure beneath wave (10 partitions): [%0.0f %0.0f %0.0f %0.0f %0.0f %0.0f %0.0f %0.0f %0.0f %0.0f] Pa\n\n", P.')

%% Utility Functions

% Calculate wave parameters for a wave in a specified water depth using
% multiple equations and dispersion
% @param H_0:       deepwater wave height in meters (m)
% @param T:         wave period in seconds (s)
% @param theta_0:   deep water shore angle in degrees (°)
% @param h:         design water depth in meters (m)
% @param rho:       design water density in kg/m^3 (Default: 1025)
% @param g:         gravitational acceleration in m/s^2 (Default: 9.81)

% @return L:        design wavelength in meters (m)
% @return C:        design wave celerity in meters/second (m/s)
% @return n:        design wave group modulation
% @return Cg:       design wave group celerity in meters/second (m/s)
% @return theta:    design wave shore angle in degrees (°)
% @return H:        design wave height in meters (m)
% @return F:        design wave energy flux in Joules/sq. meter (J/m^2)
function [L, C, n, Cg, theta, Ks, Kr, H, E, F] = wave_params(H0, T, h, theta0, rho, g)
    arguments
        % TODO: Argument validation
        H0
        T
        h
        theta0 = 0;    % Deg
        rho = 1025;     % kg/m^3
        g = 9.81;       % m/s^2
    end
    
    [L0, L, ~, ~, ~, kh] = dispersion(T, h, g); % Calculate basic design wave parameters
    
    C0 = L0/T;
    C = L/T;                                % Caluclate wave celerity
    n = (1/2*(1+2*kh/sinh(2*kh)));          % Calculate design group wave modulation
    Cg = n*C;                               % Calculate design group wave celerity
    theta = asind(C0/C*sind(theta0));      % Calculate design shore angle using Snell's law
    Ks = sqrt(L0/(2*T*Cg));                 % Calculate shoaling coefficient
    Kr = sqrt(cosd(theta0)/cosd(theta));   % Calculate refraction coefficient
    H = H0 * Ks * Kr;                      % Calculate design wave height
    E = 1/8 * rho * g * H^2;                % Calculate average wave energy per unit area
    F = E * Cg;                             % Calculate design energy flux
end

% Calculates the particle trajectories, velocities, and pressures beneath a
% given design wave in 10 equal partitions from z=0 to z=-h.

% @param H: the design wave height in meters (m)
% @param T: the design wave period in seconds (s)
% @param h: the design water depth in meters (m)
% @param g: gravitation acceleration in meters/second (m/s) (Default: 9.81)
% @param rho: fluid density in kg/m^3 (Default: 1025)

% @return Sx: Array of maximum particle x-displacements at each water depth partition in m
% @return Sz: Array of maximum particle z-displacements at each water depth partition in m
% @return u: Array of maximum horizontal velocities at each water depth partition in m/s
% @return w: Array of maximum vertical velocities at each water depth partition in m/s
% @return Kp: Array of pressure reduction factors through the water column
% @return P: Array of pressures at each water depth partition in Pa
% @return z: Array of water depths in m
function [Sx, Sz, u, w, Kp, P, z] = wave_particles(H, T, h, g, rho)
    arguments
        % TODO: Argument validation
        H
        T
        h
        g = 9.81; % m/s^2
        rho = 1025; % kg/m^3
    end
    
    [~,~,~,k] = dispersion(T, h, g);    % Calculate wave number with dispersion equation
    sigma = 2*pi/T;                     % Calculate angular frequency
    z = linspace(0, -h, 10); % Generate 10 equal partitions through the water column
    for i=1:length(z)                                           % Iterate through each depth
        Sx(i) = H/2 * cosh(k*(h+z(i)))/sinh(k*h);              % Calculate max x-displacement
        Sz(i) = H/2 * sinh(k*(h+z(i)))/sinh(k*h);               % Calculate max z-displacement
        u(i) = H/2 * g*k/sigma * cosh(k*(h+z(i)))/cosh(k*h);    % Calculate maximum x velocity
        w(i) = H/2 * g*k/sigma * sinh(k*(h+z(i)))/cosh(k*h);    % Calculate maximum z velocity
        Kp(i) = cosh(k*(h+z(i)))/cosh(k*h);                     % Calculate pressure response factor
        P(i) = rho*g*H/2*Kp(i) - rho*g*z(i);                    % Calculate hydro pressure
    end
end

% Calculates the theoretical wavelength using the dispersion equation:
% L = L0 * tanh( 2*pi*h / L )

% @param T: the wave period in seconds (s)
% @param h: the water depth in meters (m)
% @param g: the acceleration due to gravity in m/s^2 (Defaults to 9.792)
% @param tol: the accuracy tolerance of the calculation in meters (Defaults to 0.001)

% @return L_0: the deepwater wavelength in m
% @return L: the design wavelength in m
% @return k0: the deep water wave number
% @return k: the design wavenumber
% @return kh0: the deep water relative depth
% @return kh: the design relative depth 
% @return num_iter: the number of iterations it took to calculate the wavelength
function [L0, L, k0, k, kh0, kh, num_iter] = dispersion(T, h, g, tol)
    arguments
        T
        h
        g = 9.81
        tol = 0.001 % ±1 mm
    end
    
    num_iter = 0;
    L0 = g * T^2 / (2 * pi); % Define the deepwater wavelength
    k0 = 2*pi/L0;
    kh0 = k0*L0/2;
    y = 0; % Instantiate an interim temporary value
    L = L0; % Start the iteration with the deep water wavelength
    
    while(abs(L-y) > tol) % Iterate through dispersion equation until tolerance is reached
        y = L; % Save previous iteration of L
        L = L0 * tanh(2*pi*h / L); % Define next iteration of L
        L = abs((L+y) / 2); % Find average between calculations to get towards tolerance
        num_iter = num_iter + 1; % Increment number of iterations
    end % Return the final L value as the theoretical wavelength
    k = 2*pi/L;                     % Calculate wave number
    kh = k*h;
end