%% Loading Environment Description
% 
%     This example is a full 12 degree-of-freedom quadcopter model around
% the hover flight condition. This flight condition is extremely conducive
% to linearization. 
% 
%         

%     The system is described as a linear dynamical system with
% 
%         xdot(t) = A*x(t) + Ad*x(t - tau(t)) + B*u(t)
% 
%     where:
%         x = linearized state vector of [x,y,z,vx,vy,vz,phi,theta,psi,omegax,omegay,omegaz]^T 
%         u = control scalar for torque which is an input noise.
%         A  = state transition matrix.
%         Ad = delay contribution. From a feedback law suffering time
%              delays.
%         B  = Control matrix

% System parameters
g = 9.81;           % gravitational constant [m/s^2]
m = 0.915;          % quadcopter mass [kg]
lx = 0.13;          % rotor distance arm from body x-axis [m]
ly = 0.13;          % rotor distance arm from body y-axis [m]
Jxx = 4.9e-2;       % Inertia around x-axis [kg*m^2]
Jyy = 4.9e-2;       % Inertia around y-axis [kg*m^2]
Jzz = 6.9e-2;       % Inertia around z-axis [kg*m^2]
kdx = 0.397/0.915;  % linear drag coefficient for x axis [Ns/m]
kdy = 0.319/0.915;  % linear drag coefficient for y axis [Ns/m]
kdz = 0.231/0.915;  % linear drag coefficient for z axis [Ns/m]
kT = 2.269e-8;      % Rotor thrust coefficient [N/rpm^2]
kQ = 0.05*kT;        % Rotor torque coefficient [N*m/rpm^2](UNKNOWN)

% Artificial dampening 
kl = 1;     % [1/s]
km = 2;     % [1/s]
kn = 0.01;  % [1/s]

% Nominal States
P = 0;      % Roll rate   [rad/s]
Q = 0;      % Pitch rate  [rad/s]
R = 0;      % Yaw rate    [rad/s]
U = 0;      % X-body velocity [m/s]
V = 0;      % Y-body velocity [m/s] 
W = 0;      % Z-body velocity [m/s]
Phi = 0;    % Roll angle  [rad]
Theta = 0;  % Pitch angle [rad]

% State matrix [-sT; cTsP; cTcP]
A22 = -[kdx, 0, 0; 0, kdy, 0; 0, 0, kdz]/m + ...  % Drag component
      [0, R, -Q; -R, 0, P; Q, -P, 0];  % Rotational frame component
A23 = g*[                   0,          -cos(Theta), 0;...
          cos(Theta)*cos(Phi), -sin(Theta)*sin(Phi), 0;...
         -cos(Theta)*sin(Phi), -sin(Theta)*cos(Phi), 0];
A34 = [0, -W, V; W, 0, -U; -V, U, 0];  % Rotational frame component
A42 = -[kl, 0, 0; 0, km, 0; 0, 0, kn];  % Artificially included for stability.
A44 = [                0, R*(Jyy - Jzz)/Jxx, Q*(Jyy - Jzz)/Jxx;...
       R*(Jzz - Jxx)/Jyy,                 0, P*(Jzz - Jxx)/Jyy;...
       Q*(Jxx - Jyy)/Jzz, P*(Jxx - Jyy)/Jzz,                 0];

A = [zeros(3),   eye(3), zeros(3), zeros(3);...  % d/dt(p)
     zeros(3),      A22,      A23, zeros(3);...  % d/dt(v)
     zeros(3), zeros(3), zeros(3),   eye(3);...  % d/dt(Theta)
     zeros(3),      A42, zeros(3),     A44];     % d/dt(Omega)

% Control Matrix
B = [zeros(3,4);...                % Linear integrator
     zeros(1,4);...                % ddx
     zeros(1,4);...                % ddy
     -kT/m*ones(1,4);...           % ddz
     zeros(3,4);...                % Angular integrator
     kT*lx/Jxx*[1, -1, -1, 1];...  % ddphi
     kT*lx/Jyy*[1, 1, -1, -1];...  % ddtheta
     kQ/Jzz*[-1, 1, -1, 1]];       % ddpsi

% Delay matrix
[K, ~, ~] = lqr(A, B, 2*ones(12), ones(4), zeros(12,4));
% xvec = [x, y, z, vx, vy, vz, phi, theta, psi, omegax, omegay, omegaz]
% kp = 0.5;
% kv = -0.5;
% kap = 0.25;
% kav = 0.05;
% K = 0.25*[1, 1, 1, 1; 1, -1, -1, 1; 1, 1, -1, -1; -1, 1, -1, 1]*...
%     [kp*eye(3), kv*eye(3), kap*eye(3), kav*eye(3); zeros(1,6), kap*[0, 0, 1], kav*[0, 0, 1]];
Ad = -B*K;