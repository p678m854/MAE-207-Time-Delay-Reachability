%% Loading Environment Description
% 
%     This example is a simple attitude stabilization of a quadcopter which
% is a 6 degree-of-freedom system. Only the roll, pitch and yaw angles are
% considered. This file has the stabilization point around the hover
% condition where roll and pitch are both 0 and the yaw is arbitrary. The
% quadcopter is the Blackbird quadcopter used by MIT's Aero group. The
% vehicle parameters are found in the "The Blackbird Dataset: A large-scale 
% dataset for UAV perception in aggressive flight" paper at
% https://arxiv.org/abs/1810.01987 mostly. Drag parameters are from
% "Pre-integrated dynamics factors and a dynamical agile visual-inertial 
% dataset for UAV perception" at 
% https://dspace.mit.edu/handle/1721.1/118667
% 
%                                                  _
%          _      _                               / \
%         / \    / \                              \_/
%         \_/    \_/                               |
%            \   /           time            _     |     _
%              X        ~~~~~~~~~~~~~~~~~>  / \____|____/ \
%          _ /   \ _                        \_/    |    \_/
%         / \     / \                              |
%         \_/     \_/                              |
%                                                 / \
%                                                 \_/
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
kT = 2.269e-8;      % Rotor thrust coefficient [N/rpm^2]
kQ = 0.05*kT;        % Rotor torque coefficient [N*m/rpm^2](UNKNOWN)

% Artificial dampening 
kl = 1;     % [1/s]
km = 2;     % [1/s]
kn = 0.01;  % [1/s]

% Nominal States (all zero for hover condition)
% State matrix
A = [zeros(3),   eye(3);...  % d/dt(Theta)
     zeros(3), -[kl, 0, 0; 0, km, 0; 0, 0, kn]];    % d/dt(Omega)

% Control Matrix
B = [zeros(3,4);...                % Angular integrator
     kT*lx/Jxx*[1, -1, -1, 1];...  % ddphi
     kT*lx/Jyy*[1, 1, -1, -1];...  % ddtheta
     kQ/Jzz*[-1, 1, -1, 1]];       % ddpsi

% Delay matrix
simplesys = ss(A(4:end, 4:end), B(4:end,:), eye(3), zeros(3,4));
[K,~,~] = lqi(simplesys, [125*eye(3), zeros(3); zeros(3), 200*eye(3)], eye(4),zeros(6,4));
K = [K(:,4:end), K(:,1:3)];
Ad = B*K;

% Initial conditions
E = pi/12.*eye(1);  % Inital bound

% Delay and Noise parameters
tau_m        =  0.05;  % Minimum time delay
tau_M        =  0.15;  % Maximum time delay
d_m          = -100.0;  % Minimum delay change rate
d_M          =  0.1;   % Maximum delay change rate
mu_scalar    =  0.25;  % Norm bound on derivative of initial condition function
omega_scalar =  0.1;   % Norm bound on input disturbance