%% Loading Environment Description
% 
%     This environment is a pendulum linearized around an angle. System is
% parameterized by the gravitational value (g), mass of the pendulum (m),
% friction coefficient (k), and the length of the string/rod (l). The time
% delay in the system arises from delays in a PD controller with specified
% gains kp and kd.
% 
%         /////////////////////////////////////////////////////
%         ------------------------+----------------------------
%                                 |\
%                                   \
%                                 |  \
%                                     \
%                                 |    \
%                                       \
%                                 |----->\_
%                                  theta / \
%                                 |      \_/
% 
%     The system is described as a linear dynamical system with
% 
%         xdot(t) = A*x(t) + Ad*x(t - tau(t)) + B*u(t)
% 
%     where:
%         x = linearized state vector of [delta theta, delta omega]^T 
%         u = control scalar for torque which is an input noise.
%         A  = state transition matrix.
%         Ad = delay contribution. From a feedback law suffering time
%              delays.
%         B  = Control matrix

% System parameters
theta0 = 0;  % [rad]
l = 1;       % [m]
m = 1;       % [kg]
g = 9.81;    % [m/s^2]
k = 0.1;     % [kg*m/s]

% PD controller where error is e = theta - theta0
kp = 0.1;
kd = 0.2;

% State Transition Matrix
A = [              0,  1;
     -g/l*cos(theta0), -k];

% State delay contribution
Ad = [  0,     0;
      -kp,   -kd];

% Disturbance matrix
B = [0; 1/(m*l*l)];

% Initial conditions
E = pi/36.*eye(1);  % Inital bound

% Delay and Noise parameters
tau_m        =  0;  % Minimum time delay
tau_M        =  0.1;  % Maximum time delay
d_m          = -100.0;  % Minimum delay change rate
d_M          =  0.1;  % Maximum delay change rate
mu_scalar    =  0.25;  % Norm bound on derivative of initial condition function
omega_scalar =  0.1;  % Norm bound on input disturbance