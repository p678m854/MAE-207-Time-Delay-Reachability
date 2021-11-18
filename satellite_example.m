%% Loading Environment Description
% 
%     This environment describes a satelitte with two solar panels at 
% angles relative to a common axis of rotation through the body of the 
% satelitte. This system is from "On the reachable set bounding of 
% uncertain dynamic systemswith time-varying delays and disturbances" by
% O.M. Kwon, S.M. Lee, and Ju H. Park
% (https://doi.org/10.1016/j.ins.2011.04.045).
% 
%     +--------------------+         Body           +--------------------+
% /-  |                    |  +-----------------+   |                    |
%-|-- |      Panel 1       +--- Rotational axis ----+       Panel 2      |
% \-> |                    |  +-----------------+   |                    |
%     +--------------------+                        +--------------------+
% 
%     The system is described as a linear dynamical system with
% 
%         xdot(t) = A*x(t) + Ad*x(t - tau(t)) + B*omega(t)
% 
%     where:
%         x  = state vector of [theta1, theta2, omega1, omega2]^T
%         omega  = scalar for torque controller noise to move 1 panel.
%         A  = state transition matrix.
%         Ad = delay contribution. From a feedback law suffering time
%              delays.
%         B  = Control matrix

% State Transition Matrix
A = [      0,         0,    1.0000,         0;
           0,         0,         0,    1.0000;
     -0.0900,    0.0900,   -0.0040,    0.0040;
      0.0900,   -0.0900,    0.0040,   -0.0040];

% State delay contribution
Ad = [      0,         0,         0,         0;
            0,         0,         0,         0;
      -3.3092,   -0.7443,   -2.5909,   -8.0395;
            0,         0,         0,         0];

% Disturbance matrix
B = [0; 0; 1; 0];

% Initial conditions
E = 10*eye(2);  % Inital bound

% Delay and Noise parameters
tau_m        =  0.2;  % Minimum time delay
tau_M        =  0.3;  % Maximum time delay
d_m          = -0.1;  % Minimum delay change rate
d_M          =  0.1;  % Maximum delay change rate
mu_scalar    =  0.1;  % Norm bound on derivative of initial condition function
omega_scalar =  1.0;  % Norm bound on input disturbance
