%% Loading Environment Description
% 
%     This environment is the example of the paper "On reachable sets for linear systems with dealy and bounded peak 
%  inputs", E. Fridman, U. Shaked, June 2003.
% 
%     The system is described as a linear dynamical system with
% 
%         xdot(t) = A*x(t) + Ad*x(t - tau(t)) + B*u(t)
% 
%     where:
%         x = state vector 
%         u = control scalar for torque which is an input noise
%         A = state matrix (name in the Fridman_paper: A0)
%         Ad = delay contribution (name in the Fridman_paper: Ad)
%         B  = Control matrix
%         tau_M = max. time delay (name in the Fridman_paper: h)
%         omega_scalar = Bounds on the disturbance vector omega such that
%                        omega(t)^T omega(t) <= omega_scalar^2
%                        (name in the Fridman_paper: wbar) 
%                             wbar = omega_scalar^2

% scalar parameter
rho = 0; % rho < |0.2|

% State Matrix
A = [-2,        0;
      0, -0.9+rho];

% State delay Matrix
Ad = [-1,          0;
      -1, -1+0.5*rho];
 
% Control Matrix
B = [-0.5;
        1];

% Max. time delay
tau_M = 0.7;

% Scalar
omega_scalar = 1;