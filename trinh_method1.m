function [beta0] = trinh_method1(A, Ad, B, ...
                            E,...
                            tau_m, tau_M,...
                            d_m, d_M,...
                            mu_scalar,...
                            omega_scalar...
                            )
% 
% TRINH_METHOD1 Determines a forward reachable set of a linear system with
%               time delays from an initial point x0.
% 
% Usage: beta0 = trinh_method1(A, Ad, B, x0, tau_m, tau_M, d_m, d_M, ...
%                              mu_scalar, omega_scalar)
% 
% Description: Given a linear system with a single time delay described by
%              the differential equation:
% 
%     d/dt(x) = Ax + Ad x(t - tau(t)) + B u + omega(t)
% 
%             the trinh_method1 solves a set of linear equations to find a
%             bounding ball with radius beta0 for the forward
%             reachability set.
% 
% Intputs:
% 
%     A            = nxn state matrix
%     Ad           = nxn state matrix for the delayed states
%     B            = nxp input matrix
%     E            = semi-positive definite matrix used to define the
%                    initial condition based on the initial condition
%                    function phi such that phi(s)^T E phi(s) <= 1 for all
%                    s in [-tau_M, 0].
%     tau_m, tau_M = scalar bounds on the delay such that the inequalities
%                    0 <= tau_m <= tau(t) <= tau_M are observed.
%     d_m, d_M     = scalar bounds on the derivative of the delay with the
%                    inequalities d_m <= d/dt(tau(t)) <= d_M <= 1
%     mu_scalar    = There is an initial condition function phi(s) that is
%                    C^1 continuous on the time interval [-tau_M, 0] and
%                    maps to R^n. The bounds on its derivative is
%                    sup_{s\in[-\tau_M, 0]} phi'(s)^T phi'(s) <= mu_scalar
%     omega_scalar = Bounds on the disturbance vector omega such that
%                    omega(t)^T omega(t) <= omega_scalar^2
% 
% Outputs:
%     beta0 = the radius of the bounding ball for all forward rechable
%             sets.
% 
% Authors:
%     Patrick McNamee
% 
% Date:
%     October 19, 2021
% 
% References:
%     1. "On backwards and forwards reachable sets bounding for perturbed 
%        time-delay systems", H. Trihn, Phan T. Nam, Pubudu N. Pathirana,
%        and H. P. Le, June 2015.
% 
% Requires:
%     Robost Control Toolbox
%
% TODO:
%     1. Write body of function -- PM

%% Input checking
assert(length(size(A))  == 2, "A must be a matrix.")
assert(length(size(Ad)) == 2, "Ad must be a matrix.")
assert(length(size(B))  == 2, "B must be a matrix.")
assert(length(size(E))  == 2, "E must be a matrix.")

n = size(A, 1);
p = size(B, 2);

%% Constants
Gamma1 = ...
    [ei_generator(n,p,1) - ei_generator(n,p,3);...
     sqrt(3)*(ei_generator(n,p,1) + ...
              ei_generator(n,p,3) + ...
              (-2)*ei_generator(n,p,5))...
     ];
Gamma2 = ...
    [ei_generator(n,p,2) - ei_generator(n,p,4);...
     sqrt(3)*(ei_generator(n,p,2) + ...
              ei_generator(n,p,4) + ...
              (-2)*ei_generator(n,p,7));...
     ei_generator(n,p,3) - ei_generator(n,p,2);...
     sqrt(3)*(ei_generator(n,p,3) + ...
              ei_generator(n,p,2) + ...
              (-2)*ei_generator(n,p,6))...
     ];

mu0 = min(eig(E));

%% Linear Matrix Inequality Variable Setup
setlmis([]);  % initialize the linear matrix inequality setup

% Scalar variables
alpha;
beta0 = lmivar(1, [1,0]);
[beta1, ~, sbeta1] = lmivar(1, [1, 0]);
[beta2, ~, sbeta2] = lmivar(1, [1, 0]);
[beta3, ~, sbeta3] = lmivar(1, [1, 0]);
[beta4, ~, sbeta4] = lmivar(1, [1, 0]);
q1 = lmivar(1, [1, 0]);  % Scalar
q2 = lmivar(1, [1, 0]);
q3 = lmivar(1, [1, 0]);
r1 = lmivar(1, [1, 0]);
r2 = lmivar(1, [1, 0]);

% Matrix variables
Q1 = lmivar(1, [n, 1]);  % Symmetric nxn matrix of full elements
Q2 = lmivar(1, [n, 1]);  % nxn
Q3 = lmivar(1, [n, 1]);  % nxn
[R1, ~, sR1] = lmivar(1, [n, 1]);  % nxn
[R2, ~, sR2] = lmivar(1, [n, 1]);  % nxn
P111 = lmivar(1, [2*n, 1]);  % 2n x 2n symmetric
P112 = lmivar(1, [2*n, 1]);  % 2n x 2n symmetric
P12 = lmivar(2, [2*n, 2*n]);  % Full 2n x 2n
P22 = lmivar(1, [2*n, 1]);  % 2n x 2n symmetric
[X, ~, sX] = lmivar(2, [2*n, 2*n]);  % Full 2n x 2n


% Constructed Matrix Variables
beta12 = lmivar(3,[sbeta1*eye(n), zeros(n); zeros(n), sbeta2*eye(n)]);
beta34 = lmivar(3,[sbeta3*eye(n), zeros(n); zeros(n), sbeta4*eye(n)]);
R1tilde = lmivar(3, [sR1, zeros(n); zeros(n) sR1]);
R2tilde = lmivar(3, [sR2, zeros(n); zeros(n) sR2]);


%% Linear Matrix Inequality Constraint Setup
% Eq 7 [[(tau_M - tau_m)P_11^i, P_12];[*, P_22]] < diag{beta1 In, ... beta}
% i = 1
lmiterm([1 1 1 P111],(tau_M-tau_m),1);  % LHS (tau_M - tau_m)P_11^1
lmiterm([1 1 2 P12], 1, 1);             % LHS P_12
lmiterm([1 2 2 P22], 1, 1);             % LHS P_22
lmiterm([-1 1 1 beta12], 1, 1);         % RHS diag{beta1 I_n, beta2 I_n}
lmiterm([-1 2 2 beta34], 1, 1);         % RHS diag{beta3 I_n, beta4 I_n}

% Eq. 7 i = 2
lmiterm([1 1 1 P112],(tau_M-tau_m),1);  % LHS (tau_M - tau_m)P_11^1
lmiterm([1 1 2 P12], 1, 1);             % LHS P_12
lmiterm([1 2 2 P22], 1, 1);             % LHS P_22
lmiterm([-1 1 1 beta12], 1, 1);         % RHS diag{beta1 I_n, beta2 I_n}
lmiterm([-1 2 2 beta34], 1, 1);         % RHS diag{beta3 I_n, beta4 I_n}

% Eq 11.

%% Linear Matrix Inequality Solving
% TODO: use gevp from MATLAB's robust control toolbox.


end

%% Auxiliary functions
function xprime = He(x)
    xprime = (x + x')/2;
end

function ei = ei_generator(n, p, i)
    if i < 8
        ei = [zeros(n, (i-1)*n), eye(n), zeros(n, (7-i)*n + p)];
    else
        ei = [zeros(p, 7*n), eye(p)];
    end
end