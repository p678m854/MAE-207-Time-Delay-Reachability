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
%     d/dt(x) = Ax + Ad x(t - tau(t)) + B omega(t)
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
Ac = [A, Ad, zeros(n,5*n), B]';
e1 = ei_generator(n, p, 1);
e2 = ei_generator(n, p, 2);
e3 = ei_generator(n, p, 3);
e4 = ei_generator(n, p, 4);
e5 = ei_generator(n, p, 5);
e6 = ei_generator(n, p, 6);
e7 = ei_generator(n, p, 7);
e8 = ei_generator(n, p, 8);
Gamma_1 = ...
    [ei_generator(n,p,1) - ei_generator(n,p,3),...
     sqrt(3)*(ei_generator(n,p,1) + ...
              ei_generator(n,p,3) + ...
              (-2)*ei_generator(n,p,5))...
     ];
Gamma_2 = ...
    [ei_generator(n,p,2) - ei_generator(n,p,4),...
     sqrt(3)*(ei_generator(n,p,2) + ...
              ei_generator(n,p,4) + ...
              (-2)*ei_generator(n,p,7)),...
     ei_generator(n,p,3) - ei_generator(n,p,2),...
     sqrt(3)*(ei_generator(n,p,3) + ...
              ei_generator(n,p,2) + ...
              (-2)*ei_generator(n,p,6))...
     ];

mu0 = min(eig(E));


%% Linear Matrix Inequality Solving
    function [beta0_or_inf] = lmi_portion(alpha,beta0)
%
% Usage: [beta0_or_inf] = lmi_portion(alpha, beta0, beta1, beta2, beta3, beta4, q1, q2, ...
%     q3, r1, r2)
% 
% Description: Internal loop that given 11 scalar values solves an LMI.
% 
% Input:
%     alpha = parameter associated with Lyapunov-Krasovskii stability.
%     beta0 = radius of the bounding ball
% 
% Output:
%     beta0_or_inft = beta0 if solvable LMI with parameters, otherwise
%                     infinite for outerloop optimization
    
    %% Linear Matrix Inequality Variable Setup
    setlmis([]);  % initialize the linear matrix inequality setup
    
    % Scalar constants
    [beta1, ~, beta1s] = lmivar(1, [1, 1]);
    [beta2, ~, beta2s] = lmivar(1, [1, 1]);
    [beta3, ~, beta3s] = lmivar(1, [1, 1]);
    [beta4, ~, beta4s] = lmivar(1, [1, 1]);
    q1 = lmivar(1, [1, 1]);
    q2 = lmivar(1, [1, 1]);
    q3 = lmivar(1, [1, 1]);
    r1 = lmivar(1, [1, 1]);
    r2 = lmivar(1, [1, 1]);

    % Matrix constants
    beta12 = lmivar(3, [beta1s*eye(n), zeros(n); zeros(n), beta2s*eye(n)]);
    beta34 = lmivar(3, [beta3s*eye(n), zeros(n); zeros(n), beta4s*eye(n)]);
    
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
    [~, ~, sX] = lmivar(2, [2*n, 2*n]);  % Full 2n x 2n
    
    % Constructed Matrix Variables
    [R1tilde, ~, ~] = lmivar(3, [sR1, zeros(n); zeros(n) sR1]);
    [~, ~, sR2tilde] = lmivar(3, [sR2, zeros(n); zeros(n) sR2]);
    Theta = lmivar(3, [sR2tilde, sX; sX', sR2tilde]);
    
    %% Linear Matrix Inequality Constraint Setup
    % Eq 6
    lmiterm([1 1 1 Q1], 1, 1);
    lmiterm([-1 1 1 q1], eye(n), 1);
    lmiterm([2 1 1 Q2], 1, 1);
    lmiterm([-2 1 1 q2], eye(n), 1);
    lmiterm([3 1 1 Q3], 1, 1);
    lmiterm([-3 1 1 q3], eye(n), 1);
    lmiterm([4 1 1 R1], 1, 1);
    lmiterm([-4 1 1 r1], eye(n), 1);
    lmiterm([5 1 1 R2], 1, 1);
    lmiterm([-5 1 1 r2], eye(n), 1);
    
    % Eq 7 [[(tau_M - tau_m)P_11^i, P_12];[*, P_22]] < diag{beta1 In, ... beta}
    % i = 1
    lmiterm([6 1 1 P111],(tau_M-tau_m),1);  % LHS (tau_M - tau_m)P_11^1
    lmiterm([6 1 2 P12], 1, 1);             % LHS P_12
    lmiterm([6 2 2 P22], 1, 1);             % LHS P_22
    lmiterm([-6 1 1 0], beta12);         % RHS diag{beta1 I_n, beta2 I_n}
    lmiterm([-6 2 2 0], beta34);         % RHS diag{beta3 I_n, beta4 I_n}
    
    % Eq. 7 i = 2
    lmiterm([7 1 1 P112],(tau_M-tau_m),1);  % LHS (tau_M - tau_m)P_11^1
    lmiterm([7 1 2 P12], 1, 1);             % LHS P_12
    lmiterm([7 2 2 P22], 1, 1);             % LHS P_22
    lmiterm([-7 1 1 0], beta12);            % RHS diag{beta1 I_n, beta2 I_n}
    lmiterm([-7 2 2 0], beta34);            % RHS diag{beta3 I_n, beta4 I_n}
    
    % Eq 8 [I/beta0^2, zeros; zeros]< [[(tau_M - tau_m)P_11^i, P_12];[*, P_22]]
    % i = 1
    lmiterm([-8 1 1 P111], (tau_M - tau_m), 1);
    lmiterm([-8 1 2 P12], 1, 1);
    lmiterm([-8 2 2 P22], 1, 1);
    lmiterm([8 1 1 0], ones(2*n)/(beta0^2));
    
    % i = 2
    lmiterm([-9 1 1 P112], (tau_M - tau_m), 1);
    lmiterm([-9 1 2 P12], 1, 1);
    lmiterm([-9 2 2 P22], 1, 1);
    lmiterm([9 1 1 0], ones(2*n)/(beta0^2));
    
    % Eq 9
    lmiterm([-10 1 1 Theta], 1, 1);
    
    % Eq 10. Sigma(tau, dtau) - alpha/omega^2 e8 e'_8 <0
    %     This one is long but from Lemma 4, one only needs to look at 4
    %     points. However it is long and tedius so I'll put it in a for loop.
    lmit_i = 10;  % current term index for keeping track

    for tau=[tau_m, tau_m]
        for dtau=[d_m, d_M]
            % Update which constrain term this is
            lmit_i = lmit_i+1;
    
            % He operator
            lmiterm([lmit_i 1 1 P111],...
                (tau_M - tau)*[e1, tau_m*e5],...
                [Ac, e1-e3]','s');  % P111 term
            lmiterm([lmit_i 1 1 P112],...
                (tau - tau_m)*[e1, tau_m*e5],...
                [Ac, e1-e3]','s');  % P112 term
            lmiterm([lmit_i 1 1 -P12],...
                (tau - tau_m)*[e1, tau_m*e5],...
                [Ac, e1-e3]','s');  % P12^T term
            lmiterm([lmit_i 1 1 P22],...
                [(tau - tau_m)*e6, (tau_M - tau)*e7],...
                [e3 - (1-dtau)*e2, (1-dtau)*e2 - e4]', 's');  % P12 term
            lmiterm([lmit_i 1 1 P22],...
                [(tau - tau_m)*e6, (tau_M - tau)*e7],...
                [e3 - (1-dtau)*e2, (1-dtau)*e2 - e4]', 's');  % P22 term
    
            % (dt[-P^1_11 + P^2_11] + alpha P) term
            % P11 terms
            lmiterm([lmit_i 1 1 P111], -[e1, tau_m*e5], [e1, tau_m*e5]');
            lmiterm([lmit_i 1 1 P112], [e1, tau_m*e5], [e1, tau_m*e5]');
            lmiterm([lmit_i 1 1 P111], -[e1, tau_m*e5], [e1, tau_m*e5]');
            lmiterm([lmit_i 1 1 P111],...
                alpha*(tau_M - tau)*[e1, tau_m*e5], [e1, tau_m*e5]');
            lmiterm([lmit_i 1 1 P112],...
                alpha*(tau - tau_m)*[e1, tau_m*e5], [e1, tau_m*e5]');
            % P12 terms
            lmiterm([lmit_i 1 1 -P12],...
                alpha*[e1, tau_m*e5],...
                [(tau - tau_m)*e6, (tau_M - tau)*e7]', 's');
            % P22 terms
            lmiterm([lmit_i 1 1 P22],...
                [(tau - tau_m)*e6, (tau_M - tau)*e7],...
                [(tau - tau_m)*e6, (tau_M - tau)*e7]');
    
            % e_1(Q1 + Q2 + Q3)e_1^T
            lmiterm([lmit_i 1 1 Q1], e1, e1');
            lmiterm([lmit_i 1 1 Q2], e1, e1');
            lmiterm([lmit_i 1 1 Q3], e1, e1');
    
            % e3 Q1 e^{-alpha t_m} e3^T
            lmiterm([lmit_i 1 1 Q1], -exp(-alpha*tau_m)*e3, e3');
    
            % e2 Q2 e^{-alpha t_M} e2^T
            lmiterm([lmit_i 1 1 Q2], -(1-dtau)*exp(-alpha*tau_M)*e2, e2');
    
            % e4 Q3 e^{-alpha t_M} e4^T
            lmiterm([lmit_i 1 1 Q3], -exp(-alpha*tau_M)*e4, e4');
    
            % Ac[tau_m^2 R1 + (tau_M - tau_m)^2 R2 Ac
            lmiterm([lmit_i 1 1 R1], (tau_m^2)*Ac, Ac');
            lmiterm([lmit_i 1 1 R2], ((tau_M - tau_m)^2)*Ac, Ac');
    
            % e^{-alpha tau_m} Gamma_1 Rtilde_1 Gamma_1^T
            lmiterm([lmit_i 1 1 R1tilde], -exp(-alpha*tau_m)*Gamma_1, Gamma_1');
    
            % e^{-alpha tau_M} Gamma_2 Theta Gamma_2^T
            lmiterm([lmit_i 1 1 Theta], -exp(-alpha*tau_M)*Gamma_2, Gamma_2');
    
            % Constant term (alpha/omega^2*e8 e8^T)
            lmiterm([lmit_i 1 1 0], alpha/(omega_scalar^2)*e8*e8');
        end
    end

    % Eq. 11
    lmit_i = lmit_i + 1;
    lmiterm([lmit_i 1 1 beta1], 1, mu0*mu0);
    lmiterm([lmit_i 1 1 beta2], tau_m*tau_m, mu0*mu0);
    lmiterm([lmit_i 1 1 beta3], (tau_M - tau_m)^2, mu0*mu0);
    lmiterm([lmit_i 1 1 beta4], (tau_M - tau_m)^2, mu0*mu0);
    lmiterm([lmit_i 1 1 q1], (1 - exp(-alpha*tau_m))/alpha, mu0*mu0);
    lmiterm([lmit_i 1 1 q2], (1 - exp(-alpha*tau_M))/alpha,mu0*mu0);
    lmiterm([lmit_i 1 1 q3], (1 - exp(-alpha*tau_M))/alpha,mu0*mu0);
    lmiterm([lmit_i 1 1 r1], tau_m*(tau_m*alpha + exp(-alpha*tau_m) - 1),...
        ((tau_M - tau_m)*alpha + exp(-alpha*tau_M) - exp(-alpha*tau_m))*...
        (mu_scalar/alpha)^2);
    lmiterm([lmit_i 1 1 r2], tau_m*(tau_M - tau_m),...
        ((tau_M - tau_m)*alpha + exp(-alpha*tau_M) - exp(-alpha*tau_m))*...
        (mu_scalar/alpha)^2);
    lmiterm([-lmit_i 1 1 0], 1);

    % Get the LMI problem and solve
    lmisys = getlmis;
    [tmin, ~] = feasp(lmisys);
    epsilon = 0;  % Acceptable tolerance for LMI solving 
    if tmin < epsilon
        beta0_or_inf = beta0;
    else
        % beta0_or_inf = inf
        beta0_or_inf = beta0*exp(tmin);
    end
    end

%% Outer Loop optimization problem over LMI constraint parameters
% Wrapper functions for objective function and constraints
objf = @(x) lmi_portion(x(1), x(2));

% Set up optimization variable and problem
x = optimvar('x',1,2, 'LowerBound',0);
paramprob = optimproblem(...
    'Objective',fcn2optimexpr(objf, x, 'OutputSize', [1, 1]));

% Initial optimization guess
reA = real(eig(A));
alpha0 = max(reA(reA < 0));   % Largest negative eigenvalue
beta0 = (max(eig(E)))^(0.5);  % radius associated largest orthogonal diagonal.

x0.x = [alpha0, beta0];  % Decay rate (alpha) guess

% Solve problem
[sol,fval,exitflag,output] = solve(paramprob,x0);

beta0 = sol.x(2);

end

%% Auxiliary Functions
function ei = ei_generator(n, p, i)
    if i < 8
        ei = [zeros(n, (i-1)*n), eye(n), zeros(n, (7-i)*n + p)]';
    else
        ei = [zeros(p, 7*n), eye(p)]';
    end
end
