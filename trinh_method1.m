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
Ac = [A, Ad, zeros(n,5*n), B]';
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


%% Linear Matrix Inequality Solving
    function [beta0_or_inf] = lmi_portion(...
    alpha,...
    beta0, beta1, beta2, beta3, beta4, ...
    q1, q2, q3,...
    r1, r2)
%
% Usage: [beta0_or_inf] = lmi_portion(alpha, beta0, beta1, beta2, beta3, beta4, q1, q2, ...
%     q3, r1, r2)
% 
% Description: Internal loop that given 11 scalar values solves an LMI.
% 
% Input:
%     alpha = parameter associated with Lyapunov-Krasovskii stability.
%     beta0 = radius of the bounding ball
%     beta1 = 
%     beta2 = 
%     beta3 = 
%     beta4 = 
%     q1    = parameter associated with control of Q1 matrix variable
%     q2    = parameter associated with control of Q2 matrix variable
%     q3    = parameter associated with control of Q3 matrix variable
%     r1    = parameter associated with control of R1 matrix variable
%     r2    = parameter associated with control of R2 matrix variable
% 
% Output:
%     beta0_or_inft = beta0 if solvable LMI with parameters, otherwise
%                     infinite for outerloop optimization
    
    %% Linear Matrix Inequality Variable Setup
    setlmis([]);  % initialize the linear matrix inequality setup
    
    % Matrix constants
    beta12 = [beta1*eye(n), zeros(n); zeros(n), beta2*eye(n)];
    beta34 = [beta3*eye(n), zeros(n); zeros(n), beta4*eye(n)];
    
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
    [R1tilde, ~, sR1tilde] = lmivar(3, [sR1, zeros(n); zeros(n) sR1]);
    [R2tilde, ~, sR2tilde] = lmivar(3, [sR2, zeros(n); zeros(n) sR2]);
    Theta = lmivar(3, [sR2tilde, sX; sX', sR2tilde]);
    
    %% Linear Matrix Inequality Constraint Setup
    % Eq 6
    lmiterm([1 1 1 Q1]);
    lmiterm([-1 1 1 0], q1*eye(n));
    lmiterm([2 1 1 Q2]);
    lmiterm([-2 1 1 0], q2*eye(n));
    lmiterm([3 1 1 Q3]);
    lmiterm([-3 1 1 0], q3*eye(n));
    lmiterm([4 1 1 R1]);
    lmiterm([-4 1 1 0], r1*eye(n));
    lmiterm([5 1 1 R2]);
    lmiterm([-5 1 1 0], r2*eye(n));
    
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
    lmiterm([-7 2 2 0], beta34, 1);         % RHS diag{beta3 I_n, beta4 I_n}
    
    % Eq 8 [I/beta0^2, zeros; zeros]< [[(tau_M - tau_m)P_11^i, P_12];[*, P_22]]
    % i = 1
    lmiterm([-8 1 1 P111], (tau_M - tau_m));
    lmiterm([-8 1 2 P12], 1);
    lmiterm([-8 2 2 P22], 1);
    lmiterm([8 1 1 0], ones(n)/(beta0^2));
    
    % i = 2
    lmiterm([-9 1 1 P112], (tau_M - tau_m));
    lmiterm([-9 1 2 P12], 1);
    lmiterm([-9 2 2 P22], 1);
    lmiterm([9 1 1 0], ones(n)/(beta0^2));
    
    % Eq 9
    lmiterm([-10 1 1 Theta], 1);
    
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
            lmiterm([lmit_1 1 1 -P12],...
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
            lmiterm([lmit_i 1 1 R1tilde] -exp(-alpha*tau_m)*Gamma_1, Gamma_1');
    
            % e^{-alpha tau_M} Gamma_2 Theta Gamma_2^T
            lmiterm([lmit_i 1 1 Theta] -exp(-alpha*tau_M)*Gamma_2, Gamma_2');
    
            % Constant term (alpha/omega^2*e8 e8^T)
            lmiterm([lmit_i 1 1 0], alpha/(omega_scalar^2)*e8, e8');
        end
    end

    % Get the LMI problem and solve
    lmisys = getlmis;
    [tmin,~] = feasp(lmisys);
    if tmin < 0
        beta0_or_inf = beta0;
    else
        beta0_or_inf = inf;
    end
    end

%% Constrain function
    % Eq 11.
    function [value] = kellipse(alpha, beta1, beta2, beta3, beta4,...
                                q1, q2, q3, r1, r2)
        K1 = beta1 + (tau_m^2)*beta2 + ((tau_M - tau_m)^2)*(beta3 + beta4) +...
            (q1*(1-exp(-alpha*tau_m)) + (q2+q3)*(1-exp(-alpha*tau_M)))/alpha;
        K2 = (r1*tau_m*(tau_m*alpha + exp(-alpha*tau_m) - 1) + ...
            r2*(tau_M - tau_m)*((tau_M - tau_m)*alpha +...
            exp(-alpha*tau_M) - exp(-alpha*tau_m)))/(alpha^2);
        value = K1*mu0^2 + K2*mu_scalar^2;
    end

%% Outer Loop optimization problem over LMI constraint parameters
% Wrapper functions for objective function and constraints
objf = @(x) lmi_portion(x(1), ...
                        x(2), x(3), x(4), x(5), x(6),...
                        x(7), x(8), x(9),...
                        x(10), x(11));
keconstr = @(x) kellipse(x(1), ...
                         x(2), x(3), x(4), x(5), x(6),...
                         x(7), x(8), x(9),...
                         x(10), x(11));

% Set up optimization variable and problem
x = optimvar('x',1,11);
paramprob = optimproblem(...
    'Objective',fcn2optimexpr(objf, x),...
    'Constraints',fcn2optimexpr(keconstr,x) <= 1 ...
    );

% Add in optimization constraints
paramprob.Constraints.alpha = x(1) > 0;
paramprob.Constraints.beta0 = x(2) > 0;
paramprob.Constraints.beta1 = x(3) > 0;
paramprob.Constraints.beta2 = x(4) > 0;
paramprob.Constraints.beta3 = x(5) > 0;
paramprob.Constraints.beta4 = x(6) > 0;
paramprob.Constraints.q1 = x(7) > 0;
paramprob.Constraints.q2 = x(8) > 0;
paramprob.Constraints.q3 = x(9) > 0;
paramprob.Constraints.r1 = x(10) > 0;
paramprob.Constraints.r2 = x(11) > 0;

% Solve problem
x0.x = ones(1,11);
[sol,fval,exitflag,output] = solve(prob,x0);

beta0 = sol.x(2);

end

%% Auxiliary Functions
function ei = ei_generator(n, p, i)
    if i < 8
        ei = [zeros(n, (i-1)*n), eye(n), zeros(n, (7-i)*n + p)];
    else
        ei = [zeros(p, 7*n), eye(p)];
    end
end
