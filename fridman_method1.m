function [P, F, delta]= fridman_method1(A0, A1, B, ...
                                            h, wbar...
                                            ) 
% 
% FRIDMAN_METHOD1   Determines the bounding ellipsoid of a forward reachable
%                   set with time delay from an initial point x0
% 
% Usage: [P, F, delta] = fridman_method1(A0, A1, B, ...
%                                             h, wbar, ...
%                                             ) 
% 
% Description: Given a linear system with a single time delay described by
%              the differential equation:
% 
%     d/dt(x) = Ax + Ad x(t - tau(t)) + B w
% 
% Intputs:  
%     
%     A0      = nxn state matrix
%     A1      = nxn state matrix for the delayed states
%     B       = nxp input matrix
%     h       = scalar bounds on the delay such that the inequalities
%               0 <= tau(t) <= tau_M are observed.
%     wbar    = scalar bounds on the input vector w such that
%               w(t)^T w(t) <= wbar
%         
%           
% Outputs:  
%     P       = nxn symmetric matrix for the bounding ellipsoid
%     F       = nxn matrix
%     delta   = scalar
%     beta0   = scalar radius of the geometric avarage
% 
% Authors:
%     Florian Oberneder
% 
% Date:
%     October 28, 2021
% 
% References:
%     1. "On reachable sets for linear systems with dealy and bounded peak 
%        inputs", E. Fridman, U. Shaked, June 2003.
% 
% Requires:
%     Robost Control Toolbox
% 

%% Input checking
assert(length(size(A0))  == 2, "A0 must be a matrix.")
assert(length(size(A1)) == 2, "A1 must be a matrix.")
assert(length(size(B))  == 2, "B must be a matrix.")

n = size(A0, 1);
p = size(B, 2);


%% Inner loop
    function [tmin,xfeas]=lmi_prob(lambda,gamma_0,gamma_1,gamma_2,gamma_3)
    %% Linear Matrix Inequality Variable Setup
    setlmis([]); % initialize the linear matrix inequality setup

    % Matrix Variables
    P = lmivar(1,[n 1]); % Symmetric nxn matrix of full elements
    W = lmivar(2,[n n]); % Full nxn; W = P*F
    delta = lmivar(1,[1 0]); % Scalar

    % Matrix Constants
    I_p = eye(p,p);     % pxp Identity Matrix
    I_n = eye(n,n);     % nxn Identity Matrix

    %% Linear Matrix Inequality Constraint Setup
    % Definition LMI 1 (Eq 15a): 
    % LMI 1 Diagonal Entries
    lmiterm([1 1 1 W],1,1,'s')
    lmiterm([1 1 1 P],1,A0,'s')
    lmiterm([1 1 1 P],(lambda*wbar + gamma_0 + gamma_3*h*wbar + h*gamma_1 ...
                        + h*gamma_2),1)
    lmiterm([1 2 2 0],-lambda*I_p)
    lmiterm([1 3 3 P],-gamma_0,1)
    lmiterm([1 4 4 P],-h*gamma_1,1)
    lmiterm([1 5 5 P],-h*gamma_2,1)
    lmiterm([1 6 6 0],-h*gamma_3*I_p)

    % LMI 1 Upper Entries              
    lmiterm([1 1 2 P],1,B)
    lmiterm([1 1 3 P],1,A1)
    lmiterm([1 1 3 W],-1,1)
    lmiterm([1 1 4 W],h,A0)
    lmiterm([1 1 5 W],h,A1)
    lmiterm([1 1 6 W],h,B)

    % Definition LMI 2 (Eq 19): 
    lmiterm([-2 1 1 delta],1,I_n)
    lmiterm([-2 1 2 0],I_n)
    lmiterm([-2 2 2 P],1,1)

    %% Define the LMI System
    lmisys = getlmis;

    %% Solving the LMI
    target = 0;
    options = [0 0 0 0 0];
    [tmin,xfeas] = feasp(lmisys,options,target);
    
    %% Postprocessing
    P = dec2mat(lmisys,xfeas,P);
    W = dec2mat(lmisys,xfeas,W);
    delta = dec2mat(lmisys,xfeas,delta);
    F = P\W;

    end

%% Outer Loop optimization problem over LMI tuning parameters
% Set up objective function
objf = @(x) lmi_prob(x(1), x(2), x(3), x(4), x(5));
options = optimset('Display','final');

% Initial optimization guess of the paper
 x0 = [0.25; 0.2; 1; 1; 0.17];
% x0 =[lambda; gamma_0; gamma_1; gamma_2; gamma_3];

% Solve Problem
[x,fval,exitflag,output] = fminsearch(objf,x0,options)

end


