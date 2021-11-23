function [P, F, delta, tmin, xfeas]= fridman_method1(A, Ad, B, ...
                                            h, utilde, ...
                                            lambda, gamma_0, ...
                                            gamma_1, gamma_2, ...
                                            gamma_3 ...
                                            ) 

% FRIDMAN_METHOD1   Determines the bounding ellipsoid of a forward reachable
%                   set with time delay from an initial point x0
% 
% Usage: [P, F, delta, tmin, xfeas]= fridman_method1(A, Ad, B, h, utilde, ...
%                              lambda, gamma_0, gamma_1, gamma_2, gamma_3)
% 
% Description: Given a linear system with a single time delay described by
%              the differential equation:
%
%     d/dt(x) = Ax + Ad x(t - tau(t)) + B u
% 
% Intputs:  - Matrix A        nxn
%           - Matrix Ad       nxn
%           - Matrix B        nxp
%           - Scalar h        1x1
%           - Scalar utilde   1x1
%           - Scalar lambda   1x1
%           - Scalar gamma_0  1x1
%           - Scalar gamma_1  1x1
%           - Scalar gamma_2  1x1
%           - Scalar gamma_3  1x1
%           
% Outputs:  - Matrix P        nxn
%           - Matrix F        nxn
%           - Scalar delta    nxn
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
% TODO:
%     1. Write body of function -- FO
%     2. Get the function working -- FO
%     3. Update names and description -- FO  
%

%% Input checking
assert(length(size(A))  == 2, "A must be a matrix.")
assert(length(size(Ad)) == 2, "Ad must be a matrix.")
assert(length(size(B))  == 2, "B must be a matrix.")

n = size(A, 1);
p = size(B, 2);

%% Linear Matrix Inequality Variable Setup
setlmis([]); % initialize the linear matrix inequality setup

% Matrix Variables
P = lmivar(1,[n 1]); % Symmetric nxn matrix of full elements
W = lmivar(2,[n n]); % Full nxn; W = P*F
delta = lmivar(1,[1 0]); % Scalar

% Matrix Constants
I_p = eye(p,p);       % pxp Identity Matrix
I_n = eye(n,n);     % nxn Identity Matrix

%% Linear Matrix Inequality Constraint Setup
% Definition LMI 1 (Eq 15a): 
% LMI 1 Diagonal Entries
lmiterm([1 1 1 W],1,1,'s')
lmiterm([1 1 1 P],1,A,'s')
lmiterm([1 1 1 P],(lambda*utilde + gamma_0 + gamma_3*h*utilde + h*gamma_1 ...
                    + h*gamma_2),1)
lmiterm([1 2 2 0],-lambda*I_p)
lmiterm([1 3 3 P],-gamma_0,1)
lmiterm([1 4 4 P],-h*gamma_1,1)
lmiterm([1 5 5 P],-h*gamma_2,1)
lmiterm([1 6 6 0],-h*gamma_3*I_p)

% LMI 1 Upper Entries              
lmiterm([1 1 2 P],1,B)
lmiterm([1 1 3 P],1,Ad)
lmiterm([1 1 3 W],-1,1)
lmiterm([1 1 4 W],h,A)
lmiterm([1 1 5 W],h,Ad)
lmiterm([1 1 6 W],h,B)


% % Definition LMI 2 (Eq 19): 
% lmiterm([-2 1 1 P],1,1)
% lmiterm([2 1 2 0],I_n)
% % lmiterm([-2 2 1 0],I_n)
% % lmiterm([-2 2 2 P],1,1)


% Definition LMI 2: 
lmiterm([-2 1 1 delta],1,I_n)
lmiterm([-2 1 2 0],I_n)
lmiterm([-2 2 2 P],1,1)

%% Define the LMI System
lmisys = getlmis;

%% Solving the LMI
target = 0;
options = [0 0 0 0 0];
[tmin,xfeas] = feasp(lmisys,options,target);
%[tmin,xfeas] = gevp(lmisys,1);

%% Postprocessing
P = dec2mat(lmisys,xfeas,P);
W = dec2mat(lmisys,xfeas,W);
delta = dec2mat(lmisys,xfeas,delta);
F = P\W;
end
