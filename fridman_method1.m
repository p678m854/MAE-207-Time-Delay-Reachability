function [P, F, delta, tmin, xfeas]= fridman_method1(A, Ad, B, h, ut, ...
                               lambda, gamma_0, gamma_1, gamma_2, gamma_3) 

% FRIDMAN_SHAKED_METHOD1
% 
% Usage: [P, F, delta, tmin, xfeas]= fridman_method1(A, Ad, B, h, ut, ...
%                              lambda, gamma_0, gamma_1, gamma_2, gamma_3)
% 
% Description: 
% 
%     d/dt(x) = Ax + Ad x(t - tau(t)) + B u
% 
% Intputs:  - Matrix A        nxn
%           - Matrix Ad       nxn
%           - Matrix B        nxp
%           - Scalar h        1x1
%           - Scalar ut       1x1
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
% TODO:
%     1. Write body of function -- FO
%     2. Get the function working  ¯\_(ツ)_/¯ -- FO
%
%%  Define Dimensions

n = size(A,1);
p = size(B,2);

%% Define Variables
% Set on the LMI-System
setlmis([]);

% LMI-variables
P = lmivar(1,[n 1]); 
W = lmivar(2,[n n]); % W = P*F
delta = lmivar(1,[1 0]);

% Other constants
I = eye(p,p);
I_n = eye(n,n);

%% Define LMI inputs
% Definition LMI 1: 
% LMI 1 Diagonal Entries
lmiterm([1 1 1 W],1,1,'s')
lmiterm([1 1 1 P],1,A,'s')
lmiterm([1 1 1 P],(lambda*ut + gamma_0 + gamma_3*h*ut + h*gamma_1 ...
                    + h*gamma_2),1)
lmiterm([1 2 2 0],-lambda*I)
lmiterm([1 3 3 P],-gamma_0,1)
lmiterm([1 4 4 P],-h*gamma_1,1)
lmiterm([1 5 5 P],-h*gamma_2,1)
lmiterm([1 6 6 0],-h*gamma_3*I)

% LMI 1 Upper Entries              
lmiterm([1 1 2 P],1,B)
lmiterm([1 1 3 P],1,Ad)
lmiterm([1 1 3 W],-1,1)
lmiterm([1 1 4 W],h,A)
lmiterm([1 1 5 W],h,Ad)
lmiterm([1 1 6 W],h,B)

% LMI 1 Lower Entries
lmiterm([1 2 1 P],B',1)
lmiterm([1 3 1 P],Ad',1)
lmiterm([1 3 1 W],-1,1)
lmiterm([1 4 1 W],h*A',1)
lmiterm([1 5 1 W],h*Ad',1)
lmiterm([1 6 1 W],h*B',1)

% Definition LMI 2: 
lmiterm([-2 1 1 delta],1,I_n)
lmiterm([-2 1 2 0],I_n)
lmiterm([-2 2 1 0],I_n)
lmiterm([-2 2 2 P],1,1)

%% Define the LMI System

lmis = getlmis;

%% Solving the LMI

[tmin,xfeas] = feasp(lmis);

%% Postprocessing
P = dec2mat(lmis,xfeas,P);
W = dec2mat(lmis,xfeas,W);
delta = dec2mat(lmis,xfeas,delta);
F = P\W;
end
