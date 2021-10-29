%% Main Document
% Input document for the functions:

rho = 0; % rho < |0.2|
A = [-2 0; 0, -0.9+rho];
Ad = [-1, 0; -1, -1+0.5*rho];
B = [-0.5; 1];
h = 0.7;
ut = 1;
lambda = 0.25;
gamma_0 = 0.2;
gamma_1 = 1;
gamma_2 = 1;
gamma_3 = 0.17;

[P, F, delta, tmin, xfeas]=fridman_method1(A, Ad, B, h, ut, ...
                               lambda, gamma_0, gamma_1, gamma_2, gamma_3) 

