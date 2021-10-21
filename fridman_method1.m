function [beta0] = fridman_method1(A, Ad, B, ...
                            x0,...
                            tau_m, tau_M,...
                            d_m, d_M,...
                            mu_scalar,...
                            omega_scalar,...
                            )

% FRIDMAN_METHOD1 TBD
% 
% Usage: beta0 = fridman_method1(A, Ad, B, x0, tau_m, tau_M, d_m, d_M, ...
%                              mu_scalar, omega_scalar)
% 
% Description: 
% 
%     d/dt(x) = Ax + Ad x(t - tau(t)) + B u
% 
% Intputs:
% Outputs:
% 
% Authors:
%     Florian Oberneder
% 
% Date:
%     October 21, 2021
% 
% References:
%     1. "On reachable sets for linear systems with dealy and bounded peak 
%        inputs", E. Fridman, U. Shaked, June 2003.
% 
% TODO:
%     1. Write body of function -- FO
% 
% end
