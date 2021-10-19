function [beta0] = trinh_method1(A, Ad, B, ...
                            x0,...
                            tau_m, tau_M,...
                            d_m, d_M,...
                            mu_scalar,...
                            omega_scalar,...
                            )

TRINH_METHOD1 Determines a forward reachable set of a linear system with
              time delays from an initial point x0.

Usage: beta0 = trinh_method1(A, Ad, B, x0, tau_m, tau_M, d_m, d_M, ...
                             mu_scalar, omega_scalar)

Description: 

    d/dt(x) = Ax + Ad x(t - tau(t)) + B u

Intputs:
Outputs:

Authors:
    Patrick McNamee

Date:
    October 19, 2021

References:
    1. "On backwards and forwards reachable sets bounding for perturbed 
       time-delay systems", H. Trihn, Phan T. Nam, Pubudu N. Pathirana,
       and H. P. Le, June 2015.

TODO:
    1. Write body of function -- PM

end