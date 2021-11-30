function q = nam_method(A,Ak,B,omega_bar)
% 
% nam_method Determines a forward reachable set of a linear system with
%               time delays from an initial point x0. Returns a bounding
%               box.
% 
% Usage: bounding_box = nam_method(A,Ak,B,omega_bar)
% 
% Description: Given a linear system with a single time delay described by
%              the differential equation:
% 
%     d/dt(x) = Ax + (\Sum_k Ad_k x(t - tau_k(t))) + B omega(t)
% 
%             the nam_method finds a bounding box for the forward
%             reachability set. It does assume that the system is a priori
%             known to be stable under time delays.
% 
% Intputs:
% 
%     A            = nxn state matrix
%     Ak           = nxn state matrices for the k delayed states. Third axis
%                    represents the delayed state index.
%     B            = nx1 input matrix
%     omega_scalar = Bounds on the disturbance scalar omega such that
%                    |omega(t)| <= omega_bar for all t >= 0
% 
% Outputs:
%     q = Bounding box on the possible state vector where |x_i(t)| <= q_i 
%         for all t >= 0.
% 
% Authors:
%     Patrick McNamee
% 
% Date:
%     November 30, 2021
% 
% References:
%     1. "Reachable set bounding for nonlinear perturbed time-delay systems:
%        The smallest bound", P.T. Nam, Pubudu N. Pathirana, and H. Trihn, 
%        November 2014.


sumA = A;
for k=1:size(Ak,3)
    sumA = sumA + abs(Ak(:,:,k));
end

q = -sumA\abs(B)*omega_bar;

end