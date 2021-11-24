function [sol] = get_sol(usecase,method)
%
% GET_SOL   Is used to solve a defined problem for the choosen usecase
%   
% Usage: sol,time = get_sol(usecase,method)
% 
% Description
% 
% Inputs:
%
% Outputs:
%
% Authors: 
%   Florian Oberneder
%
% References:
%
% Requires:
%
% TODO:
%   1. Write body of the function -- FO
%   2. Specify the description -- FO

%% Usecase selection
switch usecase
   case 'pendulum_example'
      run('pendulum_example.m');
   case 'satellite_example'
      run('satellite_example.m');
   case 'simple_quadcopter_example'
      run('simple_quadcopter_example.m');
   case 'quadcopter_example'
      run('quadcopter_example.m');
   case 'fridman_example'
      run('fridman_example.m');
   otherwise
      disp('Error! No usecase selected!')
end

fprintf('Usecase: %s is selected \n',usecase);

%% Method selection
switch method
   case 'trinh_method1'
      disp('Method: TRINH_METHOD1 is selected');
      [beta0] = trinh_method1(A, Ad, B, ...
                            E,...
                            tau_m, tau_M,...
                            d_m, d_M,...
                            mu_scalar,...
                            omega_scalar...
                            );
      sol = beta0;                                  
   case 'fridman_method1'
      disp('Method: FRIDMAN_METHOD1 is selected');
      [P, F, delta]= fridman_method1(A, Ad, B, ...
                            tau_M, sqrt(omega_scalar)...
                            );
      sol = P;
   otherwise
      disp('Error! No method selected!')
end

end

