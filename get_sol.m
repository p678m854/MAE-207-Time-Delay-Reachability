function [sol,time] = get_sol(usecase,method,opt)
%
% GET_SOL   Is used to solve a defined problem for the choosen usecase
%   
% Usage: sol,time = get_sol(usecase,method,opt)
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
      load('pendulum_example.mat');
   case 'satellite_example'
      load('satellite_example.mat');
   case 'quadcopter_example'
      load('quadcopter_example.mat');
   otherwise
      disp('Please select a usecase')
end

fprintf('usecase %s is selected',usecase);

%% Method selection
switch method
   case 'trinh_method1'
      disp('method TRINH_METHOD1 is selected');
      statements
   case 'fridman_method1'
      disp('method FRIDMAN_METHOD1 is selected');
      statements
   otherwise
      disp('Please select a method')
end

fprintf('method %s is selected',usecase);
end

