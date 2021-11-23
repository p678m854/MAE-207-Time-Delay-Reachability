%% Execution script
%
%       This script is select the usecase, select the method and 
%   to run the execution of the selected method

%% Usecase selection
% Select the Usecase for execution
% Options:
%     - 2D usecase: pendulum_example
%     - 4D usecase: satellite_example
%     - 10D usecase: quadcopter_example
% Uncomment to select the usecase:

usecase = 'pendulum_example';
% usecase = 'satellite_example';
% usecase = 'quadcopter_example';

%% Method selection
% Select a Mehtod for execution
% Options.
%     - trinh_method1
%     - fridman_method1
% Uncomment to select the method:

 method = 'trinh_method1';
% method = 'fridman_method1';

%% Execution
%
tic
    [sol,time] = get_sol(usecase,method,opt);
toc

fprintf('The method %s has an execution time of %f seconds',method,timeit);
