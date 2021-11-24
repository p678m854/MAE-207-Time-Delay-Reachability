%% Execution script
%
%       This script is select the usecase, select the method and 
%   to run the execution of the selected method

%% Usecase selection
% Select the Usecase for execution
% Options:
%     - 2D usecase: pendulum_example
%     - 4D usecase: satellite_example
%     - 9D usecase: simple_quadcopter_example
%     - 12D usecase: quadcopter_example
% Additional Options for code verification:
%     - 2D usecase: fridman_example
% Uncomment to select the usecase:

% usecase = 'pendulum_example';
% usecase = 'satellite_example';
% usecase = 'simple_quadcopter_example';
% usecase = 'quadcopter_example';

% usecase = 'fridman_example';

%% Method selection
% Select a Mehtod for execution
% Options.
%     - trinh_method1
%     - fridman_method1
% Uncomment to select the method:

% method = 'trinh_method1';
method = 'fridman_method1';

%% Execution
%
tic
    [sol] = get_sol(usecase,method);
time = toc;

fprintf('The Method %s has an execution time of %f seconds for Usecase %s \n',method,time,usecase);
