clear, clc,

n_methods = 2;
n_usecases = 3;

method_solutions = {};

for i = 2:n_methods
    if i == 1
        method = 'fridman_method1';
    elseif i == 2
        method = 'trinh_method1';
    elseif i == 3
        method = 'nam_method';
    else
        error("Unspecified error")
    end

    profile_sol_struct = {};

    for j = 1:n_usecases
        if j == 1
            usecase = 'pendulum_example';
        elseif j == 2
            usecase = 'satellite_example';
        elseif j == 3
            usecase = 'simple_quadcopter_example';
        else
            error("Unspecified use case.")
        end

        sol_struct = {};

        tic,
            sol = get_sol(usecase, method);
        run_time = toc;

        sol_struct.n = 2*j;  % works for j=1:3
        sol_struct.sol = sol;
        sol_struct.runtime = run_time;
        
        profile_sol_struct = vertcat(profile_sol_struct, sol_struct);
    end

    method_solutions = vertcat(method_solutions, {{method}, profile_sol_struct});
end

save('profileresults.mat', "method_solutions");
