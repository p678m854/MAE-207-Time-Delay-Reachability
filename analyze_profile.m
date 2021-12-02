clear, clc

load("profileresults.mat")

for i = 1:size(method_solutions, 1)
    fprintf("Results for '%s'\n", string(method_solutions(i,1)));
    profile_sols = method_solutions(i,2);
    profile_sols = profile_sols{1,1};

    fprintf("\tn\ttime [s]\tr [~]\n")
    for j = 1:size(profile_sols, 1)
        sols = profile_sols(j,1);
        sols = sols{1,1};

        n = sols.n;
        runtime = sols.runtime;
        r = sols.sol;
        if (size(r,1) > 1 || size(r,2) > 1) && (size(r,1) ~= size(r,2))
            r = sum(abs(r), 2);
            np = sum(r > 0);
            r = prod(r(r > 0)).^(1/np);
        elseif size(r,1) == size(r,2) && size(r,1) > 1  % assuming symmetric
            r = ((pi^(n/2))/factorial(n/2)*sqrt(det(r)))^(1/n);
        end
        fprintf("\t%i\t%f\t%f\n", n, runtime, r)
    end
end

