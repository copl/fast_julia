#
# test the linear program solver on a simple problem
#

include("../src/homogeneous_algorithm/homogeneous_algorithm.jl")
include_print("testing_tools.jl")

#simple_tests();

# 54
# 59
#


file_name = "Problems/STANDGUB.mat";
#for file_name in ["Problems/AFIRO.mat", "Problems/BLEND.mat", "Problems/CAPRI.mat", "Problems/BRANDY.mat"]
A, b, c = get_netlib_problem(file_name);
#Q = -speye(length(c));
Q = tridiagonal(length(c),0.0,1.0);
#Q = spzeros(length(c),length(c));

println("Solving ", file_name, " with the homogeneous algorithm")
println(size(A,2), " variables and ", size(A,1), " constraints")
println("Non-zeros: ", length(nonzeros(A)))
settings = class_settings();

println("=================== Linear system solver is matlab ldl ==================")
settings.linear_system_solver = linear_solver_MATLAB();
settings.linear_system_solver.options.sym = 2;

settings.newton_solver = class_newton_hsd();
trivial_test(A, b, c, Q, 1, file_name, settings, true);

settings.newton_solver = class_newton_ip();
trivial_test(A, b, c, Q, 1, file_name, settings, true);

using Ipopt
println("=========================================================================")
println("Calls IPOPT")
solve_with_JuMP(A, b, c, Q, IpoptSolver());
#end






