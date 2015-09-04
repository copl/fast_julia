#
# test the linear program solver on a simple problem
#

include("homogeneous_algorithm.jl")

#simple_tests();

file_name = "Problems/QAP12.mat";

A, b, c = get_netlib_problem(file_name);
println("Solving ", file_name, " with the homogeneous algorithm")
println(size(A,2), " variables and ", size(A,1), " constraints")
println("Non-zeros: ", length(nonzeros(A)))
settings = class_settings();

println("=================== Linear system solver is MUMPS chol (normal system) ==================")
settings.linear_system_solver = linear_solver_MUMPS();
settings.linear_system_solver.options.normal = true;
settings.linear_system_solver.options.sym = 1;
lp_test_homogeneous_algorithm(A, b, c, settings);

using Mosek
println("=========================================================================")
println("Calls Mosek")
solve_with_JuMP(A, b, c, MosekSolver());

using Gurobi
println("=========================================================================")
println("Calls Gurobi")
solve_with_JuMP(A, b, c, GurobiSolver(Method=2));

#using Ipopt
#println("=========================================================================")
#println("Calls IPOPT")
#solve_with_JuMP(A, b, c, IpoptSolver());








