include("homogeneous_algorithm.jl")

#simple_tests();

file_name = "Problems/QAP8.mat";

A, b, c = get_netlib_problem(file_name);
println("Solving ", file_name, " with the homogeneous algorithm")
settings = class_settings();

println("=================== Linear system solver is julia LU factor =============")
settings.linear_system_solver = :solver_julia
lp_test_homogeneous_algorithm(A, b, c, settings);

println("=================== Linear system solver is MUMPS LU ====================")
settings.linear_system_solver = :solver_MUMPS_LU
lp_test_homogeneous_algorithm(A, b, c, settings);

println("=================== Linear system solver is MUMPS LDL ===================")
settings.linear_system_solver = :solver_MUMPS_LDL
lp_test_homogeneous_algorithm(A, b, c, settings);

println("=========================================================================")

solve_with_JuMP(A, b, c)

