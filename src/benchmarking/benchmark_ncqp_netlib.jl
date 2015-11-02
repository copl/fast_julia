#
# test the linear program solver on a simple problem
#
using Ipopt, KNITRO
#simple_tests();

# 54
# 59
#

function test_problem(name::String)
	file_name = name * ".mat";
	A, b, c = get_netlib_problem("/home/oliver/Documents/Programming_Projects/fast_julia/src/benchmarking/Problems", file_name);
	#Q = -speye(length(c));
	Q = tridiagonal(length(c),0.0,1.0);
	#Q = spzeros(length(c),length(c));
	println("Solving ", file_name, " with the homogeneous algorithm")
	println(size(A,2), " variables and ", size(A,1), " constraints")
	println("Non-zeros: ", length(nonzeros(A)))
	settings = class_settings();

	println("=================== Linear system solver is matlab ldl ==================")
	settings.linear_system_solver = linear_solver_MATLAB();

	settings.newton_solver = class_homogeneous_newton();
	trivial_test(A, b, c, Q, 1, file_name, settings, true);

	#settings.newton_solver = class_newton_ip();
	#trivial_test(A, b, c, Q, 1, file_name, settings, true);

	println("=========================================================================")
	println("Calls IPOPT")
	solve_with_JuMP(A, b, c, Q, IpoptSolver(max_iter=10,print_level=3));
end
