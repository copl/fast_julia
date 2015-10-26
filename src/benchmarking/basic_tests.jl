#
# test the linear program solver on a simple problem
#

include("../src/homogeneous_algorithm/homogeneous_algorithm.jl")
include_print("benchmarking_tools.jl")

function validate_functions()
	#################
	# simple tests to check that basic elements are working
	#################

	# quadratic program
	qp = class_quadratic_program();

	qp.n = 5;
	qp.m = 2;
	qp._H = spzeros(5,5);
	qp.A = spzeros(2, 5);
	qp.A[1, 1:2] = 1;
	qp.A[2, 3:5] = 1;
	qp.b = [1; 1];
	qp.q = rand(5);
	
	# settings
	settings = class_settings();

	# variables
	vars = class_variables(5, 2);
	vars.check_positive();

	# qp and vars
	validate_dimensions(qp, vars);

	# residuals
	res = class_residuals();
	res.update(qp, vars);

	# newton_solver
	newton_solver = class_newton_solver(qp, settings);
	newton_solver.update_system(vars);
	newton_solver.compute_direction(vars,res,0.5);
	
	# maximum step
	println("0.5 = ", maximum_step(1.0, 1.0, -2.0))

	# solver
	homogeneous_algorithm(qp, vars, settings)
end


function trivial_lp(settings::class_settings)
	# trivial_test(A, b, c,  correct_status, problem_name, verbose)
	trivial_test(sparse([1.0 1.0]), [1.0], [1.0, 1.0], 1, "LP-OPTIMAL-A", settings);
	trivial_test(sparse([1.0 -1.0]), [1.0], [1.0, 5.0], 1, "LP-OPTIMAL-B", settings);
	trivial_test(sparse([1.0 5.0]), [10.0], [1.0, -1.0], 1, "LP-OPTIMAL-C", settings);
	trivial_test(sparse([1.0 1.0]), [0.0], [1.0, -1.0], 1, "LP-OPTIMAL-D", settings);
	
	trivial_test(sparse([1.0 1.0]), [-1.0], [1.0, 1.0], 2, "LP-INFEASIBLE-A", settings);
	trivial_test(sparse([1.0 5.0]), [-2.0], [-1.0, 1.0], 2, "LP-INFEASIBLE-B", settings);
	trivial_test(sparse([1.0 1.0]), [-5.0], [0.0, 0.0], 2, "LP-INFEASIBLE-C", settings);
	trivial_test(sparse([1.0 1.0]), [-0.001], [3.0, 1.0], 2, "LP-INFEASIBLE-D", settings);
	
	trivial_test(spzeros(0,1), EMPTY_ARRAY, [-1.0], 3, "LP-UNBOUNDED-A", settings);
	trivial_test(sparse([1.0 -5.0]), [0.0], [0.0, -1.0], 3, "LP-UNBOUNDED-B", settings);
	trivial_test(sparse([1.0 -3.0]), [1.0], [0.5, -2.0], 3, "LP-UNBOUNDED-C", settings);
	trivial_test(sparse([1.0 -1.0]), [1.0], [0.0, -0.01], 3, "LP-UNBOUNDED-D", settings);
end


function trivial_qp(settings::class_settings)
	Q1 = sparse([ [1.0 -1.0]; [-1.0 1.0] ]);
	Q2 = sparse(ones(2,2));
	Q3 = sparse([ [1.0 0.0]; [0.0 0.0] ]);
		
	# trivial_test(A, b, c, Q,  correct_status, problem_name, verbose)
	trivial_test(sparse([1.0 1.0]), [1.0], [1.0, 1.0], speye(2), 1, "QP-OPTIMAL-A", settings);
	trivial_test(sparse([1.0 -1.0]), [0.0], [-1.0, -1.0], Q3, 1, "QP-OPTIMAL-B", settings);
	trivial_test(sparse([1.0 -3.0]), [1.0], [0.5, -2.0], Q2, 1, "QP-OPTIMAL-C", settings);
	trivial_test(sparse([1.0 1.0]), [0.0], [1.0, -1.0], Q1, 1, "QP-OPTIMAL-D", settings);

	trivial_test(sparse([1.0 1.0]), [-1.0], [1.0, 1.0], speye(2), 2, "QP-INFEASIBLE-A", settings);
	trivial_test(sparse([1.0 5.0]), [-1.0], [-1.0, -1.0], Q3, 2, "QP-INFEASIBLE-B", settings);
	trivial_test(sparse([1.0 3.0]), [-1.0], [0.5, -2.0], Q2, 2, "QP-INFEASIBLE-C", settings);
	trivial_test(sparse([1.0 1.0]), [-0.0001], [1.0, -1.0], Q1, 2, "QP-INFEASIBLE-D", settings);

	
	trivial_test(sparse([0.0 0.0]), [0.0], [0.0, -1.0], Q3, 3, "QP-UNBOUNDED-A", settings);
	trivial_test(sparse([1.0 0.0]), [5.0], [0.0, -1.0], Q3, 3, "QP-UNBOUNDED-B", settings);
	trivial_test(sparse([1.0 -1.0]), [1.0], [0.5, -2.0], Q1, 3, "QP-UNBOUNDED-C", settings);
	trivial_test(sparse([1.0 -1.0]), [0.0], [0.0, -0.01], Q1, 3, "QP-UNBOUNDED-D", settings);
end

function trivial_ncqp(settings::class_settings)
	Q1 = sparse([ [-1.0 1.0]; [1.0 -1.0] ]);
	Q2 = sparse(ones(2,2));
	Q3 = sparse([ [0.0 1.0]; [1.0 0.0] ]);
	
	
	# trivial_test(A, b, c, Q,  correct_status, problem_name, verbose)
	if false
		trivial_test(sparse([1.0 1.0]), [1.0], [1.0, 1.0], speye(2), 1, "NCQP-OPTIMAL-A", settings);
		trivial_test(sparse([1.0 1.0]), [1.0], [-1.0, -1.0], Q3, 1, "NCQP-OPTIMAL-B", settings);
		trivial_test(sparse([1.0 -3.0]), [1.0], [0.5, -2.0], Q2, 1, "NCQP-OPTIMAL-C", settings);
		trivial_test(sparse([1.0 1.0]), [0.0], [1.0, -1.0], Q1, 1, "NCQP-OPTIMAL-D", settings);

		trivial_test(sparse([1.0 1.0]), [-1.0], [1.0, 1.0], speye(2), 2, "NCQP-INFEASIBLE-A", settings);
		trivial_test(sparse([1.0 5.0]), [-1.0], [-1.0, -1.0], Q3, 2, "NCQP-INFEASIBLE-B", settings);
		trivial_test(sparse([1.0 3.0]), [-1.0], [0.5, -2.0], Q2, 2, "NCQP-INFEASIBLE-C", settings);
		trivial_test(sparse([1.0 1.0]), [-0.0001], [1.0, -1.0], Q1, 2, "NCQP-INFEASIBLE-D", settings);


		trivial_test(sparse([0.0 0.0]), [0.0], [0.0, -1.0], Q3, 3, "NCQP-UNBOUNDED-A", settings);
		trivial_test(sparse([1.0 0.0]), [5.0], [0.0, -1.0], Q3, 3, "NCQP-UNBOUNDED-B", settings);
		trivial_test(sparse([1.0 -1.0]), [1.0], [0.5, -2.0], Q1, 3, "NCQP-UNBOUNDED-C", settings);
	end
	trivial_test(sparse([1.0 -1.0]), [1.0], [0.0, -0.01], Q1, 3, "NCQP-UNBOUNDED-D", settings, true);

end


settings = class_settings();

settings.newton_solver = class_newton_hsd();
trivial_ncqp(settings);

settings.newton_solver = class_newton_ip();
trivial_ncqp(settings);






#trivial_lp(settings);
#trivial_qp(settings);








