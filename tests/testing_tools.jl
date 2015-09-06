function simple_tests()
	#################
	# simple tests
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
	qp.c = rand(5);
	
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


function get_netlib_problem(file_name::String)
	file = matopen(file_name)
	A= sparse(read(file, "A"))
	b = read(file,"b")
	b = b[:]
	c = read(file,"c")
	c = c[:];
	
	return A, b, c
end


using MAT
function lp_test_homogeneous_algorithm(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, c::Array{Float64,1}, settings = class_settings())
	qp = class_quadratic_program(A, b, c);

	vars = class_variables(qp.n,qp.m);
	
	println(qp.n , " variables, ", qp.m, " constraints")
	homogeneous_algorithm(qp, vars, settings)
end

function qp_test_homogeneous_algorithm(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, c::Array{Float64,1}, Q::SparseMatrixCSC{Float64,Int64}, settings = class_settings())
	qp = class_quadratic_program(A, b, c, Q);

	vars = class_variables(qp.n,qp.m);
	
	println(qp.n , " variables, ", qp.m, " constraints")
	homogeneous_algorithm(qp, vars, settings)
end

function tridiagonal(k,di,offdi)
	Q = spzeros(k,k)
	for j = 1:(k-1)
		Q[j,j+1] = offdi	
		Q[j,j] = di
		Q[j+1,j] = offdi
	end

	return Q;
end


using JuMP, Ipopt

function solve_with_JuMP(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, c::Array{Float64,1}, solver=IpoptSolver(max_iter=300)) 
	solve_with_JuMP(A, b, c, spzeros(length(c),length(c)), solver) 
end

function solve_with_JuMP(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, c::Array{Float64,1}, Q::SparseMatrixCSC{Float64,Int64}, solver=IpoptSolver(max_iter=300)) 
	model = Model(solver=solver)
	
	n, k = size(A)
	
	@defVar(model, x[1:k] >= 0)
	@setObjective( model, Min, sum{c[i]*x[i], i=1:k} + sum{0.5*x[j]*Q[i,j]*x[i], i=1:k, j=1:k})

	
	#@addConstraint(model, A*x .== b)
	@defConstrRef constr[1:n]
	for i = 1:n
		constr[i] = @addConstraint(model, sum{A[i,j]*x[j], j=1:k} == b[i])
	end
		
	status = solve(model)
	
	if status == :Optimal
		return 1, getValue(x)[:], getDual(constr)[:]
	elseif status == :Infeasible
		return 2
	elseif status == :Unbounded
		return 3
	else
		return 0
	end
end
