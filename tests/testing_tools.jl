EMPTY_ARRAY = spzeros(0,1)*[0.0];

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
function test_homogeneous_algorithm(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, c::Array{Float64,1}, settings::class_settings)
	return test_homogeneous_algorithm(A, b, c, spzeros(length(c),length(c)), settings);
end

function test_homogeneous_algorithm(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, c::Array{Float64,1}, Q::SparseMatrixCSC{Float64,Int64}, settings::class_settings)
	try	
		qp = class_quadratic_program(A, b, c, Q);

		vars = class_variables(qp.n,qp.m);
		if settings.verbose
			println(qp.n , " variables, ", qp.m, " constraints")
		end
		
		return homogeneous_algorithm(qp, vars, settings)
	catch e
		println("ERROR in test_homogeneous_algorithm")
		throw(e)
	end
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

function trivial_test(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, c::Array{Float64,1}, correct_status::Int64, problem_name = "", settings::class_settings = class_settings(), verbose = false)
	return trivial_test(A, b, c, spzeros(length(c),length(c)), correct_status, problem_name, settings)
end

function trivial_test(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, c::Array{Float64,1}, Q::SparseMatrixCSC{Float64,Int64}, correct_status::Int64, problem_name = "", settings::class_settings = class_settings(), verbose = false)
	try
		
		
		#settings.linear_system_solver = linear_solver_MATLAB();
		#settings.linear_system_solver.options.sym = 0;
		settings.verbose = verbose;
		
		status, vars, iterations = test_homogeneous_algorithm(A, b, c, Q, settings);
	
		if status == correct_status
			println("Problem ", problem_name, " solved successfully in ", iterations, " iterations")
		else
			println("ERROR - Problem ", problem_name, " not solved sucessfully.")
			println("Expected status ", correct_status, " got status ", status)
			#println("x_scaled = ", vars.x()/vars.tau())		
		end

		return status, vars, iterations
	catch e
		println("ERROR in trivial_test")
		println("===================================")
		throw(e)
	end
end

