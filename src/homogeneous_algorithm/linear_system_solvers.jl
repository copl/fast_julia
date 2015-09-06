using MUMPS, MATLAB

#### TO DO:
#
# - clean up normal equation stuff


######################################################################
# Each instance of abstract_linear_system_solver should contain 
#	initialize::Function
# 	ls_factor::Function
# 	ls_solve!::Function
# i.e. see class_mumps_linear_system_solver
#######################################################################


# standard function for solving either the whole system or normal equations

# solve whole system
function construct_K(qp::class_quadratic_program, my_settings::class_settings)
	return sparse([[qp._H qp.A']; [qp.A -settings.diagonal_modification*speye(qp.m,qp.m)]])
end

function update_K!(K::SparseMatrixCSC{Float64,Int64}, vars::class_variables, qp::class_quadratic_program)
	x = vars.x();
	s = vars.s();

	H = qp.H(vars);
	for i = 1:(qp.n)
		K[i,i] = H[i,i] + s[i]/x[i];
	end
end

# solve normal equations A*A'

function form_normal_equations(vars::class_variables, qp::class_quadratic_program)
	H = qp.H(vars);
	diH = diag(H);
	@assert( sum(diH .>  0) == length(nonzeros(H)) ) # check that nothing is on the diagonals of H

	invH = spdiagm( 1./(diH + vars.s() ./ vars.x() ) );
	
	GLOBAL_timer.start("Factor/Schur")
	normal_matrix = qp.A * invH * qp.A' + settings.diagonal_modification*speye(qp.m,qp.m);
	GLOBAL_timer.stop("Factor/Schur")

	return normal_matrix, invH;
end

function standard_ls_initialize(qp, my_settings, is_normal)
	if is_normal
		return speye(0);
	else
		return construct_K(qp, my_settings);
	end
end


function standard_ls_before_factor!(Sparse_Matrix, vars, qp, is_normal)
	#this.SparseMatrix, this.invH = before_factor(this.SparseMatrix, vars, qp, this.options.normal);
	if is_normal
		return form_normal_equations(vars, qp);
	else
		update_K!(Sparse_Matrix, vars,  qp);
		return Sparse_Matrix, speye(0)
	end
end

function standard_ls_before_solve(my_sol::Array{Float64}, my_rhs::Array{Float64}, invH::SparseMatrixCSC{Float64,Int64}, qp::class_quadratic_program, is_normal::Bool)
	if is_normal
		n = qp.n
		m = qp.m
		temp_rhs = A * invH * my_rhs[1:n] - my_rhs[(1 + n):(m + n)];
		temp_sol = my_sol[(n + 1):(n + m)];
		
		return temp_sol, temp_rhs
	else
		return my_sol, my_rhs
	end
end

function standard_ls_after_solve!(my_sol::Array{Float64}, temp_sol::Array{Float64}, my_rhs::Array{Float64}, invH::SparseMatrixCSC{Float64,Int64}, qp::class_quadratic_program, is_normal::Bool) 
	if is_normal
		my_sol[1:qp.n] = invH*(my_rhs[1:qp.n] - qp.A' * temp_sol);
		my_sol[(qp.n + 1):(qp.n + qp.m)] = temp_sol;
	end
	# don't need to do anything if equations are not normal
	return my_sol
end

#########################
# compute the inertia
#########################

function compute_inertia(diag_vec::Array{Float64}, num_vars::Int64, num_constraints::Int64)
	try
		@assert(length(diag_vec) == num_vars + num_constraints)
		pos_eigs = sum(diag_vec .> 0.0);
		neg_eigs = sum(diag_vec .< 0.0);
	
		if pos_eigs == num_vars && neg_eigs == num_constraints
			return true
		elseif pos_eigs > num_vars || neg_eigs + pos_eigs != num_vars + num_constraints
			println("Warning: numerical instability in LDL factorization")
			if pos_eigs > num_vars
				println("more positive eigenvalues than variables")
				return false
			else
				println("zero eigenvalues")
				return false
			end
		else
			return false
		end
	catch e
		println("ERROR in compute_inertia")
		throw(e)
	end
end

# norm_rhs = normal_rhs(rhs,qp)
# dy = factor \ norm_rhs
# dx = solve_normal(dy, b) 
# dx = H\ A'*

#########################################
# Linear solver options
#########################################

abstract abstract_linear_solver_options

type standard_linear_solver_options <: abstract_linear_solver_options
	sym::Int64
	normal::Bool

	function standard_linear_solver_options()
		this = new();

		this.sym = 0; # use LU factorization by default
		this.normal = false;

		return this;	
	end
end

##########################################################################
# THE DIFFERENT LINEAR SOLVERS
##########################################################################

####
#### add any new linear solver here ...
####

abstract abstract_linear_system_solver

type linear_solver_MUMPS <: abstract_linear_system_solver
	initialize::Function
	ls_factor::Function
	ls_solve!::Function

	SparseMatrix::SparseMatrixCSC{Float64,Int64}	
	factor::MUMPSfactorization{Float64}
	
	# options
	options::standard_linear_solver_options
	
	qp::class_quadratic_program
	invH::SparseMatrixCSC{Float64,Int64}		

	function linear_solver_MUMPS()
		this = new();
		this.options = standard_linear_solver_options();
		
		
		this.initialize = function(qp::class_quadratic_program, my_settings::class_settings)	
			this.qp = qp;	
			this.invH = speye(0);	
			this.SparseMatrix = standard_ls_initialize(qp, my_settings, this.options.normal);
		end

		this.ls_factor = function(vars::class_variables, qp::class_quadratic_program)
			#if isdefined(this,:factor) 	
			#	destroyMUMPS(this.factor);
			#end
			this.SparseMatrix, this.invH = standard_ls_before_factor!(this.SparseMatrix, vars, qp, this.options.normal);
			
			GLOBAL_timer::class_algorithm_timer
			GLOBAL_timer.start("Factor/MUMPs")
			this.factor = factorMUMPS(this.SparseMatrix, this.options.sym);
			GLOBAL_timer.stop("Factor/MUMPs")

			return 1; # inertia
		end
		
		this.ls_solve! = function(my_rhs::Array{Float64}, my_sol::Array{Float64}, qp::class_quadratic_program)
			is_normal = this.options.normal;

			temp_sol, temp_rhs = standard_ls_before_solve(my_sol, my_rhs, this.invH, this.qp, is_normal)

			GLOBAL_timer::class_algorithm_timer
			GLOBAL_timer.start("Solve/MUMPs")
			applyMUMPS!(this.factor, temp_rhs, temp_sol);
			GLOBAL_timer.stop("Solve/MUMPs")

			standard_ls_after_solve!(my_sol, temp_sol, my_rhs, this.invH, this.qp, is_normal) 			
		end
		
		return this;
	end
end

type linear_solver_julia <: abstract_linear_system_solver
	initialize::Function
	ls_factor::Function
	ls_solve!::Function

	SparseMatrix::SparseMatrixCSC{Float64,Int64}	
	factor


	# options
	options::standard_linear_solver_options

	qp::class_quadratic_program
	invH::SparseMatrixCSC{Float64,Int64}

	function linear_solver_julia()
		this = new();
		this.options = standard_linear_solver_options()
		this.options.sym = 0; # use LU factorization by default
		
		this.initialize = function(qp::class_quadratic_program, my_settings::class_settings)
			this.invH = speye(0);
			this.SparseMatrix = standard_ls_initialize(qp, my_settings, this.options.normal);
		end

		this.ls_factor = function(vars::class_variables,qp::class_quadratic_program)	
			this.SparseMatrix, this.invH = standard_ls_before_factor!(this.SparseMatrix, vars, qp, this.options.normal);
			
			if this.options.sym == 0
				this.factor = lufact(this.SparseMatrix)
			elseif this.options.sym == 1
				this.factor = cholfact(this.SparseMatrix)
			else
				error("this.options.sym = " * string(this.options.sym) * " not supported")
			end

			return 1; # inertia
		end
		
		this.ls_solve! = function(my_rhs::Array{Float64}, my_sol::Array{Float64}, qp::class_quadratic_program)
			is_normal = this.options.normal;

			temp_sol, temp_rhs = standard_ls_before_solve(my_sol, my_rhs, this.invH, qp, is_normal)

			GLOBAL_timer::class_algorithm_timer
			GLOBAL_timer.start("Solve/julia")
			temp_sol[1:length(temp_sol)] = this.factor \ temp_rhs; #::UmfpackLU{Float64,Int64}
			GLOBAL_timer.stop("Solve/julia")

			standard_ls_after_solve!(my_sol, temp_sol, my_rhs, this.invH, qp, is_normal) 
		end
		
		return this;
	end
end


type linear_solver_matlab <: abstract_linear_system_solver
	initialize::Function
	ls_factor::Function
	ls_solve!::Function

	SparseMatrix::SparseMatrixCSC{Float64,Int64}


	# options
	options::abstract_linear_solver_options

	function linear_solver_matlab()
		this = new();
		this.options = standard_linear_solver_options()
		this.options.sym = 0; # use LU factorization by default

		this.initialize = function(qp::class_quadratic_program, my_settings::class_settings)
			this.SparseMatrix = construct_K(qp, my_settings);
		end

		this.ls_factor = function(vars::class_variables, qp::class_quadratic_program)	
			update_K!(this.SparseMatrix, vars, qp);
			
			sparse_matrix = this.SparseMatrix
			@mput sparse_matrix
			
			@matlab begin
				mat_L, mat_D, mat_P, mat_S = ldl(sparse_matrix,1e-6);
			end

			@mget mat_D
			
			return compute_inertia(diag(mat_D), qp.n, qp.m)
		end
		
		this.ls_solve! = function(my_rhs::Array{Float64}, my_sol::Array{Float64}, qp::class_quadratic_program)
			@mput my_rhs

			@matlab begin
				mat_solution = mat_S * ( mat_P* (mat_L'\ (mat_D\ (mat_L\ (mat_P' * ( mat_S * my_rhs ) ) ) ) ) );
			end
		
			@mget mat_solution

			my_sol[1:length(my_sol)] =  mat_solution;
		end
		
		return this;
	end
end



