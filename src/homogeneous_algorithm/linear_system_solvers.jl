using MUMPS


######################################################################
# Each instance of abstract_linear_system_solver should contain 
# 	ls_factor::Function
# 	ls_solve!::Function
# i.e. see class_mumps_linear_system_solver
#######################################################################

abstract abstract_linear_system_solver

##########################################################################
# THE DIFFERENT LINEAR SOLVERS
##########################################################################

####
#### add any new linear solver here ...
####

type linear_solver_MUMPS <: abstract_linear_system_solver

	ls_factor::Function
	ls_solve!::Function

	# options
	sym::Int64 # 0, 1, 2

	function linear_solver_MUMPS()
		this = new();
		this.sym = 0; # use LU factorization by default
		
		this.ls_factor = function(SparseMatrix::SparseMatrixCSC{Float64,Int64}, my_settings::class_settings)
			return factorMUMPS(SparseMatrix, this.sym);
		end
		
		this.ls_solve! = function(my_factor::MUMPSfactorization{Float64}, my_rhs::Array{Float64}, my_sol::Array{Float64}, my_settings::class_settings)
			applyMUMPS!(my_factor, my_rhs, my_sol);
		end
		
		return this;
	end
end

type linear_solver_julia <: abstract_linear_system_solver

	ls_factor::Function
	ls_solve!::Function

	# options
	sym::Int64 # 0, 1, 2

	function linear_solver_julia()
		this = new();
		this.sym = 0; # use LU factorization by default
		
		this.ls_factor = function(SparseMatrix::SparseMatrixCSC{Float64,Int64}, my_settings::class_settings)
			return lufact(SparseMatrix);
		end
		
		this.ls_solve! = function(my_factor, my_rhs::Array{Float64}, my_sol::Array{Float64}, my_settings::class_settings)
			my_sol[1:length(my_sol)] = my_factor \ my_rhs; #::UmfpackLU{Float64,Int64}
		end
		
		return this;
	end
end


##########################################################################
### CODE THAT CALLS DIFFERENT LINEAR SOLVERS
##########################################################################

function ls_factor(SparseMatrix::SparseMatrixCSC{Float64,Int64}, my_settings::class_settings)
	try
		my_settings.linear_system_solver.ls_factor(SparseMatrix, my_settings);
	catch e
		println("ERROR in ls_factor")
		throw(e)
	end
end

function ls_solve!(my_factor, my_rhs::Array{Float64}, my_sol::Array{Float64}, my_settings::class_settings)
	try
		my_settings.linear_system_solver.ls_solve!(my_factor, my_rhs, my_sol, my_settings);
	catch e
		println("ERROR in ls_factor")
		throw(e)
	end
end


