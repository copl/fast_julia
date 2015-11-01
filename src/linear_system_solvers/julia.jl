type linear_solver_JULIA <: abstract_linear_system_solver
	_SparseMatrix::SparseMatrixCSC{Float64,Int64}
	_factor # TO DO, give type

	# options
	sym::Int64

  function linear_solver_JULIA()
      return new();
  end
end

function ls_factor(solver::linear_solver_JULIA, n::Int64, m::Int64)
			if solver.sym == 0
				 solver._factor = lufact(solver._SparseMatrix);
			elseif solver.sym == 1
				solver._factor = cholfact(solver._SparseMatrix);
			else
				error("this.options.sym = " * string(solver.sym) * " not supported")
			end

			return 1; # inertia
end

function ls_solve!(solver::linear_solver_JULIA, my_rhs::Array{Float64,1}, my_sol::Array{Float64,1})
  #GLOBAL_timer::class_algorithm_timer
  #GLOBAL_timer.start("Solve/julia")
  my_sol[1:length(my_sol)] = solver._factor \ my_rhs; #::UmfpackLU{Float64,Int64}
  #GLOBAL_timer.stop("Solve/julia")
end
