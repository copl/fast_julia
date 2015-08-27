using MUMPS

function ls_factor(SparseMatrix::SparseMatrixCSC{Float64,Int64}, my_settings::class_settings)
	try
		if 	my_settings.linear_system_solver == :solver_MUMPS_LU
			return factorMUMPS(SparseMatrix, 0);
		elseif  my_settings.linear_system_solver == :solver_MUMPS_LDL
			return factorMUMPS(SparseMatrix, 2);
		elseif my_settings.linear_system_solver == :solver_julia
			return lufact(SparseMatrix)
		else
			error("Symbol: `", my_settings.linear_system_solver, "' does not correspond to an implemented linear equation solver. Check the settings.")
		end
	catch e
		println("ERROR in ls_factor")
		throw(e)
	end
end

function ls_solve!(my_factor, my_rhs, my_sol , my_settings::class_settings)
	try
		if my_settings.linear_system_solver in [:solver_MUMPS_LU, :solver_MUMPS_LDL]
			applyMUMPS!(my_factor, my_rhs, my_sol);
		elseif my_settings.linear_system_solver == :solver_julia
			my_sol[1:length(my_sol)] = my_factor \ my_rhs;
		else
			error("Symbol: `", my_settings.linear_system_solver, "' does not correspond to an implemented linear equation solver. Check the settings.")
		end
	catch e
		println("ERROR in ls_factor")
		throw(e)
	end
end
