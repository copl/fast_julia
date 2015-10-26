type class_homogeneous_direction <: abstract_newton_direction
	delta::Float64
	direction::class_variables
	residuals::class_residuals

	K::SparseMatrixCSC{Float64,Int64}

	linear_system_solver::abstract_linear_system_solver
end

function initialize!(newton_direction::class_homogeneous_direction,...)

end


function update_system!(vars::class_variables, qp::class_quadratic_program, settings::class_settings)

end



