function predictor_corrector(newton_solver::abstract_newton_direction, vars::class_variables, settings::class_settings)
	try
		# predictor
		gamma = 0.0;
		compute_newton_direction!(newton_solver, vars, gamma, 1 - gamma); # eta_P, eta_D, eta_G ???

		predictor_vars, = line_search(vars, newton_solver.direction);

		intial_mu = mu(newton_solver,vars);
		predictor_mu = mu(newton_solver,predictor_vars);

		# corrector
		reduction = predictor_mu/intial_mu;
		gamma = (reduction)^3 #*0.8;
		compute_newton_direction!(newton_solver, vars, gamma, 1 - gamma);
		vars, alpha = line_search(vars, newton_solver.direction);

		return vars, alpha, gamma
	catch e
		println("ERROR in predictor_corrector")
		throw(e)
	end
end

function simple_gamma_strategy(newton_solver::abstract_newton_direction, vars::class_variables, settings::class_settings)
	try
		gamma = 0.8;
		compute_newton_direction!(newton_solver, vars, gamma, 1.0);
		vars, alpha = line_search(vars, newton_solver.direction);

		return vars, alpha, gamma
	catch e
		println("ERROR in simple_gamma_strategy")
		throw(e)
	end
end

function hybrid_mu_strategy(newton_solver::abstract_newton_direction, vars::class_variables, settings::class_settings, used_delta::Float64)
	try
		if used_delta > settings.delta_min
			vars, alpha, gamma = simple_gamma_strategy(newton_solver, vars, settings)
		else
			vars, alpha, gamma = predictor_corrector(newton_solver, vars, settings)
		end
		return vars, alpha, gamma
	catch e
		println("ERROR in hybrid_mu_strategy")
		throw(e)
	end
end



function residuals_gamma_strategy(newton_solver::abstract_newton_direction, vars::class_variables, qp::class_quadratic_program)
	# try to balance residuals
	### not implemented
end
