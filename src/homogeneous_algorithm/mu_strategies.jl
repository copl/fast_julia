function predictor_corrector(newton_solver::abstract_newton_solver, vars::class_variables, qp::class_quadratic_program, settings::class_settings)
	try
		# predictor
		gamma = 0.0;
		newton_solver.compute_direction(vars, qp, gamma, 1 - gamma); # eta_P, eta_D, eta_G ???
		predictor_vars, = line_search(vars, newton_solver.direction);

		intial_mu = newton_solver.mu(vars);
		predictor_mu = newton_solver.mu(predictor_vars);
		 
		# corrector
		reduction = predictor_mu/intial_mu;			
		gamma = (reduction)^3 #*0.8;	
		newton_solver.compute_direction(vars, qp, gamma, 1 - gamma);
		vars, alpha = line_search(vars, newton_solver.direction);

		return vars, alpha, gamma
	catch e
		println("ERROR in predictor_corrector")
		throw(e)
	end
end

function simple_gamma_strategy(newton_solver::abstract_newton_solver, vars::class_variables, qp::class_quadratic_program, settings::class_settings)
	try
		gamma = 0.8;	
		newton_solver.compute_direction(vars, qp, gamma, 1.0);
		vars, alpha = line_search(vars, newton_solver.direction);

		return vars, alpha, gamma
	catch e
		println("ERROR in simple_gamma_strategy")
		throw(e)
	end
end

function hybrid_mu_strategy(newton_solver::abstract_newton_solver, vars::class_variables, qp::class_quadratic_program, settings::class_settings, used_delta::Float64)
	try
		if used_delta > settings.delta_min
			vars, alpha, gamma = simple_gamma_strategy(newton_solver, vars, qp, settings)
		else
			vars, alpha, gamma = predictor_corrector(newton_solver, vars, qp, settings)
		end

		return vars, alpha, gamma
	catch e
		println("ERROR in hybrid_mu_strategy")
		throw(e)
	end
end



function residuals_gamma_strategy(newton_solver::abstract_newton_solver, vars::class_variables, qp::class_quadratic_program)
	# try to balance residuals	
	### not implemented
end
