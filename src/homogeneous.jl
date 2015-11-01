# make a module !!!
# move to homogeneous_algorithm file
# the ideas is that testing can call load code

function include_print(filename::String)
	include(filename)
	println("loaded ", filename)
end

function print_if(statement::String, print_this::Bool)
	if print_this
		println(statement)
	end
end


function display_progress(itr::Int64, alpha::Float64, gamma::Float64, residuals::class_homogeneous_residuals, vars::class_variables, direction::class_variables, delta::Float64, ls_num_trials::Int64, num_facs::Int64, settings::class_settings)
	try
		output = @sprintf("%s %2.1e %2.1e | %2.1e %2.1e | %2.1e %2.1e %2.1e %2.1e %2.1e | %2.1e %2.1e %s %s \n",
                      rpad(string(itr),3), alpha, gamma, tau(vars), kappa(vars), residuals.scaled_mu, residuals.r_G_norm, residuals.r_P_norm, residuals.r_D_norm, 0.0, delta, alpha * norm(direction._v,1), rpad(string(ls_num_trials),2), rpad(string(num_facs),2) )

		if settings.verbose
			print(output)
		end
	catch e
		println("ERROR in display_progress")
		throw(e)
	end
end



function terminate_algorithm(vars::class_variables, residuals::class_homogeneous_residuals, settings::class_settings)
	try
		kappa_tau_off = true
		if residuals.scaled_mu < settings.gap_tol
			if kappa(vars)/tau(vars) < settings.kappa_tau_tol || kappa_tau_off
				if residuals.primal_norm < settings.primal_tol && residuals.dual_norm < settings.dual_tol
					return 1;
				end
			end

			if tau(vars) / kappa(vars) < settings.kappa_tau_tol || kappa_tau_off
				#println( -(residuals.b' * vars.y())[1], " ", (vars.x()' * residuals.c)[1])
				#println("res:", residuals.primal_infeas_norm, " ", residuals.dual_infeas_norm)
				if residuals.primal_infeas_sign == 1 && residuals.primal_infeas_norm < settings.primal_infeas_tol
					return 2;
				elseif residuals.dual_infeas_sign == 1 && residuals.dual_infeas_norm < settings.dual_infeas_tol && residuals.func_val < -settings.unbounded_value
					return 3;
				end
			end
		end

		return 0
	catch e
		println("ERROR in terminate_algorithm")
		throw(e)
	end


end

function homogeneous_algorithm(qp::class_quadratic_program, vars::class_variables, settings::class_settings)
	alpha = 0.0;
	try
		#GLOBAL_timer::class_algorithm_timer
		#GLOBAL_timer.reset();
		#GLOBAL_timer.start("Total")
		#GLOBAL_timer.start("Intial")

		status = 0;

		validate_dimensions(qp,vars)
		#newton_solver = class_newton_solver2(qp, vars, settings);
		newton_solver = settings.newton_solver;

    initialize_newton!(newton_solver, qp, vars, settings);


		#GLOBAL_timer.stop("Intial")

		it = 0;
		gamma = 0.0;
		num_trials = 0;
		total_factorizations = 0;

		print_if("It | alpha | gamma  || tau   | kappa  ||  mu  |  gap  | primal | dual | f(x/tau)|| delta norm(d) ls# del#", settings.verbose)
		display_progress(it, alpha, gamma, newton_solver.residuals, vars, newton_solver.direction, newton_solver.delta, num_trials, 0, settings);

		for it = 1:settings.max_it
			#update_newton!(newton_solver, vars, settings);
      #used_delta = 0.0;
      #num_facs = 1;
      used_delta, num_facs = ipopt_style_inertia_correction!(newton_solver, vars, settings)
			total_factorizations += num_facs;
			#used_delta = newton_solver.update_system(vars)

      #vars, alpha, gamma = simple_gamma_strategy(newton_solver, vars, settings)
			vars, alpha, gamma = hybrid_mu_strategy(newton_solver, vars, settings, used_delta)

			update_residuals!(newton_solver.residuals, qp, vars, newton_solver);
			display_progress(it, alpha, gamma, newton_solver.residuals, vars, newton_solver.direction, used_delta, num_trials, num_facs, settings);

			status = terminate_algorithm(vars, newton_solver.residuals, settings);
			if status != 0
				print_if("Termination criteron met", settings.verbose)
				print_if("status = " * string(status), settings.verbose)
				break
			end


		end

		if status == 0
			print_if("MAXIMUM ITERATIONS REACHED", settings.verbose)
		end

		#GLOBAL_timer.stop("Total")

		#if settings.verbose
		#	GLOBAL_timer.print_stats("Total");
		#end

		println("Total factorizations: ", total_factorizations)

		return status, vars, it, total_factorizations
	catch e
		println("alpha = ", alpha)
		println("ERROR in homogeneous_algorithm")
		throw(e)
	end
end

