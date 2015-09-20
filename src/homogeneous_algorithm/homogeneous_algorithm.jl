# to do:
# check primal & dual termination criteron

function include_print(filename::String)
	include(filename)
	println("loaded ", filename)
end


include_print("homogeneous_settings.jl")
include_print("homogeneous_core.jl")
include_print("newton_solver.jl")
include_print("line_search.jl")
include_print("linear_system_solvers.jl")
include_print("eigenvalues.jl")
include_print("timer.jl")
println("starting code ...")

function print_if(statement::String, print_this::Bool)
	if print_this
		println(statement)
	end
end

function display_progress(itr::Int64, alpha::Float64, gamma::Float64, residuals::class_residuals, vars::class_variables, direction::class_variables, delta::Float64, num_trials::Int64, settings::class_settings)
	try
		output = @sprintf("%s %2.1e %2.1e | %2.1e %2.1e | %2.1e %2.1e %2.1e %2.1e %2.1e | %2.1e %2.1e %s \n", rpad(string(itr),3), alpha, gamma, vars.tau(), vars.kappa(), residuals.scaled_mu, residuals.r_G_norm_scaled, residuals.r_P_norm_scaled, residuals.r_D_norm_scaled, residuals.func_val, delta, norm(direction.v,1), rpad(string(num_trials),2) )

		if settings.verbose
			print(output)
		end
	catch e
		println("ERROR in display_progress")
		throw(e)
	end
end



function terminate_algorithm(vars::class_variables, residuals::class_residuals, settings::class_settings)
	try
		kappa_tau_off = true
		if residuals.scaled_mu < settings.gap_tol	
			if vars.kappa()/vars.tau() < settings.kappa_tau_tol || kappa_tau_off
				if residuals.primal_norm < settings.primal_tol && residuals.dual_norm < settings.dual_tol
					return 1;
				end
			end

			if vars.tau()/vars.kappa() < settings.kappa_tau_tol || kappa_tau_off
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


#############################################################
# TO DO:
# - update solvers to allow inertia input
# - add normal equations solving option
# - put strategies into seperate file
#############################################################


#############################################################
# DIFFERENT STRATEGIES
##############################################################

function inertia_correction(newton_solver::abstract_newton_solver, vars::class_variables, qp::class_quadratic_program, settings::class_settings)
	try

		MAX_IT = 50;
		j = 0;	
		used_delta = 0.0;

		newton_solver.delta = max(newton_solver.delta, settings.delta_min );
		
		for j = 1:MAX_IT
			inertia = newton_solver.update_system(vars, qp)
			
			if inertia == 1
				used_delta = newton_solver.delta;
				
				smallest_delta, eigen_vector = newton_solver.minimum_delta(qp, newton_solver.delta);
				#println(smallest_delta)
				#println(full(newton_solver.K))
				#suggested_delta = min(1.000,newton_solver.delta * 0.999) #0.5*newton_solver.delta + 0.5*smallest_delta;

				newton_solver.delta = newton_solver.delta * settings.delta_decrease;
				#newton_solver.delta = max(suggested_delta, settings.delta_min)
				break
			elseif inertia == 0
				println(full(newton_solver.K))
				
				if newton_solver.delta <= settings.delta_min;
					newton_solver.delta = settings.delta_start;
				else
					newton_solver.delta = newton_solver.delta * settings.delta_increase;
				end
			elseif inertia == -1
				error("numerical stability issues when computing KKT system !!!")
			else
				error("inertia_corection")
			end
		end
		
		if j == MAX_IT
			error("maximum iterations for inertia_corection reached")
		else
			return used_delta
		end

		
	catch e
		println("ERROR in inertia_correction")
		throw(e)
	end
end

function predictor_corrector(newton_solver::abstract_newton_solver, vars::class_variables, qp::class_quadratic_program)
	try
		# predictor
		gamma = 0.0;
		newton_solver.compute_direction(vars, qp, gamma, 1 - gamma);
		predictor_vars, = line_search(vars, newton_solver.direction);

		intial_mu = newton_solver.mu(vars);
		predictor_mu = newton_solver.mu(predictor_vars);
		 
		# corrector			
		gamma = (predictor_mu/intial_mu)^3;	
		newton_solver.compute_direction(vars, qp, gamma, 1 - gamma);
		vars, alpha = line_search(vars, newton_solver.direction);

		return vars, alpha, gamma
	catch e
		println("ERROR in predictor_corrector")
		throw(e)
	end
end

function simple_gamma_strategy(newton_solver::abstract_newton_solver, vars::class_variables, qp::class_quadratic_program)
	try
		gamma = 0.9;	
		newton_solver.compute_direction(vars, qp, gamma, 1.0);
		vars, alpha = line_search(vars, newton_solver.direction);

		return vars, alpha, gamma
	catch e
		println("ERROR in simple_gamma_strategy")
		throw(e)
	end
end

function homogeneous_algorithm(qp::class_quadratic_program, vars::class_variables, settings::class_settings)
	alpha = 0.0;
	try	
		GLOBAL_timer::class_algorithm_timer
		GLOBAL_timer.reset();
		GLOBAL_timer.start("Total")
		GLOBAL_timer.start("Intial")
		
		status = 0;
		
		validate_dimensions(qp,vars)
		#newton_solver = class_newton_solver2(qp, vars, settings);
		newton_solver = settings.newton_solver;
		newton_solver.intialize(qp, vars, settings);
		

		GLOBAL_timer.stop("Intial")
		 
		it = 0;
		gamma = 0.0;
		num_trials = 0;

		print_if("It | alpha | gamma  || tau   | kappa  ||  mu  |  gap  | primal | dual | f(x/tau) || delta norm(d) trial#", settings.verbose)
		display_progress(it, alpha, gamma, newton_solver.residuals, vars, newton_solver.direction, newton_solver.delta, num_trials, settings);
		
		for it = 1:settings.max_it
			used_delta = inertia_correction(newton_solver, vars, qp, settings)
			#used_delta = newton_solver.update_system(vars)
			
			vars, alpha, gamma = predictor_corrector(newton_solver, vars, qp)
			#vars, alpha, gamma = simple_gamma_strategy(newton_solver, vars, qp)
			
			newton_solver.residuals.update(qp, vars, newton_solver);
			
			display_progress(it, alpha, gamma, newton_solver.residuals, vars, newton_solver.direction, used_delta, num_trials, settings);

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
		
		GLOBAL_timer.stop("Total")
		
		if settings.verbose
			GLOBAL_timer.print_stats("Total");
		end

		return status, vars, it
	catch e
		println("alpha = ", alpha)
		println("ERROR in homogeneous_algorithm")
		throw(e)
	end
end

