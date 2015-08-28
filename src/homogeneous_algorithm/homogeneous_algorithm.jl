# to do:
# check primal & dual termination criteron

function include_print(filename::String)
	include(filename)
	println("loaded ", filename)
end


include_print("homogeneous_settings.jl")
include_print("linear_system_solvers.jl")
include_print("homogeneous_core.jl")
include_print("testing_tools.jl")
println("starting code ...")

function display_progress(itr::Int64, alpha::Float64, gamma::Float64, residuals::class_residuals, vars::class_variables, settings)
	try
		@printf("%s %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %i\n", rpad(string(itr),3), alpha, gamma, vars.tau(), vars.kappa(), residuals.scaled_mu, residuals.r_G_norm_scaled, residuals.r_P_norm_scaled, residuals.r_D_norm_scaled, 2)
	catch e
		println("ERROR in display_progress")
		throw(e)
	end
end



function terminate_algorithm(vars::class_variables,residuals::class_residuals,settings::class_settings)
	if residuals.scaled_mu < settings.gap_tol	
		if vars.kappa()/vars.tau() < settings.kappa_tau_tol
			if residuals.primal_norm < settings.primal_tol && residuals.dual_norm < settings.dual_tol
				return 1;
			end
		elseif vars.tau()/vars.kappa() < settings.kappa_tau_tol
			if residuals.primal_infeas_norm < settings.primal_infeas_tol
				return 2;
			elseif residuals.dual_infeas_norm < settings.dual_infeas_tol
				return 3;
			end
		end
	end

	return 0
end

function round_num(num,k=1)
	return round(num*10^k)/10^k
end

type class_algorithm_timer
	total::Float64
	intial::Float64
	factor::Float64
	solve::Float64
	ls::Float64
	res::Float64
	disp::Float64

	show_percent::Function

	function class_algorithm_timer()
		this = new();
		
		this.total = 0.0;
		this.intial = 0.0;
		this.factor = 0.0;
		this.solve = 0.0;
		this.res = 0.0;
		this.disp = 0.0;

		this.show_percent = function(val)
			return round_num(100*val/this.total)
		end

		return this;
	end
end

function homogeneous_algorithm(qp::class_quadratic_program, vars::class_variables, settings::class_settings)
	alpha = 0.0;
	try	
		mytimer = class_algorithm_timer();
		
		tic()
		
		status = 0;
		
		validate_dimensions(qp,vars)
		newton_solver = class_newton_solver(qp,settings);
		residuals = class_residuals();
		residuals.update(qp, vars);

		mytimer.intial = toq()
		
		tic()

		println("It | alpha | gamma | tau  | kappa  |  mu  |  gap  | primal | dual | trial#")
		display_progress(0,0.0,0.0,residuals,vars,settings);
		
		for it = 1:settings.max_it
			mytimer.factor += @elapsed newton_solver.update_system(vars)
			
			# predictor
			gamma = 0.0;
			mytimer.solve += @elapsed newton_solver.compute_direction(vars, residuals, gamma);
			mytimer.ls += @elapsed predictor_vars, = line_search(vars, newton_solver.direction);
			 
			# corrector			
			gamma = (predictor_vars.mu()/vars.mu())^3;	
			mytimer.solve += @elapsed newton_solver.compute_direction(vars, residuals, gamma);
			mytimer.ls += @elapsed vars, alpha = line_search(vars, newton_solver.direction);
			
			mytimer.res += @elapsed residuals.update(qp, vars);
			mytimer.disp += @elapsed  display_progress(it,alpha,gamma,residuals,vars,settings);

			status = terminate_algorithm(vars,residuals,settings);
			if status != 0
				println("Termination criteron met")
				println("status = ", status)
				break
			end
		end

		if status == 0
			println("MAXIMUM ITERATIONS REACHED")
		end
		
		mytimer.total = toq() + mytimer.intial;

		println("======= Time statistics =======")
		println("Total time = ", round_num(mytimer.total,3))		
		println("Intialization time = ",mytimer.show_percent(mytimer.intial),"%")
		println("Factorization time = ",mytimer.show_percent(mytimer.factor),"%")
		println("Solve time = ", mytimer.show_percent(mytimer.solve),"%")
		println("Line search = ", mytimer.show_percent(mytimer.ls),"%")
		println("Residual computation = ", mytimer.show_percent(mytimer.res),"%")
		println("Display time = ", mytimer.show_percent(mytimer.disp),"%")
	catch e
		println("alpha = ", alpha)
		println("ERROR in homogeneous_algorithm")
		throw(e)
	end
end

