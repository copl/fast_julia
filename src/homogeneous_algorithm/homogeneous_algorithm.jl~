# to do:
# check primal & dual termination criteron

function include_print(filename::String)
	include(filename)
	println("loaded ", filename)
end


include_print("homogeneous_settings.jl")
include_print("homogeneous_core.jl")
include_print("class_newton_solver.jl")
include_print("line_search.jl")
include_print("linear_system_solvers.jl")
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



type class_algorithm_timer
	
	start_times::Dict{String,Float64}
	total_times::Dict{String,Float64}

	start::Function
	stop::Function
	reset::Function
	print_stats::Function

	function class_algorithm_timer()
		this = new();
		this.start_times = Dict{String,Float64}()
		this.total_times = Dict{String,Float64}()

		function round_num(num,k=1)
			return round(num*10^k)/10^k
		end

		function show_percent(val,total)
			return round_num(100*val/total)
		end

		this.start = function(str)
			this.start_times[str] = time();
		end

		this.stop = function(str)
			if ~(str in keys(this.total_times))
				this.total_times[str] = 0.0;
			end
			time_passage = (time() - this.start_times[str]);
			this.total_times[str] = this.total_times[str] + time_passage;
		end

		this.reset = function()
			this.total_times = Dict{String,Float64}()
		end

		this.print_stats = function(standardize_by)
			println("======= Time statistics =======")
			if standardize_by != ""
				total = this.total_times[standardize_by]
				println( standardize_by, " = ", round_num(total, 3) )
				for label in keys(this.total_times)
					if label != standardize_by
						println( label, " = ", show_percent(this.total_times[label],total), " %" )
					end
				end
			else
				for label in keys(this.total_times)
					println( label, " = ", round_num(this.total_times[label], 3) )
				end
			end			
		end

		return this;
	end
end

global GLOBAL_timer = class_algorithm_timer()

#############################################################
# TO DO:
# - update solvers to allow inertia input
# - add normal equations solving option
# - add matlab LDL solver
#
# - update mpc_julia to exactly follow ipopt strategy
#############################################################

function inertia_correction(newton_solver::class_newton_solver, vars::class_variables, qp::class_quadratic_program)
	j = 0;	
	for j = 1:100
		inertia = newton_solver.update_system(vars)
		if inertia == 1
			qp.delta = min(qp.delta/3, 1e-8);
			break
		elseif inertia == 0
			if qp.delta == 1e-8;
				qp.delta = 1e-4;
			else
				qp.delta = qp.delta*8;
			end
		elseif inertia == -1
			error("numerical stability issues when computing KKT system !!!")
		else
			error("inertia_corection")
		end
	end
	if j == 100
		error("maximum iterations for inertia_corection reached")
	end
end

function predictor_corrector(newton_solver::class_newton_solver, vars::class_variables, residuals::class_residuals)
	# predictor
	gamma = 0.0;
	newton_solver.compute_direction(vars, residuals, gamma);
	predictor_vars, = line_search(vars, newton_solver.direction);
	 
	# corrector			
	gamma = (predictor_vars.mu()/vars.mu())^3;	
	newton_solver.compute_direction(vars, residuals, gamma);
	vars, alpha = line_search(vars, newton_solver.direction);

	return vars, alpha, gamma
end

function homogeneous_algorithm(qp::class_quadratic_program, vars::class_variables, settings::class_settings)
	alpha = 0.0;
	try	
		GLOBAL_timer.reset();
		GLOBAL_timer.start("Total")
		GLOBAL_timer.start("Intial")

		status = 0;
		
		validate_dimensions(qp,vars)
		newton_solver = class_newton_solver(qp,settings);
		residuals = class_residuals();
		residuals.update(qp, vars);

		GLOBAL_timer.stop("Intial")
		 

		println("It | alpha | gamma | tau  | kappa  |  mu  |  gap  | primal | dual | trial#")
		display_progress(0,0.0,0.0,residuals,vars,settings);
		
		for it = 1:settings.max_it
			inertia_correction(newton_solver, vars, qp)
			#newton_solver.update_system(vars)
			
			vars, alpha, gamma = predictor_corrector(newton_solver, vars, residuals)
			
			residuals.update(qp, vars);
			display_progress(it, alpha, gamma, residuals, vars, settings);

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
		
		GLOBAL_timer.stop("Total")
		
		GLOBAL_timer.print_stats("Total");
	catch e
		println("alpha = ", alpha)
		println("ERROR in homogeneous_algorithm")
		throw(e)
	end
end

