function line_search(vars::class_variables, direction::class_variables)
	try
		GLOBAL_timer::class_algorithm_timer
		GLOBAL_timer.start("Line-search")
		
		vars = deepcopy(vars)
		alpha_max = maximum_step(vars, direction);
		alpha = 0.9*min(1.0,0.9*alpha_max);	
		vars.v += alpha*direction.v;
	
		vars.check_positive();

		GLOBAL_timer.stop("Line-search")

		return vars, alpha
	catch e
		println("ERROR in line_search")
		throw(e)
	end
end

function maximum_step(vars::class_variables, direction::class_variables)
	try
		alpha = Inf;
		alpha = maximum_step(alpha, vars.x(), direction.x());
		alpha = maximum_step(alpha, vars.s(), direction.s());
		alpha = maximum_step(alpha, vars.tau(), direction.tau());
		alpha = maximum_step(alpha, vars.kappa(), direction.kappa());
		
		return alpha;
	catch e
		println("ERROR in maximum_step")
		throw(e)
	end
end

function maximum_step(alpha::Float64, array_point::Array{Float64}, array_direction::Array{Float64})
	try
		@assert(length(array_point) == length(array_direction))	
		for i in 1:length(array_point)
			alpha = maximum_step(alpha, array_point[i], array_direction[i]);
		end
	
		return alpha
	catch e
		println("ERROR in maximum_step(::Float64,::Array,::Array)")
		throw(e)
	end
end

function maximum_step(alpha::Float64, var::Float64, dir::Float64)
	try
		candidate_alpha = -var/dir;
		if candidate_alpha >= 0
			alpha = min(alpha, candidate_alpha)
		end

		return alpha
	catch e
		println("ERROR in maximum_step(::Float64,::Float64,::Float64)")
		throw(e)
	end
end
