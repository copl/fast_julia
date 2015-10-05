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
				
				#toler = 1e-6
				#eig_min, eigen_vector, err = newton_solver.minimum_delta(qp, newton_solver.delta, toler);
				
				#println(full(newton_solver.K))
				#suggested_delta = min(1.000,newton_solver.delta * 0.999) #0.5*newton_solver.delta + 0.5*smallest_delta;
				
				#if err < toler
				#newton_solver.delta = max(-eig_min * 1.001, settings.delta_min)
				#else
				newton_solver.delta = newton_solver.delta * settings.delta_decrease;
				#end
				
				break
			elseif inertia == 0
				#println(full(newton_solver.K))
				
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
