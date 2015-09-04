type class_newton_solver
	direction::class_variables
	
	update_system::Function
	compute_direction::Function

	function class_newton_solver(qp::class_quadratic_program,settings::class_settings)
		this = new();
		
		linear_system_solver = settings.linear_system_solver;
		linear_system_solver.initialize(qp, settings);
		
		this.direction = class_variables(qp.n,qp.m);
		rhs = zeros(qp.n + qp.m);
		h = zeros(qp.n + qp.m);
		p = zeros(qp.n + qp.m);

		function compute_h_vector(vars::class_variables)
			try
				rhs[1:(qp.n)] = qp.g(vars) - qp.H(vars) * vars.x()/vars.tau();
				rhs[(qp.n + 1):(qp.m + qp.n)] = -qp.b;

				GLOBAL_timer::class_algorithm_timer
				GLOBAL_timer.start("Solve")
				linear_system_solver.ls_solve!(rhs, h, qp);
				GLOBAL_timer.stop("Solve")
			catch e
				println("ERROR in class_newton_solver.compute_h_vector")
				throw(e)
			end
		end

		this.update_system = function(vars::class_variables)
			try
				GLOBAL_timer::class_algorithm_timer
				GLOBAL_timer.start("Factor")
				linear_system_solver.ls_factor(vars, qp);
				GLOBAL_timer.stop("Factor")

				compute_h_vector(vars)

				
				return 1 # inertia is correct
			catch e
				println("ERROR in class_newton_solver.update_system")
				throw(e)
			end
		end
		
		this.compute_direction = function(vars::class_variables, residuals::class_residuals, gamma::Float64)
			try			
				eta = 1 - gamma;
			
				r_G = residuals.r_G;
				r_D = residuals.r_D;
				r_P = residuals.r_P;
				x = vars.x();
				s = vars.s();	
				tau = vars.tau();
				kappa = vars.kappa();
				mu = vars.mu();
				g = qp.g(vars);

				xs = -x .* s + gamma * mu * ones(qp.n);
				tk = -tau*kappa + gamma * mu;

				rhs[1:(qp.n)] = eta * r_D + xs ./ x;
				rhs[(qp.n + 1):(qp.m + qp.n)] = eta * r_P;
				
				GLOBAL_timer::class_algorithm_timer
				GLOBAL_timer.start("Solve")
				linear_system_solver.ls_solve!(rhs, p, qp);
				GLOBAL_timer.stop("Solve")

				v_tmp = [g' + x'/tau * qp.H(vars) qp.b'];
				numerator_dtau = eta*r_G + 1/tau * tk + (v_tmp * p)[1];
				denominator_dtau = kappa/tau + (1/tau)^2 * (x' * qp.H(vars) * x)[1] + (v_tmp * h)[1];
				dtau = numerator_dtau/denominator_dtau;

				dir = this.direction;
				dir.v[dir.tau_ind] = dtau;
				dir.v[1:(qp.n + qp.m)] = p - h * dtau; # (dx,dy)

				dir.v[dir.s_ind] = (xs - s.*dir.x()) ./ x;
				dir.v[dir.kappa_ind] = (tk - kappa*dir.tau()) / tau;
			catch e
				println("ERROR in class_newton_solver.compute_direction")
				throw(e)
			end
		end

		return this;
	end
end

