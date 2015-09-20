##############################################
# RESIDUALS
##############################################
type class_residuals
	a
	H
	h
	c
	b
	func_val::Float64
	
	r_G::Float64
	r_D::Array{Float64,1}
	r_P::Array{Float64,1}

	r_D_norm_scaled::Float64
	r_G_norm_scaled::Float64
	r_P_norm_scaled::Float64

	mu::Float64
	scaled_mu::Float64

	primal_norm::Float64
	dual_norm::Float64
	
	primal_infeas_sign::Int64
	dual_infeas_sign::Int64
	primal_infeas_norm::Float64
	dual_infeas_norm::Float64

	update::Function

	function class_residuals()
		this = new();
		
		this.update = function(qp::class_quadratic_program, vars::class_variables, newton_solver::abstract_newton_solver)
			try 
				this.mu = newton_solver.mu(vars);
	
				x_scaled = vars.x()/vars.tau();
				
				this.a = qp.A * x_scaled - qp.b;				
					
				this.H = qp.H(vars);
				this.h = this.H * x_scaled;
				this.c = qp.q + this.H * x_scaled;
				this.b = this.a - qp.A * x_scaled;
				this.func_val = qp.f(vars);

				g = this.c - qp.A' * vars.y() / vars.tau();
				
				this.r_D = vars.s() - vars.tau() * g;
				this.r_G = (vars.kappa() + g' * vars.x() + this.a' * vars.y())[1]; 
				this.r_P = - vars.tau() * this.a;


				# other stuff
				scale = vars.kappa() + vars.tau();
				this.r_D_norm_scaled = norm(this.r_D,1)/scale;
				this.r_G_norm_scaled = norm(this.r_G,1)/scale;
				this.r_P_norm_scaled = norm(this.r_P,1)/scale;

				this.scaled_mu = this.mu/(vars.kappa() + vars.tau());
				
				this.primal_norm = norm(this.r_P,1)/vars.tau();
				this.dual_norm = norm(this.r_D,1)/vars.tau();
				
				primal_infeas_obj = -(this.b' * vars.y())[1];
				dual_infeas_obj = -(vars.x()' * this.c)[1];
				this.primal_infeas_sign = sign(primal_infeas_obj)
				this.dual_infeas_sign = sign(dual_infeas_obj)
			
				this.primal_infeas_norm = norm(vars.s() + qp.A' * vars.y(),1)/abs(primal_infeas_obj);
 				this.dual_infeas_norm = norm(qp.A * vars.x(),1)/abs(dual_infeas_obj);
			catch e
				println("ERROR in class_residuals.update")
				throw(e)
			end
		end
		
		return this;	
	end

end

##############################################
# NEWTON SYSTEM SOLVER
##############################################

abstract abstract_newton_solver


# DELETE THIS IT IS OLD !!!
type class_newton_solver <: abstract_newton_solver
	direction::class_variables
	
	update_system::Function
	compute_direction::Function
	maximum_delta_decrease::Function

	function class_newton_solver(qp::class_quadratic_program, settings::class_settings)
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
				inertia = linear_system_solver.ls_factor(vars, qp);
				GLOBAL_timer.stop("Factor")

				compute_h_vector(vars) # move ???
				
				return inertia # is inertia is correct?
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
				
				# is this direction valid ?
				dir.check_direction();
			catch e
				println("ERROR in class_newton_solver.compute_direction")
				throw(e)
			end
		end

		return this;
	end
end




################################################################################################################

type class_newton_hsd <: abstract_newton_solver
	delta::Float64
	direction::class_variables
	residuals::class_residuals
	
	update_system::Function
	compute_direction::Function
	intialize::Function

	minimum_delta::Function
	mu::Function
	
	K::SparseMatrixCSC{Float64,Int64}
	K_true::SparseMatrixCSC{Float64,Int64} # temp

	linear_system_solver::abstract_linear_system_solver

	function class_newton_hsd()
		this = new();
	
		this.mu = function(vars::class_variables)
			return ( ( vars.x()' * vars.s() + vars.tau() * vars.kappa() ) / (vars.n() + 1) )[1]
		end	
		
		this.minimum_delta = function(qp::class_quadratic_program, cur_delta::Float64)
			lambda, eigen_vector = eig_min!(this.linear_system_solver, this.K, randn(qp.n + 1), qp.n + 1, qp.m, 1e-6);
			return (lambda - cur_delta), eigen_vector
		end	


		this.intialize = function(qp::class_quadratic_program, vars::class_variables, settings::class_settings)
			try
				dim = qp.n + qp.m + 1;
				this.K = spzeros(dim, dim);

				this.direction = class_variables(qp.n,qp.m);
			
				this.linear_system_solver = settings.linear_system_solver;
				this.linear_system_solver.initialize(this.K);

				this.residuals = class_residuals();
				this.residuals.update(qp, vars, this);

				this.delta = 0.0;
			catch e
				println("ERROR in class_newton_solver.intialize")
				throw(e)
			end
		end

		this.update_system = function(vars::class_variables, qp::class_quadratic_program)
			try
				res = this.residuals;
				x_scaled = vars.x() / vars.tau();
				#res.H += this.delta * speye(length(x_scaled))
				D_x = res.H + spdiagm( vars.s() ./ vars.x() ) + this.delta * speye(length(x_scaled));
				D_g = x_scaled' * res.H * x_scaled + vars.kappa() / vars.tau() + this.delta;				
				D_z = 0.1 * speye(qp.m,qp.m);
				v_1 = -res.c - res.H * x_scaled;
				v_2 = res.c - res.H * x_scaled;
				v_3 = res.a - qp.A * x_scaled;

				this.K_true = [ 
					[ D_x  	v_2 	-qp.A' 	]; 
					[ v_1' 	D_g 	-v_3' 	];  
					[ qp.A 	v_3	 D_z	] 
					];

				this.K[:,:] = [ 
					[ D_x  		res.h 	 qp.A' 	]; 
					[ res.h'	D_g 	 v_3' 	];  
					[ qp.A 		v_3	 -D_z	] 
					];

				#println(full(this.K))

				GLOBAL_timer::class_algorithm_timer
				GLOBAL_timer.start("Factor")
				inertia = this.linear_system_solver.ls_factor(qp.n + 1, qp.m);
				GLOBAL_timer.stop("Factor")
				
				return inertia
			catch e
				println("ERROR in class_newton_solver.update_system")
				throw(e)
			end
		end
		
		this.compute_direction = function(vars::class_variables, qp::class_quadratic_program, gamma::Float64, eta::Float64)
			try
				res = this.residuals;
				mu = res.mu;
				
				xs = -vars.x() .* vars.s() + gamma * mu * ones(qp.n);
				tk = -vars.tau() * vars.kappa() + gamma * mu;


				rhs = [	eta * res.r_D + xs ./ vars.x(); 
					eta * res.r_G + tk / vars.tau();
					eta * res.r_P ];

				sol = ones(qp.n + qp.m + 1);
				#linear_system_solver.ls_solve!(rhs,sol);
				fac = lufact(this.K_true)
				sol = fac \ rhs;
				
				dir = this.direction;
				dir.v[dir.x_ind] = sol[1:qp.n];
				dir.v[dir.tau_ind] = sol[qp.n + 1];
				dir.v[dir.y_ind] = sol[(qp.n + 2):(qp.m + qp.n + 1)]; # (dx,dy)

				dir.v[dir.s_ind] = (xs - vars.s() .* dir.x()) ./ vars.x();
				dir.v[dir.kappa_ind] = (tk - vars.kappa() * dir.tau()) / vars.tau();
				
				# is this direction valid ?
				dir.check_direction();
			catch e
				println("ERROR in class_newton_solver.compute_direction")
				throw(e)
			end
		end

		return this;
	end
end


################################################################################################################
# DIRECT IP METHOD
################################################################################################################


type class_newton_ip <: abstract_newton_solver
	delta::Float64
	direction::class_variables
	residuals::class_residuals
	K::SparseMatrixCSC{Float64,Int64}

	intialize::Function
	update_system::Function
	compute_direction::Function
	
	minimum_delta::Function
	mu::Function

	linear_system_solver::abstract_linear_system_solver


	function class_newton_ip()
		this = new();
		
		this.mu = function(vars::class_variables)
			return (  vars.x()' * vars.s() / vars.n()  )[1]	
		end

		this.minimum_delta = function(qp::class_quadratic_program, cur_delta::Float64)
			lambda, eigen_vector = eig_min!(this.linear_system_solver, this.K, randn(qp.n), qp.n, qp.m, 1e-6);
			#println("lambda = ", lambda)
			return (lambda - cur_delta), eigen_vector
		end

		this.intialize = function(qp::class_quadratic_program, vars::class_variables, settings::class_settings)
			try		
				vars.ipopt_start(); # set tau = 1, kappa = 0
		
				dim = qp.n + qp.m;
				this.K = spzeros(dim, dim);

				this.direction = class_variables(qp.n,qp.m);

				this.linear_system_solver = settings.linear_system_solver;
				this.linear_system_solver.initialize(this.K);

				this.residuals = class_residuals();
				this.residuals.update(qp, vars, this);

				this.delta = 0.0;
			catch e
				println("ERROR in class_newton_solver.intialize")
				throw(e)
			end
		end

		this.update_system = function(vars::class_variables, qp::class_quadratic_program)
			try
				res = this.residuals;
				x_scaled = vars.x() / vars.tau();
				res.H += this.delta * speye(length(x_scaled))
				D_x = res.H + spdiagm( vars.s() ./ vars.x() )
				D_z = 0.1 * speye(qp.m,qp.m);

				this.K[:,:] = [ 
					[ D_x  		qp.A' 	]; 
					[ qp.A 		-D_z	] 
					];


				GLOBAL_timer::class_algorithm_timer
				GLOBAL_timer.start("Factor")
				inertia = this.linear_system_solver.ls_factor(qp.n, qp.m);
				GLOBAL_timer.stop("Factor")
				
				return inertia
			catch e
				println("ERROR in class_newton_solver.update_system")
				throw(e)
			end
		end
		
		this.compute_direction = function(vars::class_variables, qp::class_quadratic_program, gamma::Float64, eta::Float64)
			try
				res = this.residuals;
				mu = res.mu;
				
				xs = -vars.x() .* vars.s() + gamma * mu * ones(qp.n);
				tk = -vars.tau() * vars.kappa() + gamma * mu;


				rhs = [	eta * res.r_D + xs ./ vars.x(); 
					eta * res.r_P ];

				sol = ones(qp.n + qp.m);
				this.linear_system_solver.ls_solve!(rhs,sol);
				#sol = this.K_true \ rhs;
				
				dir = this.direction;
				dir.v[dir.x_ind] = sol[1:qp.n];
				dir.v[dir.tau_ind] = 0.0;
				dir.v[dir.y_ind] = -sol[(qp.n + 1):(qp.m + qp.n)]; # (dx,dy)

				dir.v[dir.s_ind] = (xs - vars.s() .* dir.x()) ./ vars.x();
				dir.v[dir.kappa_ind] = 0.0;
				
				# is this direction valid ?
				dir.check_direction();
			catch e
				println("ERROR in class_newton_solver.compute_direction")
				throw(e)
			end
		end

		return this;
	end
end





