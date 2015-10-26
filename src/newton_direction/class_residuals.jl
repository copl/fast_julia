type class_residuals
	r_G::Float64
	r_D::Array{Float64,1}
	r_P::Array{Float64,1}

	r_D_norm::Float64
	r_G_norm::Float64
	r_P_norm::Float64

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
end

function update!(nlp_vals::nlp_cache, vars::class_variables)
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

				this.r_D_norm = norm(this.r_D,1);
				this.r_G_norm = abs(this.r_G);
				this.r_P_norm = norm(this.r_P,1);


				# other stuff
				scale = vars.kappa() + vars.tau();
				this.r_D_norm_scaled = this.r_D_norm/scale;
				this.r_G_norm_scaled = this.r_G_norm/scale;
				this.r_P_norm_scaled = this.r_P_norm/scale;

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
