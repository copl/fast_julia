type class_homogeneous_residuals <: abstract_residuals
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

  val_c::Float64

	primal_norm::Float64
	dual_norm::Float64

	primal_infeas_sign::Int64
	dual_infeas_sign::Int64
	primal_infeas_norm::Float64
	dual_infeas_norm::Float64

  function class_homogeneous_residuals()
      return new();
  end
end


function update_residuals!(res::class_homogeneous_residuals, nlp_eval::internal_AbstractNLPEvaluator, vars::class_variables, newt::abstract_newton_direction)
			try
        nlp_vals = newt.nlp_vals;
        update_nlp_cache!(nlp_vals, nlp_eval, vars)


				res.mu = mu(newt,vars);
        res.val_c = nlp_vals.val_c

				val_x_scaled = x_scaled(vars);
        val_s = s(vars);
        val_kappa = kappa(vars);
        val_tau = tau(vars);

				res.r_D = s(vars) - tau(vars) * nlp_vals.val_gradlag;
				res.r_G = kappa(vars) + dot(nlp_vals.val_gradlag, x(vars)) + dot(nlp_vals.val_a, y(vars));
				res.r_P = - tau(vars) * nlp_vals.val_a;

				res.r_D_norm = norm(res.r_D,1);
				res.r_G_norm = abs(res.r_G);
				res.r_P_norm = norm(res.r_P,1);

				# other stuff
				scale = kappa(vars) + tau(vars);
				res.r_D_norm_scaled = res.r_D_norm/scale;
				res.r_G_norm_scaled = res.r_G_norm/scale;
				res.r_P_norm_scaled = res.r_P_norm/scale;

				res.scaled_mu = res.mu/(kappa(vars) + tau(vars));

				res.primal_norm = norm(res.r_P,1)/tau(vars);
				res.dual_norm = norm(res.r_D,1)/tau(vars);

				primal_infeas_obj = -dot(nlp_vals.val_b, y(vars));
				dual_infeas_obj = -dot(x(vars), nlp_vals.val_gradc);
				res.primal_infeas_sign = sign(primal_infeas_obj);
				res.dual_infeas_sign = sign(dual_infeas_obj);

				res.primal_infeas_norm = norm(s(vars) + nlp_vals.val_jac_a' * y(vars),1)/abs(primal_infeas_obj);
 				res.dual_infeas_norm = norm(nlp_vals.val_jac_a * x(vars),1)/abs(dual_infeas_obj);

			catch e
				  println("ERROR in class_residuals.update")
				  throw(e)
			end
end


