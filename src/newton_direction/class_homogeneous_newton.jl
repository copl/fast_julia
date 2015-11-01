type class_homogeneous_newton <: abstract_newton_direction
	delta::Float64
	direction::class_variables
	residuals::class_homogeneous_residuals
  nlp_vals::class_nlp_cache

	K::SparseMatrixCSC{Float64,Int64}
  K_true::SparseMatrixCSC{Float64,Int64}

	linear_system_solver::abstract_linear_system_solver

  function class_homogeneous_newton()
      return new();
  end
end

###################
# methods
###################

function initialize_newton!(newt::class_homogeneous_newton, nlp_eval::internal_AbstractNLPEvaluator, vars::class_variables, settings::class_settings)
      try
				dim = n(vars) + m(vars) + 1;
				newt.K = spzeros(dim, dim);

				newt.direction = class_variables(n(vars),m(vars));

				newt.linear_system_solver = settings.linear_system_solver;
				initialize!(newt.linear_system_solver, newt.K);

        newt.nlp_vals = class_nlp_cache();

				newt.residuals = class_homogeneous_residuals();
        update_residuals!(newt.residuals, nlp_eval, vars, newt);

        newt.delta = 0.0;
			catch e
				println("ERROR in class_homogeneous_newton.intialize")
				throw(e)
			end
end


function update_newton!(newt::class_homogeneous_newton, vars::class_variables, settings::class_settings)
    try
      nlp_vals = newt.nlp_vals;

      val_x_scaled = x_scaled(vars);
      H = nlp_vals.val_hesslag_prod;
      h = H * val_x_scaled;
      c = nlp_vals.val_gradc;
      A = nlp_vals.val_jac_a;
      b = nlp_vals.val_b

      #res.H += newt.delta * speye(length(val_x_scaled))
      D_x = H + spdiagm( s(vars) ./ x(vars) ) + newt.delta * speye( n(vars) );
      D_g = val_x_scaled' * ( H + newt.delta * speye( n(vars) ) ) * val_x_scaled + kappa(vars) / tau(vars);
      #D_g = val_x_scaled' * nlp_vals.val_hesslag_prod * val_x_scaled + vars.kappa() / vars.tau() + this.delta;
      D_z = settings.diagonal_modification * speye(m(vars), m(vars));
      v_1 = -c - h;
      v_2 = c - h;
      #v_3 = a - A * val_x_scaled;

      newt.K_true = [
        [ D_x  	v_2 	-A' 	];
        [ v_1' 	D_g 	b' 	];
        [ A 	 -b	   D_z	]
        ];

      newt.K[:,:] = [
        [ D_x  		-h 	   A' 	];
        [ -h'	    D_g 	 -b' 	];
        [ A 		-b	   -D_z	]
        ];

      #GLOBAL_timer::class_algorithm_timer
      #GLOBAL_timer.start("Factor")
      inertia = ls_factor(newt.linear_system_solver,n(vars) + 1, m(vars));
      #GLOBAL_timer.stop("Factor")

      return inertia
  catch e
      println("ERROR in class_homogeneous_newton.update_newton!")
      throw(e)
  end
end

function compute_newton_direction!(newt::class_homogeneous_newton, vars::class_variables, gamma::Float64, eta::Float64)
			try
        num_vars = n(vars);
        num_const = m(vars);
        res = newt.residuals;

				mu = res.mu;

				xs = -x(vars) .* s(vars) + gamma * mu * ones(num_vars);
				tk = -tau(vars) * kappa(vars) + gamma * mu;

				rhs =
        [	eta * res.r_D + xs ./ x(vars);
				eta * res.r_G + tk / tau(vars);
				eta * res.r_P ];

				sol = ones(num_vars + num_const + 1);
				#linear_system_solver.ls_solve!(rhs,sol);
				fac = lufact(newt.K_true)
				sol = fac \ rhs;

				dir = newt.direction;
				x(dir, sol[1:num_vars] );
				tau(dir, sol[ num_vars + 1] );
				y(dir, sol[( num_vars + 2):( num_const +  num_vars + 1)] );

				s(dir,(xs - s(vars) .* x(dir)) ./ x(vars) );
				kappa(dir, (tk - kappa(vars) * tau(dir)) / tau(vars) );

				# is this direction valid ?
				check_for_wrong_vals(dir);
			catch e
				println("ERROR in class_newton_solver.compute_direction")
				throw(e)
			end
		end

function mu(newton_direction::class_homogeneous_newton, vars::class_variables)
    return (dot(s(vars),x(vars)) + dot(tau(vars),kappa(vars)))/(n(vars) + 1)
end

