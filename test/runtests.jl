# Two levels of tests
# Tests of individual functions
# Tests of code to solve simple unbounded, bounded problems etc ...

 using Base.Test
begin
  # includes
  include("../src/loadcode.jl");


  function runtests()
    begin
      include("test_linear_system_solvers.jl")

      ###############################
      # test core
      ###############################

      ###############################
      # test class_variables
      ###############################
      begin


        var = class_variables(3,3);

        x_values = [-1, 2, 2];
        var._v[var._x_ind] = x_values;

        @test x(var) == x_values
        @test_throws ErrorException check_positive(var)

        x_values = [1, NaN, 2];
        var._v[var._x_ind] = x_values;
        @test_throws ErrorException check_for_wrong_vals(var)

        x_values = [1, Inf, 2];
        var._v[var._x_ind] = x_values;
        @test_throws ErrorException check_for_wrong_vals(var)

        x_values = [1, -Inf, 2];
        var._v[var._x_ind] = x_values;
        @test_throws ErrorException check_for_wrong_vals(var)

        var = class_variables(3,3);
        @test_throws ErrorException move!(var,-2.0,var)
      end

      function qp_test()
          Q = speye(3);
          A = sparse(ones(1,3));
          b = 3*ones(1);
          c = ones(3);
          return class_quadratic_program(A, b, c, Q), A, b, c, Q;
      end

      settings = class_settings();
      qp, A, b, c, Q = qp_test();
      var = class_variables(3,1);

      # test quadratic program and validate dimensions
      begin
          x_var = x(var);
          y_var = y(var);

          @test 4.5 == internal_eval_c(qp, x_var)
          @test [0.0] == internal_eval_a(qp, x_var)
          @test A == internal_eval_jac_a(qp, x_var)
          @test Q == internal_eval_hesslag_prod(qp, x_var, y_var)
          @test c + Q * x_var + A' * y_var == internal_eval_gradlag(qp, x_var, y_var)

          vars2 = class_variables(3,3);
          @test_throws ErrorException validate_dimensions(qp, vars2)
      end

      vars = class_variables(3,1);
      nlp_vals = class_nlp_cache();
      update_nlp_cache!(nlp_vals, qp, var)
      validate_dimensions(qp, var)

        # test newton direction


        begin
            newt = class_homogeneous_newton();
            initialize_newton!(newt, qp, vars, settings)
            res = class_homogeneous_residuals();

        end


        begin
            newt = class_homogeneous_newton();

            @test 1.0 == mu(newt, vars);

            validate_dimensions(qp, vars)
            initialize_newton!(newt, qp, vars, settings)
            update_newton!(newt, vars, settings)
            compute_newton_direction!(newt, vars, 0.5, 0.5)
        end

        # test line search
        begin
            vars2= class_variables(5,3);
            direction = class_variables(5,3);
            vars2, alpha = line_search(vars2, direction);
            @test alpha < 1
        end

        # test stratagies

        # test delta heuristics
        begin
            #println("TEST IPOPT INERTIA CORRECTION")
            validate_dimensions(qp, vars)
            delta, number_of_factors = ipopt_style_inertia_correction!(newt, vars, settings)
            @test delta == settings.delta_min
            @test number_of_factors == 1
        end

         # test theta heursitics
        begin
            predictor_corrector(newt, vars, settings)
            simple_gamma_strategy(newt, vars, settings)
            hybrid_mu_strategy(newt, vars, settings, settings.delta_min + 1.0)
            hybrid_mu_strategy(newt, vars, settings, settings.delta_min)
        end




        begin
            vars = class_variables(3,1);
            A = sparse(ones(1,3));
            b = 3*ones(1);
            c = ones(3);
            qp = class_quadratic_program(A, b, c);
            status, vars = homogeneous_algorithm(qp, vars, settings)
            println(x(vars)/tau(vars))
        end


    end
  end
  runtests();
end



