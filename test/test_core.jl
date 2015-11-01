# test class_variables
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
end

function qp_test()
    Q = speye(3);
    A = sparse(ones(2,3));
    b = 3*ones(2);
    c = ones(3);
    return class_quadratic_program(A, b, c, Q), A, Q, c;
end



# test quadratic program and validate dimensions
begin
    var = class_variables(3,2);
    qp, A, Q, c = qp_test();

    x_var = x(var);
    y_var = y(var);

    @test 4.5 == internal_eval_c(qp, x_var)
    @test [0.0, 0.0] == internal_eval_a(qp, x_var)
    @test A == internal_eval_jac_a(qp, x_var)
    @test Q == internal_eval_hesslag_prod(qp, x_var, y_var)
    @test c + Q * x_var + A' * y_var == internal_eval_gradlag(qp, x_var, y_var)

    validate_dimensions(qp, var)

    nlp_vals = class_nlp_cache();
    update_nlp_cache!(nlp_vals, qp, var)
end




