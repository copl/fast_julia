include("../src/core/core.jl")
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




# test quadratic program and validate dimensions
begin
    var = class_variables(3,2);


    Q = speye(3);
    A = sparse(ones(2,3));
    b = 3*ones(2);
    c = ones(3);
    qp = class_quadratic_program(A, b, c, Q);

    x_var = x(var);
    y_var = y(var);

    @test 4.5 == internal_eval_f(qp, x_var)
    @test [0.0, 0.0] == internal_eval_g(qp, x_var)
    @test A == internal_eval_jac_g(qp, x_var)
    @test Q == internal_eval_hesslag_prod(qp, x_var, y_var)
    @test c + Q * x_var + A' * y_var == internal_eval_gradlag(qp, x_var, y_var)

    validate_dimensions(qp, var)

end



