res = class_homogeneous_residuals();
newt = class_homogeneous_newton();
vars = class_variables(3,2);
mu(newt, vars);

settings = class_settings();
qp, A, b, c = qp_test();

validate_dimensions(qp, vars)
initialize_newton!(newt, qp, vars, settings)
