
# test linear system solvers

function test_linear_solver(solver::abstract_linear_system_solver)
   mat = speye(10);
   initialize!(solver,mat);

   mat[1,1] = 2;
   @test solver._SparseMatrix == mat

   inertia = ls_factor(solver, 10, 0)
   @test inertia == 1

   sol = 1.0*zeros(10);
   rhs = 1.0*ones(10);
   ls_solve!(solver,rhs,sol);

   @test sol == [0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
end


# test julia solver
begin
  # test julia lu factor
  ls_solver_julia = linear_solver_JULIA();
  ls_solver_julia.sym = 0;
  test_linear_solver( ls_solver_julia )

  # test julia cholesky factor
  ls_solver_julia = linear_solver_JULIA();
  ls_solver_julia.sym = 1;
  test_linear_solver( ls_solver_julia )
end

# test matlab ldl solver
begin
  ls_solver_matlab = linear_solver_MATLAB();
  ls_solver_matlab.sym = 0;
  test_linear_solver( ls_solver_matlab )
end
