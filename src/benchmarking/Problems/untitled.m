A = Problem.A;
%M = Problem.A;
M =  A * A';
M = M + speye(size(M,1)) * 1e-8;

tic()
[R,p,s] = chol(M,'vector');
factor_time = toc()

rhs = rand(size(M,1),1);

tic()
x = R \ (R \ rhs);
solve_time = toc()

factor_time/solve_time