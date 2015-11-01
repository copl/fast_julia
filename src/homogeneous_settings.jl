abstract abstract_linear_system_solver # see linear_system_solver.jl
abstract abstract_newton_solver # see newton_solver.jl

type class_settings
	max_it::Int64

	kappa_tau_tol::Float64
	gap_tol::Float64

	primal_tol::Float64
	dual_tol::Float64

	primal_infeas_tol::Float64
	dual_infeas_tol::Float64
	unbounded_value::Float64

	diagonal_modification::Float64


	linear_system_solver::abstract_linear_system_solver
	newton_solver::abstract_newton_solver

	delta_min::Float64
	delta_max::Float64
	delta_start::Float64
	delta_increase::Float64
	delta_decrease::Float64

	verbose::Bool

	function class_settings()
		this = new();

		# termination criteron
		this.max_it = 1000;

		this.kappa_tau_tol = Inf; #1e-8
		this.gap_tol = 1e-8;

		this.primal_tol = 1e-8;
		this.dual_tol = 1e-8;

		this.primal_infeas_tol = 1e-8;
		this.dual_infeas_tol = 1e-8;
		this.unbounded_value = 1e4;

		# algorithm parameters
		this.diagonal_modification = 1e-6;

		# algorithm strategies
		this.linear_system_solver = linear_solver_MATLAB();
		#this.newton_solver = class_newton_hsd();
		this.newton_solver = class_newton_hsd();

		# delta parameters
		this.delta_min = 1e-8;
		this.delta_max = 1e20;
		this.delta_start = 1e-4;
		this.delta_increase = 8;
		this.delta_decrease = 1.0/3.0;

		this.verbose = true;

		return this
	end
end
