type class_settings
	max_it::Int64	

	kappa_tau_tol::Float64
	gap_tol::Float64
	
	primal_tol::Float64
	dual_tol::Float64
	
	primal_infeas_tol::Float64
	dual_infeas_tol::Float64

	
	diagonal_modification::Float64
	linear_system_solver::Symbol
	# options
	# :solver_MUMPS_LDL
	# :solver_MUMPS_LU
	# :solver_julia	

	function class_settings()
		this = new();
	
		# termination criteron
		this.max_it = 100;

		this.kappa_tau_tol = 1e-4
		this.gap_tol = 1e-8;

		this.primal_tol = 1e-6;
		this.dual_tol = 1e-6;		
				
		this.primal_infeas_tol = 10-6;
		this.dual_infeas_tol = 10-6;

		# algorithm parameters
		this.diagonal_modification = 1e-8;
		this.linear_system_solver = :solver_julia;

		return this
	end
end
