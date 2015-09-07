# TO DO:
# implement upper and lower bounds
# normalize

type class_quadratic_program
	n::Int64 # number of variables
	m::Int64 # number of constraints
	_H::SparseMatrixCSC{Float64,Int64}
	c::Array{Float64,1}
	A::SparseMatrixCSC{Float64,Int64}
	b::Array{Float64,1}

	a::Function	
	g::Function
	H::Function
	
	delta::Float64 # regularization

	function class_quadratic_program()
		this = new();
		
		#this.f = function(vars::class_variables)
		#	x_scaled = vars.x()/vars.tau();
		#	return (this.c'*x_scaled + 0.5*x_scaled'*this.H()*x_scaled)[1]
		#end

		this.a = function(vars::class_variables)
			return this.A*vars.x() - this.b;
		end

		this.g = function(vars::class_variables)
			return this.c + this.H(vars) * vars.x()/vars.tau()
		end
	
		this.H = function(vars::class_variables)
			return this._H + this.delta * speye(length(vars.x()));
		end

		return this;	
	end

	function class_quadratic_program(A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1},c::Array{Float64,1})
		this = class_quadratic_program();

		(this.m, this.n) = size(A)
		this.A = A;
		this._H = spzeros(this.n,this.n)
		this.c = c;
		this.b = b;

		return this;
	end

	function class_quadratic_program(A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1},c::Array{Float64,1},Q::SparseMatrixCSC{Float64,Int64})
		this = class_quadratic_program();

		(this.m, this.n) = size(A)
		this.A = A;
		this._H = Q;
		this.c = c;
		this.b = b;

		return this;
	end
end


type class_variables
	v::Array{Float64,1}
	
	x_ind::UnitRange{Int64}
	y_ind::UnitRange{Int64}
	s_ind::UnitRange{Int64}
	tau_ind::Int64
	kappa_ind::Int64

	check_positive::Function

	x::Function
	y::Function
	s::Function
	tau::Function
	kappa::Function

	mu::Function

	function class_variables(n::Int64, m::Int64)	
		this = new();
		# construct indicies
		this.x_ind = 1:n;
		this.y_ind = (n+1):(n+m);
		this.s_ind = (n+m+1):(2*n + m);
		this.tau_ind = 2*n + m + 1;
		this.kappa_ind = 2*n + m + 2;		

		# intialize
		this.v = zeros(2*(n+1)+m);
		this.v[this.x_ind] = 1.0;
		this.v[this.s_ind] = 1.0;
		this.v[this.tau_ind] = 1.0;
		this.v[this.kappa_ind] = 1.0;

		this.x = function()
			return this.v[this.x_ind]
		end

		this.y = function()
			return this.v[this.y_ind]
		end

		this.s = function()
			return this.v[this.s_ind]
		end

		this.tau = function()
			return this.v[this.tau_ind]
		end

		this.kappa = function()
			return this.v[this.kappa_ind]
		end

		this.mu = function()
			return ( ( this.x()' * this.s() + this.tau() * this.kappa() ) / (n + 1) )[1]
		end
		
		this.check_positive = function()
			for x_i in this.x()
				@assert(x_i > 0)		
			end

			for s_i in this.s()
				@assert(s_i > 0)		
			end

			@assert(this.tau() > 0)
			@assert(this.kappa() > 0)
		end

		return this;
	end
end

function validate_dimensions(qp::class_quadratic_program,vars::class_variables)
	try
		@assert(size(qp.H(vars)) == (qp.n,qp.n))
		@assert(size(qp.A) == (qp.m,qp.n))
		@assert(size(qp.b) == (qp.m,))
		@assert(size(qp.c) == (qp.n,))
		@assert(length(vars.x()) == qp.n)
		@assert(length(vars.y()) == qp.m)
		@assert(length(vars.v) == 2*(qp.n + 1) + qp.m) 
 	catch e
		println("ERROR validate_dimensions")
		throw(e)
	end
end

type class_residuals
	r_G::Float64
	r_D::Array{Float64,1}
	r_P::Array{Float64,1}

	r_D_norm_scaled::Float64
	r_G_norm_scaled::Float64
	r_P_norm_scaled::Float64

	scaled_mu::Float64

	primal_norm::Float64
	dual_norm::Float64
	
	primal_infeas_norm::Float64
	dual_infeas_norm::Float64

	update::Function

	function class_residuals()
		this = new();
		
		this.update = function(qp::class_quadratic_program, vars::class_variables)
			try
				g = qp.g(vars);
				this.r_D = vars.s() - qp.A' * vars.y() - vars.tau() * g;
				this.r_G = (vars.kappa() + vars.x()' * g + qp.b' * vars.y())[1]; 
				this.r_P = qp.b * vars.tau() - qp.A * vars.x();

				scale = vars.kappa() + vars.tau();
				this.r_D_norm_scaled = norm(this.r_D,1)/scale;
				this.r_G_norm_scaled = norm(this.r_G,1)/scale;
				this.r_P_norm_scaled = norm(this.r_P,1)/scale;

				this.scaled_mu = vars.mu()/(vars.kappa() + vars.tau());
				
				this.primal_norm = norm(this.r_P,1)/vars.tau();
				this.dual_norm = norm(this.r_D,1)/vars.tau();

				this.primal_infeas_norm = norm(vars.s() - qp.A' * vars.y(),1)/(qp.b' * vars.y())[1];
 				this.dual_infeas_norm = norm(qp.A * vars.x(),1)/(vars.x()' * g)[1];
			catch e
				println("ERROR in class_residuals.update")
				throw(e)
			end
		end
		
		return this;	
	end

end





