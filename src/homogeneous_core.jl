type class_quadratic_program
	n::Int64 # number of variables
	m::Int64 # number of constraints
	_H::SparseMatrixCSC{Float64,Int64}
	q::Array{Float64,1}
	A::SparseMatrixCSC{Float64,Int64}
	b::Array{Float64,1}

	a::Function	
	g::Function
	f::Function
	H::Function
	
	delta::Float64 # regularization REMOVE!!!!

	function class_quadratic_program()
		this = new();
		
		this.f = function(vars::class_variables)
			x_scaled = vars.x()/vars.tau();
			return (this.q'*x_scaled + 0.5 * x_scaled'*this.H(vars)*x_scaled)[1]
		end

		#this.a = function(vars::class_variables)
		#	return this.A*vars.x() / vars.tau() - this.b;
		#end

		#this.g = function(vars::class_variables)
		#	objective_gradient = this.c + this.H(vars) * vars.x()/vars.tau();
		#	return objective_gradient 
		#end
	
		this.H = function(vars::class_variables)
			return this._H;
		end

		return this;	
	end

	function class_quadratic_program(A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1},c::Array{Float64,1})
		this = class_quadratic_program();

		(this.m, this.n) = size(A)
		this.A = A;
		this._H = spzeros(this.n,this.n)
		this.q = c;
		this.b = b;

		return this;
	end

	function class_quadratic_program(A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1},c::Array{Float64,1},Q::SparseMatrixCSC{Float64,Int64})
		this = class_quadratic_program();

		(this.m, this.n) = size(A)
		this.A = A;
		this._H = Q;
		this.q = c;
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

	ipopt_start::Function
	check_positive::Function
	check_direction::Function
	
	x_scaled::Function
	x::Function
	y::Function
	s::Function
	tau::Function
	kappa::Function

	n::Function
	m::Function

	function class_variables(n::Int64, m::Int64)	
		this = new();
		# construct indicies
		this.x_ind = 1:n;
		this.y_ind = (n+1):(n+m);
		this.s_ind = (n+m+1):(2*n + m);
		this.tau_ind = 2*n + m + 1;
		this.kappa_ind = 2*n + m + 2;		

		# intialize
		intial_val = 1.0;
		this.v = zeros(2*(n+1)+m);
		this.v[this.x_ind] = intial_val;
		this.v[this.s_ind] = intial_val;
		this.v[this.tau_ind] = intial_val;
		this.v[this.kappa_ind] = intial_val;

		this.x_scaled = function()
			return this.x()/this.tau();
		end

		this.x = function()
			try
				return this.v[this.x_ind]
			catch e
				println(typeof(this.v))
				println("x error")
				throw(e)
			end
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

		this.n = function()
			return length(this.x_ind)
		end

		this.m = function()
			return length(this.y_ind)
		end
		
		this.check_direction = function()
			# check for infinity
			for x_i in this.x()
				@assert(x_i != -Inf)
				@assert(x_i != Inf)
			end

			for s_i in this.s()
				@assert(s_i != -Inf)
				@assert(s_i != Inf)
			end

			@assert(this.tau() != -Inf)
			@assert(this.tau() != Inf)
			@assert(this.kappa() != -Inf)
			@assert(this.kappa() != Inf)

			# check for NaNs
			for x_i in this.x()
				@assert(x_i != NaN)
			end

			for s_i in this.s()
				@assert(s_i != NaN)
			end

			@assert(this.tau() != NaN)
			@assert(this.kappa() != NaN)
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

		this.ipopt_start = function()
			this.v[this.tau_ind] = 1.0;
			this.v[this.kappa_ind] = 1e-20;
		end

		return this;
	end
end

function validate_dimensions(qp::class_quadratic_program,vars::class_variables)
	try
		@assert(size(qp.H(vars)) == (qp.n,qp.n))
		@assert(size(qp.A) == (qp.m,qp.n))
		@assert(size(qp.b) == (qp.m,))
		@assert(size(qp.q) == (qp.n,))
		@assert(length(vars.x()) == qp.n)
		@assert(length(vars.y()) == qp.m)
		@assert(length(vars.v) == 2*(qp.n + 1) + qp.m) 
 	catch e
		println("ERROR validate_dimensions")
		throw(e)
	end
end






