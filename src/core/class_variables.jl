
###############
# data
################


type class_variables
	# private variables
  _v::Array{Float64,1}
	_x_ind::UnitRange{Int64}
	_y_ind::UnitRange{Int64}
	_s_ind::UnitRange{Int64}
	_tau_ind::Int64
	_kappa_ind::Int64

	function class_variables(n::Int64, m::Int64)
		this = new();
		# construct indicies
		this._x_ind = 1:n;
		this._y_ind = (n+1):(n+m);
		this._s_ind = (n+m+1):(2*n + m);
		this._tau_ind = 2*n + m + 1;
		this._kappa_ind = 2*n + m + 2;

		# intialize
		intial_val = 1.0;
		this._v = zeros(2*(n+1)+m);
		this._v[this._x_ind] = intial_val;
		this._v[this._s_ind] = intial_val;
		this._v[this._tau_ind] = intial_val;
		this._v[this._kappa_ind] = intial_val;

		return this;
	end
end

##################################
# methods
##################################

function x(vars::class_variables)
    return vars._v[vars._x_ind]
end

function y(vars::class_variables)
    return vars._v[vars._y_ind]
end

function s(vars::class_variables)
    return vars._v[vars._s_ind]
end

function tau(vars::class_variables)
    return vars._v[vars._tau_ind]
end

function kappa(vars::class_variables)
    return vars._v[vars._kappa_ind]
end

function x_scaled(vars::class_variables)
     return x(vars)/tau(vars)
end

function n(vars::class_variables)
    return length(this._x_ind)
end

function m(vars::class_variables)
    return length(this._y_ind)
end

# point = point + alpha * direction
function move!(vars::class_variables, alpha::Float64, direction::class_variables)
    vars._v = vars._v + alpha * direction._v;
    check_positive(vars);
end

function ipopt_start()
	 this.v[this.tau_ind] = 1.0;
	 this.v[this.kappa_ind] = 1e-20;
end

# check for errors

# check x > 0
function check_positive(vars::class_variables)
      try
        for x_i in x(vars)
           @assert(x_i > 0)
        end

        for s_i in s(vars)
           @assert(s_i > 0)
        end

        @assert(tau(vars) > 0)
        @assert(kappa(vars) > 0)

    catch e
        println("ERROR in check_positive")
        throw(e)
    end

     check_for_wrong_vals(vars)
end

function check_for_wrong_vals(vars::class_variables)
      try
          # check for infinity
          for x_i in x(vars)
            @assert(x_i != -Inf)
            @assert(x_i != Inf)
          end

          for s_i in s(vars)
            @assert(s_i != -Inf)
            @assert(s_i != Inf)
          end

          @assert(tau(vars) != -Inf)
          @assert(tau(vars) != Inf)
          @assert(kappa(vars) != -Inf)
          @assert(kappa(vars) != Inf)

          # check for NaNs
          for x_i in x(vars)
            @assert(isnan(x_i))
          end

          for s_i in s(vars)
            @assert(isnan(s_i))
          end

          @assert(isnan(tau(vars)))
          @assert(isnan(kappa(vars)))
      catch e
          println("ERROR in check_for_wrong_vals")
          throw(e)
      end
end

