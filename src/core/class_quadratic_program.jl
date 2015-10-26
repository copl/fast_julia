# TO DO
# link with math programming base!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type class_quadratic_program <: internal_AbstractNLPEvaluator
  # c'*x + 0.5 * x' * Q * x
  # A*x = b
  # x >= 0

	_n::Int64 # number of variables
	_m::Int64 # number of constraints
	_A::SparseMatrixCSC{Float64,Int64}
	_b::Array{Float64,1}
  _c::Array{Float64,1}
  _Q::SparseMatrixCSC{Float64,Int64}


	function class_quadratic_program(A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1},c::Array{Float64,1})
    (m, n) = size(A)
		return new(n,m,A,b,c,spzeros(n,n));
	end

	function class_quadratic_program(A::SparseMatrixCSC{Float64,Int64},b::Array{Float64,1},c::Array{Float64,1},Q::SparseMatrixCSC{Float64,Int64})
		(m, n) = size(A)
		return new(n,m,A,b,c,Q);
	end
end

# evaluate objective
function internal_eval_f(qp::class_quadratic_program, x::Array{Float64,1})
    return dot(qp._c, x) + 0.5 * dot(x, qp._Q * x);
end

# evalutate constraints
function internal_eval_g(qp::class_quadratic_program, x::Array{Float64,1})
    return qp._A * x - qp._b;
end

# evaluate gradient of constraints
function internal_eval_jac_g(qp::class_quadratic_program, x::Array{Float64,1}) # J
    return qp._A;
end

# hessian of lagrangian
function internal_eval_hesslag_prod(qp::class_quadratic_program, x::Array{Float64,1}, y::Array{Float64,1})
    return qp._Q
end

# gradient of lagrangian
function internal_eval_gradlag(qp::class_quadratic_program, x::Array{Float64,1}, y::Array{Float64,1}) # \nabla L
    return qp._c + qp._Q * x + qp._A' * y;
end


###############
# tests
