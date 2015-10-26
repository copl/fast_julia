abstract internal_AbstractNLPEvaluator;

include("class_variables.jl")
include("class_quadratic_program.jl")

type nlp_cache
    val_f::Float64
    val_g::Array{Float64,1}
    val_jac_g::SparseMatrixCSC{Float64,Int64}
    val_hesslag_prod::SparseMatrixCSC{Float64,Int64}
    val_gradlag::Array{Float64,1}
end

function validate_dimensions( nlp_eval::internal_AbstractNLPEvaluator, vars::class_variables )
	try
    x_val = x(vars);
    y_val = x(vars);

    println( size( internal_eval_hesslag_prod(nlp_eval, x_val, y_val) ) )
    println( nlp_eval._n )

		@assert( size( internal_eval_hesslag_prod(nlp_eval, x_val, y_val) ) == ( nlp_eval._n,  nlp_eval._n) )
		@assert( size( internal_eval_jac_g(  nlp_eval, x_var) ) == ( nlp_eval._m,  nlp_eval._n) )
		@assert( size( internal_eval_g(  nlp_eval, x_var) ) == ( nlp_eval._m,))
		@assert( size( internal_eval_gradlag( nlp_eval, x_var, y_var) ) == ( nlp_eval._n,))
 	catch e
	  println("ERROR validate_dimensions")
		throw(e)
	end
end

