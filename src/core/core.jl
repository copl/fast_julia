abstract internal_AbstractNLPEvaluator;

include("class_variables.jl")
include("class_nlp_cache.jl")
include("class_quadratic_program.jl")
include("class_settings.jl")


function validate_dimensions( nlp_eval::internal_AbstractNLPEvaluator, vars::class_variables )
	try
    @assert(n(vars) == n(nlp_eval))
    @assert(m(vars) == m(nlp_eval))
 	catch e
	  println("ERROR validate_dimensions")
		throw(e)
	end
end

function validate_dimensions( nlp_val::class_nlp_cache )
    #
    # TO DO !!!!
    #

    try
        num_vars = n(nlp_eval)
        num_constraints = m(nlp_eval)

        @assert( size( internal_eval_hesslag_prod(nlp_eval, x_val, y_val) ) == (num_vars,  num_vars) )
        @assert( size( internal_eval_jac_a(  nlp_eval, x_var) ) == ( num_constraints,  num_vars) )
        @assert( size( internal_eval_a(  nlp_eval, x_var) ) == ( num_constraints ,))
        @assert( size( internal_eval_gradlag( nlp_eval, x_var, y_var) ) == ( num_vars,))
        @assert ( size( internal_eval_gradc( nlp_eval, x_val )) == ( num_vars,) )
        @assert ( size( internal_eval_b( nlp_eval, y_val )) == ( num_constraints ,) )
    catch e
        println("ERROR validate_dimensions")
        throw(e)
	  end
end
