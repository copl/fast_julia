type class_nlp_cache
    # min f(x)
    # a(x) = 0
    # x >= 0
    val_c::Float64
    val_a::Array{Float64,1}
    val_jac_a::SparseMatrixCSC{Float64,Int64}
    val_hesslag_prod::SparseMatrixCSC{Float64,Int64}
    val_gradlag::Array{Float64,1}
    val_gradc::Array{Float64,1}
    val_b::Array{Float64,1}

   function class_nlp_cache()
      return new();
  end
end


function update_nlp_cache!(nlp_vals::class_nlp_cache, nlp_eval::internal_AbstractNLPEvaluator, vars::class_variables)
    try
        val_x_scaled = x_scaled(vars);
        val_y_scaled = y_scaled(vars);

        nlp_vals.val_c = internal_eval_c( nlp_eval, val_x_scaled )
        nlp_vals.val_a = internal_eval_a( nlp_eval, val_x_scaled )
        nlp_vals.val_jac_a = internal_eval_jac_a( nlp_eval, val_x_scaled )
        nlp_vals.val_hesslag_prod = internal_eval_hesslag_prod( nlp_eval, val_x_scaled, val_y_scaled )
        nlp_vals.val_gradlag = internal_eval_gradlag( nlp_eval, val_x_scaled, val_y_scaled )

        nlp_vals.val_gradc = internal_eval_gradc( nlp_eval, val_x_scaled )
        nlp_vals.val_b = internal_eval_b( nlp_eval, val_x_scaled )
    catch e
        println("class_nlp_cache/update_nlp_cache!")
        throw(e)
    end
end
