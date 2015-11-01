include("benchmark_ncqp_netlib.jl")

using my_advanced_timer

test_timer = advanced_timer();



##############################################################




function test(R::SparseMatrixCSC{Float64,Int64},v::Array{Float64,1})
  println("cholfac")
  @time F = cholfact(R + 1e-7 * speye(size(R,1)));


  #println(typeof(F))
  #F::CholmodFactor{Float64,Int64}
  println("solve")

  @time begin
    for i = 1:50
        F \ v
    end
  end
  #println("cholfacdense")
  #@time F_dense = cholfact(eye(100))

   #println("cholsolvedense")
  #@time F_dense \ u

end

function test2(A::SparseMatrixCSC{Float64,Int64},v::Array{Float64,1})
  println("lufac")
  @time F = lufact([speye(size(A,2))]);


  #println(typeof(F))
  #F::CholmodFactor{Float64,Int64}
  println("solve")

  @time begin
    for i = 1:50
        F \ v # compare with mumps
    end
  end
  #println("cholfacdense")
  #@time F_dense = cholfact(eye(100))

   #println("cholsolvedense")
  #@time F_dense \ u

end

A,b,c = get_netlib_problem("BLEND.mat")
R = A*A';
v = ones(size(R,1));
test(R,v)

using MUMPS

