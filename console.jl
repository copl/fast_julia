include("benchmark_ncqp_netlib.jl")

using my_advanced_timer

test_timer = advanced_timer();

cd("Documents")
cd("Programming_Projects")
cd("fast_julia")
cd("tests")
print(pwd())
test_problem("BLEND")
