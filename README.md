## Julia implementation of the homogeneous algorithm
* To test run the file `src/homogeneous_algorithm/run_tests.jl'. This file runs a netlib problem using different linear equation solvers. It breaks down the time each part of the algorithm takes. It also compares the perforance to Ipopt.
* Requires MUMPS.jl https://github.com/JuliaSparse/MUMPS.jl, JuMP Pkg.add('JuMP') and IPOPT Pkg.add('Ipopt')
