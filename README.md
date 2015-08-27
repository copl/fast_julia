## Julia implementation of the homogeneous algorithm
* Requires MUMPS.jl https://github.com/JuliaSparse/MUMPS.jl, JuMP Pkg.add('JuMP') and IPOPT Pkg.add('Ipopt')
* To test run the file `src/homogeneous_algorithm/run_tests.jl'. This file runs a netlib problem using different linear equation solvers. It breaks down the time each part of the algorithm takes. It also compares the perforance to Ipopt.
* In this file, it shows how you can modify the settings to change between linear system solvers.
* To add or modify existing linear system solvers see `src/homogeneous_algorithm/linear_system_solver.jl'.
