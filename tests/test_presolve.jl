using MAT
include("testing_tools.jl")
println("loaded libraries")

A,b,c = get_netlib_problem("small_problems/BORE3D.mat")
# min c'*x
# A x = b
# x >= 0



