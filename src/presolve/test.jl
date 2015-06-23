using MAT
include("../testing_tools.jl")
println("loaded libraries")

A,b,c = get_netlib_problem("../Problems/blend.mat")
# min c'*x 
# A x = b
# x >= 0


