
println("loading libraries")


using JuMP
#using Gurobi
using Ipopt
#using KNITRO
#using SCS
using MAT
using Mosek


println("external libraries loaded")

include("testing_tools.jl")

println("internal libraries loaded")

function main()
	settings = standard_settings()
	
	#run_net_lib_problem("problems/blend.mat", settings)
	
	dir = "small_problems";
	netlib_problems = readdir(dir)
	
	successful_problems = 0
	ipopt_successful_problems = 0
	
	iter_list = zeros(1,0);
	
	for problem_name = netlib_problems	
		A, b, c = get_netlib_problem(dir * "/" * problem_name)
		
		ipopt_success = 0;
		println("Solving ", problem_name)
		try
			ipopt_success = solve_with_JuMP(A, b, c, MosekSolver(LOG=0))
		catch e
			println(e)
		end
		
		try
			status, iter = solve_net_lib_problem(A,b,c,settings)
			
			if is_problem_successful(problem_name, status, ipopt_success, None)
				successful_problems += 1
				iter_list = [iter_list iter];
			end
		catch e
			println(e)
		end
	end
	
	println("Solved ", successful_problems, " out of ", length(netlib_problems))
	println("Average iterations ", mean(iter_list))
	#println("IPOPT solved ", ipopt_successful_problems, " out of ", length(netlib_problems))
end

main();