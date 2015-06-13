
println("loading libraries")


using JuMP
#using Gurobi
using Ipopt
#using KNITRO
#using SCS
using MAT
using Mosek


println("external libraries loaded")

include("ip_algorithm.jl")
include("testing_tools.jl")

println("internal libraries loaded")

function main()
	settings = class_settings();
	
	settings.max_iter = 70;  # Total number of iterarions
	settings.max_iter_line_search = 70;
	
	settings.primal_feas_tol = 1e-10
	settings.dual_feas_tol = 1e-10
	settings.duality_gap_tol = 1e-12
	settings.primal_infeas_tol = 1e-10
	settings.dual_infeas_tol = 1e-10
	
	settings.beta1 = 10.0^(-4)
	settings.beta2 = 10.0^(-8)
	settings.beta3 = 10.0^(-4)
	settings.beta4 = 0.9
	settings.beta5 = 0.999
	settings.beta6 = 0.5
	
	#run_net_lib_problem("problems/blend.mat", settings)
	
	dir = "small_problems";
	netlib_problems = readdir(dir)
	
	successful_problems = 0
	ipopt_successful_problems = 0
	
	iter_list = zeros(1,0);
	
	for problem_name = netlib_problems	
		A, b, c = get_netlib_problem(dir * "/" * problem_name)
		
		println("Solving ", problem_name)
		try
			#ipopt_success = solve_with_JuMP(A, b, c, MosekSolver())
			#if ipopt_success
			#	ipopt_successful_problems += 1
			#end
		catch e
			println(e)
		end
		
		try
			status, iter = solve_net_lib_problem(A,b,c,settings)
			
			if is_problem_successful(problem_name, status, 1, None)
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

function get_netlib_problem(file_name)
	file = matopen(file_name)
	A= sparse(read(file, "A"))
	b = read(file,"b")
	b = b[:]
	c = read(file,"c")
	c = c[:];
	
	return A,b,c
end

function solve_with_JuMP(A, b, c, solver=IpoptSolver(max_iter=300)) #Output=0
	model = Model(solver=solver)
	
	n, k = size(A)
	
	@defVar(model, x[1:k] >= 0)
	@setObjective( model, Min, sum{ c[j]*x[j] , j=1:k} )

	for i = 1:n
		@addConstraint(model, sum{A[i,j]*x[j], j=1:k} == b[i])
	end
		
	status = solve(model)
	
	#println("Objective value: ", getObjectiveValue(model))

	return status == :Optimal
end

function solve_net_lib_problem(A,b,c,settings)
	n,k = size(A)
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(spzeros(0,k),spzeros(0,1)*[1.0],A,b,k);
	
	vars = class_variables(lp);
	
	vars, status, iter = ip_algorithm(lp, settings, vars, false);
	
	return status, iter
end

main();