
println("loading libraries")


#using JuMP
#using Gurobi
#using Ipopt
#using KNITRO
#using SCS
using MAT

println("external libraries loaded")

include("ip_algorithm.jl")
include("testing_tools.jl")

println("internal libraries loaded")

function main()
	settings = class_settings();
	
	settings.max_iter = 70;  # Total number of iterarions
	settings.max_iter_line_search = 70;
	
	settings.primal_feas_tol = 1e-8
	settings.dual_feas_tol = 1e-8
	settings.duality_gap_tol = 1e-10
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
	iter_list = zeros(1,0);
	for problem_name = netlib_problems		
		try
			println("Solving ", problem_name)
			success, iter = run_net_lib_problem(dir * "/" * problem_name, settings)

			if success
				successful_problems += 1
				iter_list = [iter_list [iter]];
			end
		catch e
			println(e)
		end
	end
	
	println("Solved ", successful_problems, " out of ", length(netlib_problems))
	println("Average iterations ", mean(iter_list))
end

function run_net_lib_problem(problem_name,settings)
	file = matopen(problem_name)
	A_bar = sparse(read(file, "A"))
	n,k = size(A_bar)	
	b_bar = read(file,"b")
	b_bar = b_bar[:]
	c = read(file,"c")
	c = c[:];
	
	#println("A matrix is size " * string((n,k)) * " and has " * string(nnz(A_bar))  * " non-zeros")
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(spzeros(0,k),spzeros(0,1)*[1.0],A_bar,b_bar,k);
	
	vars = class_variables(lp);
	
	vars, status, iter = ip_algorithm(lp, settings, vars, false);
	
	return is_problem_successful(problem_name, status, 1, None), iter
end

main();