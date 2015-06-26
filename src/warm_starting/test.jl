include("../../tests/testing_tools.jl")
using MAT

function main()
	settings = standard_settings()
	
	#run_net_lib_problem("problems/blend.mat", settings)
	
	dir = "../../tests/small_problems";
	netlib_problems = readdir(dir)
	
	successful_problems = 0
	ipopt_successful_problems = 0
		
	problem_name = "blend.mat"
	A, b, c = get_netlib_problem(dir * "/" * problem_name)
	
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(spzeros(0,length(c)),EMPTY_ARRAY,A,b,length(c));
	
	vars = class_variables(lp);
	
	# vars.x
	# vars.y_bar
	# vars.s
	# vars.tau
	# vars.kappa
	
	vars, status, number_of_iterations = ip_algorithm(lp, settings, vars);
	
	# finds the optimal solution
	
	# new perturbed problem
	
	# HERE
	# need to update vars
	#
	
	c = c + rand(length(c))/sqrt(length(c))
	lp.set_linear_objective(c);
	
	ip_algorithm(lp, settings, vars);
	
end

main();