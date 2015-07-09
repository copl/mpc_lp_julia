include("../../tests/testing_tools.jl")
using MAT
using SCS

function main()
	settings = standard_settings()
	
	#run_net_lib_problem("problems/blend.mat", settings)
	
	dir = "../../tests/Problems";
	netlib_problems = readdir(dir)
	
	successful_problems = 0
	ipopt_successful_problems = 0
		
	problem_name = "25FV47.mat" #"BRANDY.mat"
	A, b, c = get_netlib_problem(dir * "/" * problem_name)
	
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(spzeros(0,length(c)),EMPTY_ARRAY,A,b,length(c));
	
	vars = class_variables(lp);
	vars, status, original_iter = ip_algorithm(lp, settings, vars);
	

	# calls SCS solver to approximately solve (eps=1e-1) linear programs
	success, x_star, y_star = solve_with_JuMP(A, b, c, SCSSolver(eps=1e-1))
	

	x_star = max(0,x_star) # force solution to be positive
	
	# finds the optimal solution
	

	# new perturbed problem
	lambda=0.99
	vars = class_variables(lp)
    	#x_star = (vars.x/vars.tau) 
    	vars.x=lambda*x_star + (1-lambda)*ones(length(vars.x))
    	vars.s=(1-lambda)./vars.x
    	vars.y_bar=zeros(length(vars.y_bar))
    	vars.tau=1.0
    	vars.kappa=(1-lambda)
	# need to update vars
	#
	
	#K=10
	#srand(1)
	#c = c + K * rand(length(c))/(length(c))
	#lp.set_linear_objective(c);
	
	vars, status, new_iter = ip_algorithm(lp, settings, vars);
	
	println("Orginal: ",original_iter, " New: ", new_iter)
end

main();
