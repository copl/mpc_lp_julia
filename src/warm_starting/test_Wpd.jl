include("../../tests/testing_tools.jl")
using MAT
using SCS

function main()
	settings = standard_settings()
	
	#run_net_lib_problem("problems/blend.mat", settings)
	
	dir = "../../tests/small_problems";
	netlib_problems = readdir(dir)
	
	successful_problems = 0
	ipopt_successful_problems = 0
		
	problem_name = "BLEND.mat"
	A, b, c = get_netlib_problem(dir * "/" * problem_name)
	
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(spzeros(0,length(c)),EMPTY_ARRAY,A,b,length(c));
	
	vars = class_variables(lp);

	vars, status, number_of_iterations = ip_algorithm(lp, settings, vars);
	
	
	# calls SCS solver to approximately solve (eps=1e-1) linear programs
	success, x_star, y_star = solve_with_JuMP(A, b, c, SCSSolver(eps=1e-2))
	
	s_star = c - A'*y_star
	s_star = max(0,s_star) # force solution to be positive
	
	x_star = max(0,x_star) # force solution to be positive

	
	# new perturbed problem
	#x_star = (vars.x/vars.tau)
	#s_star = (vars.s/vars.tau)
	#y_star = (vars.y_bar/vars.tau)
	(num_rows,num_columns)=size(A)
	lambda=0.9
	
    	vars.x=lambda*x_star+(1-lambda)*ones(length(vars.x))
    	vars.s=lambda*s_star+(1-lambda)*ones(length(vars.s))
    	vars.y_bar=lambda*y_star
   	vars.tau=1
    	vars.kappa=1-lambda  #(vars.x'*vars.s)[1]/num_rows #
    
    
	# need to update vars
	ip_algorithm(lp, settings, vars);
	
end

main();
