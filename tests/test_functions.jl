# this file contains a series of simple ip test functions
#
# TO DO:
# - write function that records failure/sucess of tests
# - test inner workings of IP method
#		- merit function line search
# 		- 
include("testing_tools.jl")

function print_results(finite_difference, code_computed, label)
	println("Finite difference ", label, ":", finite_difference)
	println("Code computed ", label, ":", code_computed)
end


function check_r_norm_derivative(vars::class_variables, nlp::class_non_linear_program,  settings::class_settings)
	# pick a random point with x > 0 and a problem
	# compute the local approximation
	# add x + eh to compute the derivative of r_norm 
	# compare that with the actual derivative
	local_approx = class_local_approximation();
	local_approx.update_approximation(nlp, vars, settings);
	orginal_r_norm = local_approx.state.r_norm
	
	h = 1e-6
	
	vars.tau = vars.tau + h
	local_approx.update_approximation(nlp, vars, settings);
	new_r_norm = local_approx.state.r_norm
	derivative = local_approx.calculate_residual_norm_derivative(nlp, vars, settings);
	print_results((new_r_norm - orginal_r_norm)/h, derivative.dtau, "dtau")
	vars.tau = vars.tau - h;
	
	vars.kappa = vars.kappa + h
	local_approx.update_approximation(nlp, vars, settings);
	new_r_norm = local_approx.state.r_norm
	derivative = local_approx.calculate_residual_norm_derivative(nlp, vars, settings);
	print_results((new_r_norm - orginal_r_norm)/h, derivative.dkappa, "dkappa")
	vars.kappa = vars.kappa - h
	
	vars.x = vars.x + h
	local_approx.update_approximation(nlp, vars, settings);
	new_r_norm = local_approx.state.r_norm
	derivative = local_approx.calculate_residual_norm_derivative(nlp, vars, settings);
	print_results((new_r_norm - orginal_r_norm)/h, derivative.dx, "dx")
	vars.x = vars.x - h
	
	vars.y = vars.y + h
	local_approx.update_approximation(nlp, vars, settings);
	new_r_norm = local_approx.state.r_norm
	derivative = local_approx.calculate_residual_norm_derivative(nlp, vars, settings);
	print_results((new_r_norm - orginal_r_norm)/h, derivative.dy, "dy")
	vars.y = vars.y - h
	
	vars.s = vars.s + h
	local_approx.update_approximation(nlp, vars, settings);
	new_r_norm = local_approx.state.r_norm
	derivative = local_approx.calculate_residual_norm_derivative(nlp, vars, settings);
	print_results((new_r_norm - orginal_r_norm)/h, derivative.ds, "ds")
	vars.s = vars.s - h
	
	
end

settings = standard_settings()

A = -0.789*ones(1,1);
A_bar = spzeros(0,1);
	
c = [2.578];

b = -[3.14];
b_bar = EMPTY_ARRAY;

lp = class_non_linear_program();
lp.set_linear_objective(c);
lp.set_linear_constraints(sparse(A),b,sparse(A_bar),b_bar,1);


vars = class_variables(lp);
vars.randomize();

#vars, status = ip_algorithm(lp, settings, vars, false);

check_r_norm_derivative(vars, lp, settings);

