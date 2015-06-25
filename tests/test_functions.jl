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
	derivative = local_approx.calculate_residual_norm_derivative(nlp, vars, settings);
	finite_difference_derivative = local_approx.finite_difference_merit_function_derivative(nlp, vars, settings);
	
	h = 1e-7
	tol = 1e-2
	
	check_derivative_r_norm_derivative_of_point(h, tol, nlp, local_approx, vars, settings, orginal_r_norm, vars.x, derivative.dx, "dx")
	check_derivative_r_norm_derivative_of_point(h, tol, nlp, local_approx, vars, settings, orginal_r_norm, vars.x_bar, derivative.dx_bar, "dx_bar")
	check_derivative_r_norm_derivative_of_point(h, tol, nlp, local_approx, vars, settings, orginal_r_norm, vars.y, derivative.dy, "dy")
	check_derivative_r_norm_derivative_of_point(h, tol, nlp, local_approx, vars, settings, orginal_r_norm, vars.y_bar, derivative.dy_bar, "dy_bar")
	check_derivative_r_norm_derivative_of_point(h, tol, nlp, local_approx, vars, settings, orginal_r_norm, vars.s, derivative.ds, "ds")
	check_derivative_r_norm_derivative_of_point(h, tol, nlp, local_approx, vars, settings, orginal_r_norm, vars.z, derivative.dz, "dz")
	
	vars.tau = vars.tau + h
	local_approx.update_approximation(nlp, vars, settings);
	new_r_norm = local_approx.state.r_norm
	derivative = local_approx.calculate_residual_norm_derivative(nlp, vars, settings);
	if (new_r_norm - orginal_r_norm)/h -  derivative.dtau > tol
		println("dtau")
		println("ERROR, finite difference does not agree with derivative")
	end
	vars.tau = vars.tau - h;
	
	
	vars.kappa = vars.kappa + h
	local_approx.update_approximation(nlp, vars, settings);
	new_r_norm = local_approx.state.r_norm
	derivative = local_approx.calculate_residual_norm_derivative(nlp, vars, settings);
	if (new_r_norm - orginal_r_norm)/h -  derivative.dkappa > tol
		println("dkappa")
		println("ERROR, finite difference does not agree with derivative")
	end
	vars.kappa = vars.kappa - h
end


function finite_difference_r_norm(h, local_approx, nlp, vars, settings, orginal_value, point)
	result = zeros(length(point))
	
	for i in 1:length(point)
		point[i] = point[i] + h
		local_approx.update_approximation(nlp, vars, settings);
		new_value = local_approx.state.r_norm
		derivative = local_approx.calculate_residual_norm_derivative(nlp, vars, settings);
		
		result[i] = (new_value - orginal_value)/h
		
		point[i] = point[i] - h
	end
		
	return result
end


function check_derivative_r_norm_derivative_of_point(h, tol, nlp, local_approx, vars, settings, orginal_value, point, point_derivative, txt)
	result = finite_difference_r_norm(h, tol, nlp, local_approx, vars, settings, orginal_value, point)
	
	if norm(result - point_derivative,2) > tol
		println(txt)
		println("ERROR, finite difference does not agree with derivative")
	end
end

function check_merit_function_derivative(vars::class_variables, nlp::class_non_linear_program,  settings::class_settings)
	local_approx = class_local_approximation();
	orginal_merit_value = local_approx.update_approximation(nlp, vars, settings);
	derivative = local_approx.calculate_merit_function_derivative(nlp, vars, settings);
	
	derivative_vector = derivative.get_point();
	variable_vector = 
	
	h = 1e-6
	
	for i = 1:length(vars.point)
		vars.point[i] = vars.point[i] + h;
		new_merit_value = local_approx.update_approximation(nlp, vars, settings);
		
		print_results((new_merit_value - orginal_merit_value)/h, derivative.point[i], "dtau")
		vars.point[i] = vars.point[i] - h;
	end
end

settings = standard_settings()

A = -0.789*ones(1,1);
A_bar = spzeros(0,1);
	
c = [2.578];

b = -[3.14];
b_bar = EMPTY_ARRAY;

lp = class_non_linear_program();
lp.set_linear_objective(c);
lp.set_linear_constraints(sparse(A),b,sparse(A_bar),b_bar,0);


vars = class_variables(lp);
vars.randomize();

local_approx = class_local_approximation();
local_approx.update_approximation(lp, vars, settings);
orginal_r_norm = local_approx.state.r_norm
derivative = local_approx.calculate_residual_norm_derivative(lp, vars, settings);
finite_difference_derivative = local_approx.finite_difference_merit_function_derivative(lp, vars, settings);

#vars, status = ip_algorithm(lp, settings, vars, false);

if false
	check_r_norm_derivative(vars, lp, settings);

	println("===================================")

	A_bar = spzeros(0,2); 
		
	A = speye(2);

	c = [1.0, 0.0];

	b = [0.5, 0.5];
	b_bar = EMPTY_ARRAY;

	lp = class_non_linear_program();
	lp.set_linear_objective(c);

	lp.set_linear_constraints(A,b,A_bar,b_bar,0);

	vars = class_variables(lp);
	vars.randomize();

	#vars, status = ip_algorithm(lp, settings, vars, false);

	check_r_norm_derivative(vars, lp, settings);


	println("===================================")

	A = spzeros(0,2)
	A_bar = sparse(ones(1,2)); 

	c = [1.0, 0.0];

	b = EMPTY_ARRAY;
	b_bar = [0.5];

	lp = class_non_linear_program();
	lp.set_linear_objective(c);

	lp.set_linear_constraints(A,b,A_bar,b_bar,2);

	vars = class_variables(lp);
	vars.randomize();

	#vars, status = ip_algorithm(lp, settings, vars, false);

	check_r_norm_derivative(vars, lp, settings);

	#println("===================================")
	#check_merit_function_derivative(vars, lp, settings);

end