# this file contains a series of small and simple tests optimization problems
# 

println("loading libraries")


function debug_message(message)
	#println(message)
end
function debug_message_high_level(message)
	#println(message)
end

function debug_println(message,shouldPrint)
	if shouldPrint
		println(message)
	end
end
#using JuMP
#using Gurobi
#using Ipopt
#using KNITRO
#using SCS
using MAT

println("external libraries loaded")

include("ip_algorithm.jl")

println("internal libraries loaded")

function main()
	settings = class_settings();
	
	settings.max_iter = 50;  # Total number of iterarions
	settings.max_iter_line_search = 50;
	
	settings.primal_feas_tol = 1e-8
	settings.dual_feas_tol = 1e-8
	settings.duality_gap_tol = 1e-10
	settings.primal_infeas_tol = 1e-8
	settings.dual_infeas_tol = 1e-8
	
	settings.beta1 = 10.0^(-4)
	settings.beta2 = 10.0^(-8)
	settings.beta3 = 10.0^(-4)
	settings.beta4 = 0.9
	settings.beta5 = 0.999
	settings.beta6 = 0.5

	if true
	
		# test problems
		#trival_problem1(settings)
		#trival_problem2(settings)
		#trival_problem3(settings)
		
		#trival_infeasible1(settings)
		#trival_infeasible2(settings)
		#trival_unbounded1(settings)
		#trival_unbounded2(settings)
		
		#lp_feasible_and_bounded_1(settings)
		#lp_feasible_and_bounded_2(settings)
		#lp_feasible_and_bounded_3(settings)
		#lp_feasible_and_bounded_4(settings)
		
		#qp_feasible_and_bounded_1(settings)
		qp_infeasible_1(settings)
		#qp_unbounded_1(settings)
		lp_blend_problem_1(settings)
		#lp_blend_problem_2(settings)
		#lp_blend_problem_infeasible_1(settings)
		#lp_blend_problem_unbounded_1(settings)
		#qp_blend_problem_1(settings)
	end
	
	if true
		L2_problem(settings, 1.0)
		L2_problem(settings, 0.0)
		L2_problem(settings, -1.0)
		L2_problem(settings, -10.0)
		L2_problem_inequality(settings,-1.0)
	end
end

############################
# linear programs
############################

# A = [[1,1],[1,1]]

EMPTY_ARRAY = spzeros(0,1)*[0.0]

function trival_problem1(settings)
	println("trival_problem1")
	println("min x")
	println("x >= 0")
	
	A_bar = spzeros(0,1); 
	
	A = spzeros(0,1);
	
	c = [1.0];
	
	b = EMPTY_ARRAY;
	b_bar = EMPTY_ARRAY;
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(sparse(A),b,sparse(A_bar),b_bar,1);
	
	vars = class_variables(lp);
	
	if true
		vars.x = [2.0];
		vars.kappa = 5.0;
		vars.s = [3.0];
		vars.tau = 4.0;
		vars.x_scaled = [0.5];
		vars.y_scaled = [];
	end
	
	ip_algorithm(lp, settings, vars);
	
	# ans
	# x = [4.0, 0.0, 0.0, 0.0]
	
end
function trival_problem2(settings)
	println("trival_problem2")
	println("min x[1]")
	println("x[1], x[2] >= 0")

	A_bar = spzeros(0,2); 
	
	A = spzeros(0,2);
	
	c = [1, 0.0];
	
	b = EMPTY_ARRAY;
	b_bar = EMPTY_ARRAY;
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(sparse(A),b,sparse(A_bar),b_bar,2);
	
	vars = class_variables(lp);
	
	ip_algorithm(lp, settings, vars);
	
	# ans
	# x = [4.0, 0.0, 0.0, 0.0]
	
end

function trival_problem3(settings)
	println("trival_problem3")
	println("min x[1]")
	println("eye(2)*x >= 0")
	
	A_bar = spzeros(0,2); 
	
	A = speye(2);
	
	c = [1.0, 0.0];
	
	b = [0.0, 0.0];
	b_bar = EMPTY_ARRAY;
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,0);
	
	vars = class_variables(lp);
	
	println("trival_problem3")
	println("min",c,"*x")
	println(full(A),"*x >=",b)
	
	ip_algorithm(lp, settings, vars);
	
end

function trival_infeasible1(settings)
	println("trival_infeasible1")
	println("min x")
	println("x == -1")
	println("x >= 0")
	
	A_bar = speye(1); 
	
	A = spzeros(0,1);
	
	c = [1.0];
	
	b_bar = [-1.0];
	b = A*[0.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,1);
	
	vars = class_variables(lp);
	
	ip_algorithm(lp, settings, vars);
	
end

function trival_infeasible2(settings)
	println("trival_infeasible2")
	println("x <= -1, x >= 0")
	
	A_bar = spzeros(0,1); 
	
	A = -speye(1);
	
	c = [1.0];
	
	b = [1.0];
	b_bar = A_bar*[0.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,1);
	
	vars = class_variables(lp);
	
	ip_algorithm(lp, settings, vars);
	
end

function trival_unbounded1(settings)
	println("trival_unbounded1")
	println("min -x")
	println("x >= 0")
	
	A_bar = spzeros(0,1); 
	
	A = spzeros(0,1);
	
	c = [-1.0];
	
	b = A*[0.0];
	b_bar = A_bar*[0.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,1);
	
	vars = class_variables(lp);
	
	ip_algorithm(lp, settings, vars);
end

function trival_unbounded2(settings)
	println("trival_unbounded2")
	println("min -x")
	println("eye(1)*x >= 0")
	
	A_bar = spzeros(0,1); 
	
	A = speye(1,1);
	
	c = [-1.0];
	
	b = [0.0];
	b_bar = A_bar*[0.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,0);
	
	vars = class_variables(lp);
	
	ip_algorithm(lp, settings, vars);
end

function lp_feasible_and_bounded_1(settings)
	println("lp_feasible_and_bounded_1")
	println("min x[1] + 2*x[2] + 3*x[3] + 4*x[4]")
	println("x[1] + x[2] + x[3] + x[4] = 4.0")
	println("x[1], x[2], x[3], x[4] >= 0")
	
	x0 = [1; 1; 1; 1.0];

	A_bar = sparse([1.0 1.0 1.0 1.0]); 
	
	A = spzeros(0,4);
	
	c = [1, 2, 3, 4.0];
	
	b = A*x0;
	b_bar = [4.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,4);
	
	vars = class_variables(lp);
	
	ip_algorithm(lp, settings, vars);
	
	# ans
	# x = [4.0, 0.0, 0.0, 0.0]
	
end

function lp_feasible_and_bounded_2(settings)
	println("lp_feasible_and_bounded_2")
	println("min x[1] + x[2]")
	println("x[1] + x[2] >= 2")
	println("x[1], x[2] >= 0")
	
	x0 = [1; 1.0];

	A_bar = spzeros(0,2); 
	
	A = sparse([1.0 1.0]);
	
	c = [1.0, 0.0];
	
	b = [2.0];
	b_bar = A_bar*x0;
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,2);
	
	vars = class_variables(lp);
	
	ip_algorithm(lp, settings, vars);
	
	
end


function lp_feasible_and_bounded_3(settings)
	println("lp_feasible_and_bounded_3")
	println("min x[1] + 2*x[2] + 3*x[3] + 4*x[4]")
	println("x[1] + x[2] + x[3] + x[4] == 1.0")
	println("x[1] + x[2] + x[3] + x[4] >= 1.0")
	println("x[1], x[2], x[3], x[4] >= 0")
	
	x0 = [1; 1; 1; 1.0];

	A_bar = sparse([1.0 1.0 1.0 1.0]); 
	
	A = sparse([1.0 1.0 1.0 1.0]);
	
	c = [1, 2, 3, 4.0];
	
	b = [1.0];
	b_bar = [1.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,4);
	
	vars = class_variables(lp);
	
	ip_algorithm(lp, settings, vars);
	
	# ans
	# x = [4.0, 0.0, 0.0, 0.0]
	
end

function lp_feasible_and_bounded_4(settings)
	println("lp_feasible_and_bounded_4")
	println("x[1] + x[2] + x[3] = 3")
	println("eye(3)*x >= 0")
	x0 = [1; 1; 1.0];

	A_bar = sparse([1.0 1.0 1.0]); 
	
	A = speye(3);
	
	c = [1, 2, 3.0];
	
	b = [0.0, 0.0, 0.0];
	b_bar = [3.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,0);
	
	vars = class_variables(lp);
	
	
	ip_algorithm(lp, settings, vars);
	
	# ans
	# x = [4.0, 0.0, 0.0, 0.0]
	
end


function lp_blend_problem_1(settings)
	try
		println("lp_blend_problem_1")
		println("Uses equality constraints and x >= 0")
		
		file = matopen("Problems/blend.mat")
		A_bar = sparse(read(file, "A"))
		n,k = size(A_bar)	
		b_bar = read(file,"b")
		b_bar = b_bar[:]
		c = read(file,"c")
		c = c[:];
		
		println("A matrix is size " * string((n,k)) * " and has " * string(nnz(A_bar))  * " non-zeros")
		
		lp = class_non_linear_program();
		lp.set_linear_objective(c);
		lp.set_linear_constraints(spzeros(0,k),spzeros(0,1)*[1.0],A_bar,b_bar,k);
		
		vars = class_variables(lp);
		
		ip_algorithm(lp, settings, vars);
	catch	e
		println("ERROR in lp_blend_problem_1")
		throw(e)
	end
end

function lp_blend_problem_2(settings)
	try
		println("lp_blend_problem_2");
		file = matopen("Problems/blend.mat");
		A_bar = sparse(read(file, "A"));
		n,k = size(A_bar);	
		b_bar = read(file,"b");
		b_bar = b_bar[:];
		c = read(file,"c")
		c = c[:];
		
		println("A matrix is size " * string((n,k)) * " and has " * string(nnz(A_bar))  * " non-zeros")
		
		lp = class_non_linear_program();
		lp.set_linear_objective(c);
		lp.set_linear_constraints(speye(k),spzeros(k,1)*[1.0],A_bar,b_bar,0);
		
		vars = class_variables(lp);
		
		ip_algorithm(lp, settings, vars);
	catch	e
		println("ERROR in lp_blend_problem_1")
		throw(e)
	end
end


function lp_blend_problem_infeasible_1(settings)
	try
		println("lp_blend_problem_infeasible_1")
		file = matopen("Problems/blend.mat")
		A_bar = sparse(read(file, "A"))
		n,k = size(A_bar)	
		b_bar = read(file,"b")
		b_bar = b_bar[:]
		c = read(file,"c")
		c = c[:];
		
		println("A matrix is size " * string((n,k)) * " and has " * string(nnz(A_bar))  * " non-zeros")
		
		A_bar = sparse([A_bar; c']);
		b_bar = [b_bar; -35.0];
		
		lp = class_non_linear_program();
		lp.set_linear_objective(c);
		lp.set_linear_constraints(spzeros(0,k),spzeros(0,1)*[1.0],A_bar,b_bar,k);
		
		vars = class_variables(lp);
		
		ip_algorithm(lp, settings, vars);
	catch	e
		println("ERROR in lp_blend_problem_1")
		throw(e)
	end
end

function lp_blend_problem_unbounded_1(settings)
	try
		println("lp_blend_problem_unbounded_1")
		file = matopen("Problems/blend.mat")
		A_bar = sparse(read(file, "A"))
		n,k = size(A_bar)	
		b_bar = read(file,"b")
		b_bar = b_bar[:]
		#c = read(file,"c")
		c = -ones(k);
		
		println("A matrix is size " * string((n,k)) * " and has " * string(nnz(A_bar))  * " non-zeros")
		
		A_bar = sparse([A_bar]);
		b_bar = [b_bar];
		
		lp = class_non_linear_program();
		lp.set_linear_objective(c);
		lp.set_linear_constraints(spzeros(0,k),spzeros(0,1)*[1.0],A_bar,b_bar,k);
		
		vars = class_variables(lp);
		
		ip_algorithm(lp, settings, vars);
	catch	e
		println("ERROR in lp_blend_problem_1")
		throw(e)
	end
end

function lp_unbounded_1()

end

###########################################
# linear programs with convex objectives
###########################################


function qp_feasible_and_bounded_1(settings)
	println("qp_feasible_and_bounded_1")
	x0 = [1; 1; 1.0];

	A_bar = sparse([1.0 1.0 1.0]); 
	
	A = speye(0,3);
	
	c = [1.0, 2.0, 3.0];
	Q = speye(3)
	
	b = A*x0;
	b_bar = A_bar*x0;
	
	lp = class_non_linear_program()  ;
	lp.set_quadratic_objective(c,Q);
	lp.set_linear_constraints(A,b,A_bar,b_bar,3); 
	
	vars = class_variables(lp);
	
	
	ip_algorithm(lp, settings, vars);
end


function qp_infeasible_1(settings)
	println("qp_infeasible_1")
	x0 = [1; 1; 1.0];

	A_bar = sparse([1.0 1.0 1.0]); 
	
	A = speye(0,3);
	
	c = [1, 2, 3.0];
	Q = speye(3)
	
	b = A*x0;
	b_bar = [-1.0];
	
	lp = class_non_linear_program()  ;
	lp.set_quadratic_objective(c,Q);
	lp.set_linear_constraints(A,b,A_bar,b_bar,3); 
	
	vars = class_variables(lp);
	
	
	ip_algorithm(lp, settings, vars);
end

function qp_unbounded_1(settings)
	println("qp_unbounded_1")
	x0 = [1; 1; 1.0];

	A_bar = sparse([1.0 1.0 1.0]); 
	
	A = speye(0,3);
	
	c = [-1.0, 2.0, 3.0];
	Q = speye(3)
	
	b = A*x0;
	b_bar = [1.0];
	
	lp = class_non_linear_program()  ;
	lp.set_quadratic_objective(c,Q);
	lp.set_linear_constraints(A,b,A_bar,b_bar,3); 
	
	vars = class_variables(lp);
	
	ip_algorithm(lp, settings, vars);
end


function qp_blend_1(settings)
	
end

###########################################
# linear program with L2 norm constraint
###########################################

# need to modify line search !!!!

function L2_problem(settings,b)
	# min -x[2]
	# - x[1] -x[2]^2 + b == 0 (x[2]^2 <= b)
	# x[1], x[2] >= 0
	println("L2_problem with b = ",b)	
	println("i.e. min -x[2]")
	println("s.t. x[1] + x[2]^2 <= ",b)
	println("x[1],x[2] >= 0")
	
	settings.max_iter = 50;
	
	nlp = class_non_linear_program();
	nlp.set_linear_objective([0.0, -1.0]);
	
	nlp.evaluate_constraints = function(x_scaled)
		return [-x_scaled[2]^2 - x_scaled[1] + b];
	end
	
	nlp.evaluate_constraint_gradients = function(x_scaled)
		return sparse([[-1 -2*x_scaled[2]]]);
	end
	
	nlp.evaluate_constraint_lagrangian_hessian = function(x_scaled,y_scaled)
		return sparse([[0; 0] [0; -2.0]])*y_scaled[1]
	end
	
	nlp.n_1 = 2;
	nlp.n_2 = 0;
	nlp.n = nlp.n_1 + nlp.n_2;
	nlp.m_1 = 0;
	nlp.m_2 = 1;
	nlp.m = nlp.m_1 + nlp.m_2;
	
	
	vars = class_variables(nlp);
	
	ip_algorithm(nlp, settings, vars);
	
	println(vars.kappa)
	println(vars.y_scaled)
end


function L2_problem_inequality(settings,b)
	# min -x
	# x^2 <= b
	# x >= 0
	println("L2_problem with b = ",b)	
	println("i.e. min -x")
	println("s.t. x^2 <= ",b)
	println("x >= 0")
	
	settings.max_iter = 50;
	
	nlp = class_non_linear_program();
	nlp.set_linear_objective([-1.0]);
	
	nlp.evaluate_constraints = function(x_scaled)
		return [-x_scaled[1]^2 + b];
	end
	
	nlp.evaluate_constraint_gradients = function(x_scaled)
		return sparse([[-2*x_scaled[1]]]);
	end
	
	nlp.evaluate_constraint_lagrangian_hessian = function(x_scaled,y_scaled)
		return sparse([[-2.0]])*y_scaled[1]
	end
	
	nlp.n_1 = 1;
	nlp.n_2 = 0;
	nlp.n = nlp.n_1 + nlp.n_2;
	nlp.m_1 = 1;
	nlp.m_2 = 0;
	nlp.m = nlp.m_1 + nlp.m_2;
	
	
	vars = class_variables(nlp);
	
	ip_algorithm(nlp, settings, vars);
	
	println(vars.kappa)
	println(vars.y_scaled)
end

###########################################
# non-convex objective
###########################################

###########################################
# non-convex constraint
###########################################



###########################################
# main
###########################################
main()
