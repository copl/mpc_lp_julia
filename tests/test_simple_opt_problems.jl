# this file contains a series of small and simple tests optimization problems
# 

println("loading libraries")


#using JuMP
#using Gurobi
#using Ipopt
#using KNITRO
#using SCS
using MAT

println("external libraries loaded")

include("testing_tools.jl")

println("internal libraries loaded")

function main()
	settings = standard_settings()

	println("------------------------")
	println("Start tests")
	println("------------------------")
	
	if true
	
		# test problems
		trival_problem1(settings)
		trival_problem2(settings)
		trival_problem3(settings)
		
		trival_infeasible1(settings)
		trival_infeasible2(settings)
		trival_unbounded1(settings)
		trival_unbounded2(settings)
		
		#lp_feasible_and_bounded_1(settings)
		#lp_feasible_and_bounded_2(settings)
		#lp_feasible_and_bounded_3(settings)
		#lp_feasible_and_bounded_4(settings)
		
		qp_feasible_and_bounded_1(settings)
		qp_feasible_and_bounded_2(settings)
		qp_infeasible_1(settings)
		qp_unbounded_1(settings)
		qp_unbounded_and_infeasible_1(settings)
		
		#lp_blend_problem_1(settings)
		#lp_blend_problem_2(settings)
		#lp_blend_problem_infeasible_1(settings)
		#lp_blend_problem_unbounded_1(settings)
		#qp_blend_problem_1(settings)
	end
	
	if true
		L2_problem(settings, 1.0,1)
		L2_problem(settings, 0.0,1)
		L2_problem(settings, -1.0,2)
		L2_problem(settings, -10.0,2)
		L2_problem_inequality(settings,-1.0,2)
	end
	
	if true
		trival_non_convex_quadratic_problem1(settings)
		trival_non_convex_quadratic_problem2(settings)
		trival_non_convex_quadratic_problem3(settings)
	end
end


############################
# linear programs
############################

# A = [[1,1],[1,1]]


function trival_problem1(settings)
	A_bar = spzeros(0,1); 
	
	A = spzeros(0,1);
	
	c = [1.0];
	
	b = EMPTY_ARRAY;
	b_bar = EMPTY_ARRAY;
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(sparse(A),b,sparse(A_bar),b_bar,1);
	
	vars = class_variables(lp);
	
	vars, status = ip_algorithm(lp, settings, vars, false);
	details = "min x" * "\n" * "x >= 0";
	is_problem_successful("trival_problem1", status, 1, details)
	
end

function trival_problem2(settings)
	A_bar = spzeros(0,2); 
	
	A = spzeros(0,2);
	
	c = [1.0, 0.0];
	
	b = EMPTY_ARRAY;
	b_bar = EMPTY_ARRAY;
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(sparse(A),b,sparse(A_bar),b_bar,2);
	
	vars = class_variables(lp);
	
	vars, status = ip_algorithm(lp, settings, vars, false);
	details = "min x[1]" * "\n" * "x[1], x[2] >= 0";
	is_problem_successful("trival_problem2", status, 1, details)
end

function trival_problem3(settings)
	A_bar = spzeros(0,2); 
	
	A = speye(2);
	
	c = [1.0, 0.0];
	
	b = [0.0, 0.0];
	b_bar = EMPTY_ARRAY;
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,0);
	
	vars = class_variables(lp);
	
	vars, status = ip_algorithm(lp, settings, vars, false);
	details = "min" * string(c) * " *x" * "\n" * string(full(A)) * " * x >= ", string(b);
	is_problem_successful("trival_problem3", status, 1, details)
	
end

function trival_infeasible1(settings)
	A_bar = speye(1); 
	
	A = spzeros(0,1);
	
	c = [1.0];
	
	b_bar = [-1.0];
	b = A*[0.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,1);
	
	vars = class_variables(lp);
	
	vars, status = ip_algorithm(lp, settings, vars, false);
	details = "min x" * "\n" * "x == -1" * "\n" * "x >= 0"
	is_problem_successful("trival_infeasible1", status, 2, details)
end

function trival_infeasible2(settings)
	A_bar = spzeros(0,1); 
	
	A = -speye(1);
	
	c = [1.0];
	
	b = [1.0];
	b_bar = A_bar*[0.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,1);
	
	vars = class_variables(lp);
		
	vars, status = ip_algorithm(lp, settings, vars, false);
	details = "min x" * "\n" * "x <= -1" * "\n" * "x >= 0"
	is_problem_successful("trival_infeasible2", status, 2, details)
end

function trival_unbounded1(settings)
	A_bar = spzeros(0,1); 
	
	A = spzeros(0,1);
	
	c = [-1.0];
	
	b = A*[0.0];
	b_bar = A_bar*[0.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,1);
	
	vars = class_variables(lp);
	
	vars, status = ip_algorithm(lp, settings, vars, false);
	details = "min -x" * "\n" * "x >= 0"
	is_problem_successful("trival_unbounded1", status, 3, details)
end

function trival_unbounded2(settings)
	A_bar = spzeros(0,1); 
	
	A = speye(1,1);
	
	c = [-1.0];
	
	b = [0.0];
	b_bar = A_bar*[0.0];
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(A,b,A_bar,b_bar,0);
	
	vars = class_variables(lp);
	
	vars, status = ip_algorithm(lp, settings, vars, false);
	is_problem_successful("trival_unbounded1", status, 3, None)
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
	catch e
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
	A_bar = speye(0,3); 
	A = sparse([1.0 1.0 1.0]);
	b = [1.0];
	b_bar = EMPTY_ARRAY;
	
	c = [0.0, 2.0, 3.0];
	Q = spdiagm([1.0, 1.0, 1.0])
	
	qp = class_non_linear_program() ;
	qp.set_quadratic_objective(c,Q);
	qp.set_linear_constraints(A,b,A_bar,b_bar,3); 
	
	vars = class_variables(qp);
	
	vars, status = ip_algorithm(qp, settings, vars, false);
		
	is_problem_successful("qp_feasible_and_bounded_1", status, 1, None)
end

function qp_feasible_and_bounded_2(settings)
	A_bar = speye(0,3); 
	A = sparse([1.0 1.0 1.0]);
	b = [1.0];
	b_bar = EMPTY_ARRAY;
	
	c = [-1.0, 2.0, 3.0];
	Q = spdiagm([1.0, 1.0, 1.0])
	
	qp = class_non_linear_program() ;
	qp.set_quadratic_objective(c,Q);
	qp.set_linear_constraints(A,b,A_bar,b_bar,3); 
	
	vars = class_variables(qp);
	
	vars, status = ip_algorithm(qp, settings, vars, false);
		
	is_problem_successful("qp_feasible_and_bounded_2", status, 1, None)
end


function qp_infeasible_1(settings)
	A_bar = speye(0,3); 
	A = -sparse([1.0 1.0 1.0]);
	b = [1.0];
	b_bar = EMPTY_ARRAY;
	
	c = [-1.0, 2.0, 3.0];
	Q = spdiagm([1.0, 1.0, 1.0])
	
	qp = class_non_linear_program() ;
	qp.set_quadratic_objective(c,Q);
	qp.set_linear_constraints(A,b,A_bar,b_bar,3); 
	
	vars = class_variables(qp);
	
	vars, status = ip_algorithm(qp, settings, vars, false);
		
	is_problem_successful("qp_infeasible_1", status, 2, None)
end

function qp_unbounded_1(settings)
	A_bar = speye(0,3); 
	A = sparse([1.0 1.0 1.0]);
	b = [1.0];
	b_bar = EMPTY_ARRAY;
	
	c = [-1.0, 2.0, 3.0];
	Q = spdiagm([0.0, 1.0, 1.0])
	
	qp = class_non_linear_program() ;
	qp.set_quadratic_objective(c,Q);
	qp.set_linear_constraints(A,b,A_bar,b_bar,3); 
	
	vars = class_variables(qp);
	
	vars, status = ip_algorithm(qp, settings, vars, false);
		
	is_problem_successful("qp_unbounded_1", status, 3, None)
end

function qp_unbounded_and_infeasible_1(settings)
	A_bar = speye(0,3); 
	A = -sparse([1.0 1.0 1.0]);
	b = [1.0];
	b_bar = EMPTY_ARRAY;
	
	c = [-1.0, 2.0, 3.0];
	Q = spdiagm([0.0, 1.0, 1.0])
	
	qp = class_non_linear_program() ;
	qp.set_quadratic_objective(c,Q);
	qp.set_linear_constraints(A,b,A_bar,b_bar,3); 
	
	vars = class_variables(qp);
	
	vars, status = ip_algorithm(qp, settings, vars, false);
	
	correct_status = 2;
	if status == 3
		correct_status = 3; # unbounded is also a valid status
	end
	
	is_problem_successful("qp_unbounded_and_infeasible_1", status, correct_status, None)
end


function qp_blend_1(settings)
	
end

###########################################
# linear program with L2 norm constraint
###########################################

# need to modify line search !!!!

function L2_problem(settings, b, correct_status)
	# min -x[2]
	# - x[1] -x[2]^2 + b == 0 (x[2]^2 <= b)
	# x[1], x[2] >= 0
	
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
	
	vars, status = ip_algorithm(nlp, settings, vars, false);
		
	details = "min -x[2]" * "\n" * "s.t. x[1] + x[2]^2 <= " * string(b) * "\n" * "x[1],x[2] >= 0";
	
	is_problem_successful("L2_problem with b = " * string(b), status, correct_status, details)
end


function L2_problem_inequality(settings,b,correct_status)
	# min -x
	# x^2 <= b
	# x >= 0
		
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
	
	vars, status = ip_algorithm(nlp, settings, vars, false);
		
	details = "min -x" * "\n" * "s.t. x^2 <= " * string(b) * "\n" * "x >= 0";
	
	is_problem_successful("L2_problem with b = " * string(b), status, correct_status, details)
	
end

###########################################
# non-convex objective
###########################################

function trival_non_convex_quadratic_problem1(settings)
	
	
	A_bar = spzeros(0,1); 
	
	A = -sparse(ones(1,1));
		
	b = [-1.0];
	b_bar = EMPTY_ARRAY;
	
	qp = class_non_linear_program();
	
	qp.set_quadratic_objective([1.0],sparse([[-1.0]]));
	qp.set_linear_constraints(sparse(A),b,sparse(A_bar),b_bar,1);
	
	vars = class_variables(qp);
	
	vars, status = ip_algorithm(qp, settings, vars, true);
	
	
	details = "min x-x^2" * "\n" * "s.t. 1 >= x >= 0";
	
	is_problem_successful("trival_non_convex_quadratic_problem1",status,1,details)
end

function trival_non_convex_quadratic_problem2(settings)
	
	
	A_bar = spzeros(0,1); 
	
	A = -sparse(ones(1,1));
		
	b = [-1.0];
	b_bar = EMPTY_ARRAY;
	
	qp = class_non_linear_program();
	
	qp.set_quadratic_objective([0.0],sparse([[-1.0]]));
	qp.set_linear_constraints(sparse(A),b,sparse(A_bar),b_bar,1);
	
	vars = class_variables(qp);
	
	vars, status = ip_algorithm(qp, settings, vars, true);
	
	
	details = "min -x^2" * "\n" * "s.t. 1 >= x >= 0";
	
	is_problem_successful("trival_non_convex_quadratic_problem2",status,1,details)
end

function trival_non_convex_quadratic_problem3(settings)
	
	
	A_bar = spzeros(0,1); 
	
	A = -sparse(ones(1,1));
		
	b = [-1.0];
	b_bar = EMPTY_ARRAY;
	
	qp = class_non_linear_program();
	
	qp.set_quadratic_objective([0.0],sparse([[-10.0]]));
	qp.set_linear_constraints(sparse(A),b,sparse(A_bar),b_bar,1);
	
	vars = class_variables(qp);
	
	vars, status = ip_algorithm(qp, settings, vars, true);
	
	
	details = "min -10*x^2" * "\n" * "s.t. 1 >= x >= 0";
	
	is_problem_successful("trival_non_convex_quadratic_problem3",status,1,details)
end

###########################################
# non-convex constraint
###########################################



###########################################
# main
###########################################
main()
