# TO DO
# UPDATE primal feasibility to work with non-linear constraints
# add merit function
# 

println("loading libraries")

function debug_message(message)
	#println(message)
end
function debug_message_high_level(message)
	#println(message)
end

using MAT
using JuMP
#using Gurobi
using Ipopt
using KNITRO
#using SCS
println("external libraries loaded")

include("interior_point_algorithm.jl")

println("internal libraries loaded")

function main() # this function is called at the bottom of the code
	
	settings = class_settings();
	
	settings.max_iter = 100;  # Total number of iterarions
	settings.bkscale             = 0.95;  # Back scaling for line search.
	
	settings.primal_feas_tol = 1e-8
	settings.dual_feas_tol = 1e-8
	settings.duality_gap_tol = 1e-8
	settings.primal_infeas_tol = 1e-8
	settings.dual_infeas_tol = 1e-8
	
	
    #linear_feas_tol     = 1e-10;  # Threshold for feasibility of the linear constraints.
    #comp_tol            = 1e-10;  # Threshold for complementry optimality condition s^Tz
    

	# Initialize configuration variable
    

	#nonlinear_random_points_test(settings)
	lp_test(settings)
	#srand(123)
	# We are creating an instance of LP. In practice, we should read the problem data from input stream.
	
end

function lp_test(settings)
	problem_data = construct_instance_from_mat("Problems/blend.mat");
	problem_data.barrier_weights = ones(problem_data.m);# - 1.0*log(rand(problem_data.m))
	#problem_data.c = -ones(problem_data.k);
	
	#di = rand(problem_data.k)
	#di[1:10] = 0;
	#problem_data.P = spdiagm(di)
	
	#problem_data.add_inequality_constraint(ones(problem_data.k)',1000,1.0)
	#add_nonlinear_objective3(problem_data)
	#problem_data.add_inequality_constraint(problem_data.c',-30.84,1.0)
	add_nonlinear_constraints(problem_data)
	
	add_nonlinear_constraints(problem_data)
	
	
	variables = class_linear_program_variables(problem_data);
	variables, status = @time interior_point_algorithm(problem_data,settings,variables);
	
	#solve_with_JuMP(problem_data)
	
	print_status(status,variables,problem_data)
end

function print_status(status,variables,problem_data)
	if status["primal_feasible"] == true && status["dual_feasible"] == true
		xOur = variables.x/variables.tau;
		value = problem_data.currentObjectiveValue;
		println("Our objective: " * string(value))
		println("Our feasiblity: " * string(calculate_feasibility(problem_data,xOur)))
		println("Our x norm: " * string(norm(xOur,2)))
	elseif length(status) == 2
		if status["primal_feasible"] == false
			println("primal locally infeasible")
		end
		if status["dual_feasible"] == false
			println("dual locally infeasible")
		end
		if status["primal_feasible"] == false
			println("i.e. linear approximation of constraints at current point is infeasible")
		elseif status["dual_feasible"] == false
			println("i.e. problem is locally unbounded")
		end
	elseif status == -1
		println("max iterations exceeded")
	else
		println("status unknown")
	end
end

function nonlinear_random_points_test(settings)
	srand(12)


	value_list = []
	jump_value_list = []
	
	problem_data = construct_instance_from_mat("Problems/blend.mat")
	problem_data.add_inequality_constraint(ones(problem_data.k)',500,1.0)
	problem_data.add_inequality_constraint(rand(problem_data.k)',200,1.0)
	problem_data.add_inequality_constraint(rand(problem_data.k)',200,1.0)
	problem_data.add_inequality_constraint(rand(problem_data.k)',200,1.0)
	add_nonlinear_objective3(problem_data)
	
	for i = 1:50 #100
		if i > 1
			#add_nonlinear_constraints(problem_data)
			#problem_data.barrier_weights = (ones(problem_data.m)  - 3.0*log(rand(problem_data.m)))/3.0
		end
		#
		if true
			xTrue = solve_with_JuMP(problem_data)
			println(norm(min(xTrue,0)))
			
			jump_value = problem_data.evaluate_objective(xTrue);
			jump_value_list = [jump_value_list, jump_value];
			
			println("JuMP objective: " * string(jump_value))
			println("JuMP feasibility: " * string(calculate_feasibility(problem_data,xTrue)))
		end
		
		if false
			# The main function that run interior point algorithm.
			variables = class_linear_program_variables(problem_data);
			variables.x = xTrue;
			variables.s = [xTrue + 0.1, 100.0];
			variables.kappa = 0.1;
			variables, status = @time interior_point_algorithm(problem_data,settings,variables);
			
			#variables.s = xTrue + 0.1;
			#variables.kappa = 0.1
			#problem_data.update(variables)
			#problem_data.black_box_objective = false;
			#problem_data.black_box_constraints = false;
			
			#variables = class_linear_program_variables(problem_data);
			
			#variables, status = @time interior_point_algorithm(problem_data,settings,variables);
			
			#problem_data.black_box_objective = true;
			#problem_data.black_box_constraints = true;
		end
		
		variables = class_linear_program_variables(problem_data);
		variables, status = @time interior_point_algorithm(problem_data,settings,variables);
		
		if status == 1
			xOur = variables.x/variables.tau;
			value = problem_data.evaluate_objective(xOur);
			println("Our objective: " * string(value))
			value_list = [value_list, value]
			println("Our feasiblity: " * string(calculate_feasibility(problem_data,xOur)))
		elseif status == 2
				println("infeasible or unbounded")
		elseif status == 0
			println("max iterations")
		else
			println("status unknown")
		end
		print_status(status)
	end
	
	println(jump_value_list)
	println(value_list)
	#println("Smallest jump value" * string(minimum(jump_value_list)))
	println("Unweighted value: " * string(value_list[1]))
	println("Smallest value: " * string(minimum(value_list)))
	println("Smallest value jump: " * string(minimum(jump_value_list)))
	#println(variables.x/variables.tau)
	
	#println("Distance to JuMP solver solution: "* string(norm(xTrue - xOur)))
end

function get_objective_value(problem_data::class_linear_program_input,x)
	return(problem_data.c' * x + 0.5*x'*problem_data.P*x)
end

function calculate_feasibility(problem_data::class_linear_program_input,x)
	return norm([problem_data.A*x - problem_data.b, max(problem_data.G*x - problem_data.h,0)])
end

function interpret_variables(variables::class_linear_program_variables)
	if variables.tau > variables.kappa
		println("Optimal solution found")
	else
		println("Problem is infeasible or unbounded")
	end
end

function solve_with_JuMP(problem_data::class_linear_program_input)
	#model = Model(solver=GurobiSolver())
	model = Model(solver=IpoptSolver())
	#model = Model(solver=KnitroSolver()) #ms_enable=1
	#model = Model(solver=SCSSolver())
	k = problem_data.k
	m = problem_data.m
	n = problem_data.n
	
	A = problem_data.A
	G = problem_data.G
	c = problem_data.c
	b = problem_data.b
	h = problem_data.h
	
	@defVar(model, x[1:k], start = rand()*20) #rand()*2000 )
	#@setObjective(model, Min, sum{c[j]*x[j] + -0.001*x[j]*x[j] , j=1:k} ) # + -2*0.5*x[j]*x[j]
	@setNLObjective( model, Min, sum{ c[j]*x[j] + x[j]*x[j] , j=1:k} ) # + -2*0.5*x[j]*x[j] 
	# 0.000001*(-(x[j]-5.0)^2)
	
	#+ 0.01*x[j]^2
	#Afull = full(A)
	#indicies = M.rowval[a.colptr[col] : M.colptr[col+1]-1] 

	for i = 1:n
		@addConstraint(model, sum{A[i,j]*x[j], j=1:k} == b[i])
	end
	
	for i = 1:m
		@addConstraint(model, sum{G[i,j]*x[j], j=1:k} <= h[i])
	end
	
	@addNLConstraint(model, sum{x[j]*x[j], j=1:k} <= 10000)
	
	status = solve(model)

	println("Objective value: ", getObjectiveValue(model))
	
	# convert to standard array
	xJuMP = getValue(x)
	xTrue = zeros(k)
	
	for i = 1:k
		xTrue[i] = xJuMP[i]
	end
	
	return(xTrue)
end


# Constructing a random LP
function construct_instance1()
	n = 5;
	k = 10;
	m = 10;
	
	problem_data = class_linear_program_input()
	x0 = rand(k,1)
	
	A = rand(n, k);
	G = -diagm(ones(m));
	c = rand(k,1);
	h = zeros(m);
	b = A*x0;
	
	
	problem_data.A = A # equality constraints
	problem_data.G = G # inequality constraints
	problem_data.c = c
	problem_data.h = h
	problem_data.b = b
	problem_data.m = m
	problem_data.k = k
	problem_data.n = n
	
	return(problem_data)
end

function construct_instance2() # example synced with c++
	m = 2;
	n = 2; #3;
	k = 4;
	
	x0 = [1; 1; 1; 1];

	A = [1 0 0 0;
		0 1 0 0]; 
	
	G = [1 0 1 0;
		0 1 0 1];
	
	c = [1, 2, 3, 4];
	
	b = A*x0;
	h = G*x0;
	
	#[x,fmin]=linprog(c,G,h,A,b)

	problem_data = class_linear_program_input()
	problem_data.A = A
	problem_data.G = G
	problem_data.c = c
	problem_data.h = h
	problem_data.b = b
	problem_data.m = m
	problem_data.k = k
	problem_data.n = n
	
	return(problem_data)
end

function construct_instance3()
	f = open("example_problems/ex3sta1/ex3sta1.mtx")
	
	size_found = false
	A = false
	for ln in eachline(f)
		if (length(ln) > 0 && ln[1] != '%')
			sp = split(ln)
			if ~size_found
				n = int(sp[1])
				m = int(sp[2])
				nonzeros = int(sp[3])
				
				A = spzeros(n, m)
				size_found = true
			end
			
			if size_found
				i = int(sp[1])
				j = int(sp[2])
				val = float(sp[3])
				#print(i)
				#print(j)
				#print(val)
				A[i,j] = val
			end
		end
	end
end

function construct_instance_from_mat(filenamelp)
	println("Reading Matlab file for lp input " * filenamelp )

	file = matopen(filenamelp)
	A = sparse(read(file, "A"))
	n,k = size(A)
	m = k
	G = -speye(m)
	
	c = read(file,"c")
	b = read(file,"b")
	h = zeros(m,1)
	
	#max_norm = 0;
	#min_norm = 999999;
	#for j = 1:k
	#	max_norm = max(max_norm,norm(A[:,j]))
	#	min_norm = min(min_norm,norm(A[:,j]))
	#end
	#println(max_norm)
	#println(min_norm)
	
	println("A matrix is size " * string((n,k)) * " and has " * string(nnz(A))  * " non-zeros")
	
	problem_data = class_linear_program_input()
	
	problem_data.m = m; 
	problem_data.k = k
	problem_data.n = n
	
	problem_data.P = spzeros(k,k);
	problem_data.A = sparse(A)
	problem_data.G = sparse([G]) #[G; G])
	problem_data.c = c[:];
	problem_data.h = h[:]; #; h[:]];
	problem_data.b = b[:];
	
	#println(size(problem_data.h))
	
	#println(problem_data.m)
	problem_data.barrier_weights = ones(problem_data.m);
	
	return(problem_data)
end


function add_convex_objective(problem_data)
	k = problem_data.k
		
	#_c = problem_data.c;
	objective = function(x)
		return 1.0*()[1] # - 0.0001*sum(x.^2)
		#return (_c' * (x)) [1] 
	end
	
	gradient = function(x)
		xPos = max(x,0)
		return 1.0*(-ones(k) + 0.02*x  + 0.1*0.5*(x+1).^(-0.5)) #  - 2*0.0001*x
		#return _c
	end
	
	hessian = function(x)
		di = 0.02*ones(k) + -0.1*0.25*(x+1).^(-1.5); 
		return 1.0*spdiagm(di)  # - minimum([0,di]) + 0.0001 # #spdiagm( -0.25*(variables.x + 1).^(-1.5) ) #spzeros(length(variables.x),length(variables.x))#
		#return spzeros(length(x),length(x))
	end
	
	problem_data.add_black_box_objective(objective,gradient,hessian)
end

function add_nonlinear_objective(problem_data)
	k = problem_data.k
	
	#_c = problem_data.c;
	objective = function(x)
		return 1.0*(-ones(k)'*x + 0.01*sum(x.^2) + 0.1*sum((x+1).^(0.5)))[1] # - 0.0001*sum(x.^2)
		#return (_c' * (x)) [1] 
	end
	
	gradient = function(x)
		xPos = max(x,0)
		return 1.0*(-ones(k) + 0.02*x  + 0.1*0.5*(x+1).^(-0.5)) #  - 2*0.0001*x
		#return _c
	end
	
	hessian = function(x)
		di = 0.02*ones(k) + -0.1*0.25*(x+1).^(-1.5); 
		return 1.0*spdiagm(di)  # - minimum([0,di]) + 0.0001 # #spdiagm( -0.25*(variables.x + 1).^(-1.5) ) #spzeros(length(variables.x),length(variables.x))#
		#return spzeros(length(x),length(x))
	end
	
	problem_data.add_black_box_objective(objective,gradient,hessian)
end

function add_nonlinear_objective2(problem_data)
	k = problem_data.k
	
	#_c = problem_data.c;
	c = -rand(k) + problem_data.c;
	problem_data.c = c;
	
	objective = function(x)
		return 1.0*(c'*x + 0.1*sum((x+1).^(0.5)))[1] # - 0.0001*sum(x.^2)
		#return (_c' * (x)) [1] 
	end
	
	gradient = function(x)
		return 1.0*(c + 0.1*0.5*(x+1).^(-0.5)) #  - 2*0.0001*x
		#return _c
	end
	
	hessian = function(x)
		di = -0.1*0.25*(x+1).^(-1.2); 
		return 0.001*speye(k,k)#1.0*spdiagm(di)  # - minimum([0,di]) + 0.0001 # #spdiagm( -0.25*(variables.x + 1).^(-1.5) ) #spzeros(length(variables.x),length(variables.x))#
		#return spzeros(length(x),length(x))
	end
	
	problem_data.add_black_box_objective(objective,gradient,hessian)
end

function add_nonlinear_objective3(problem_data)
	k = problem_data.k
	
	#_c = problem_data.c;
	c = -rand(k) + problem_data.c;
	problem_data.c = c;
	
	objective = function(x)
		return 1.0*(sum(-(x-5.0).^2))[1] # - 0.0001*sum(x.^2)
		#return (_c' * (x)) [1] 
	end
	
	gradient = function(x)
		return -2.0*(x - 5.0) #  - 2*0.0001*x
		#return _c
	end
	
	hessian = function(x)
		di = -2; 
		return 1.0*speye(k,k) #1.0*spdiagm(di)  # - minimum([0,di]) + 0.0001 # #spdiagm( -0.25*(variables.x + 1).^(-1.5) ) #spzeros(length(variables.x),length(variables.x))#
		#return spzeros(length(x),length(x))
	end
	
	problem_data.add_black_box_objective(objective,gradient,hessian)
end

function add_nonlinear_constraints(problem_data)
	
	_h = copy(problem_data.h);
	evaluate_inequality_constraint = function(x)
		#return [problem_data.G[1:(problem_data.m-1),:]*x; ones(problem_data.k)'*x] - [_h[1:problem_data.m-1]; 10];
		return [problem_data.G[1:(problem_data.m-1),:]*x; x'*x] - [_h[1:problem_data.m-1]; 10000];  #  - 2*0.0001*x
		#  -0.01*(x'*x)[1] * ones(problem_data.m) + 
		#return _c
	end
	
	evaluate_inequality_constraint_gradients = function(x)
		#return sparse([problem_data.G[1:(problem_data.m-1),:]; ones(problem_data.k)'])
		return sparse([problem_data.G[1:(problem_data.m-1),:]; 2*x']); # - minimum([0,di]) + 0.0001 # #spdiagm( -0.25*(variables.x + 1).^(-1.5) ) #spzeros(length(variables.x),length(variables.x))#
		#-0.01*2*ones(problem_data.k)*x' + 
		#return spzeros(length(x),length(x))
	end
	
	evaluate_inequality_constraint_lagrangian = function(x,z)
		#return spzeros(problem_data.k,problem_data.k)
		return 2*speye(problem_data.k)*z[problem_data.m]; # - minimum([0,di]) + 0.0001 # #spdiagm( -0.25*(variables.x + 1).^(-1.5) ) #spzeros(length(variables.x),length(variables.x))#
		#-0.01*2*ones(problem_data.k)*x' + 
		#return spzeros(length(x),length(x))
	end
	
	problem_data.add_black_box_constraints(evaluate_inequality_constraint, evaluate_inequality_constraint_gradients, evaluate_inequality_constraint_lagrangian)
end

######################
#  run the program   #
######################
main()
