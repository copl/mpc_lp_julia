println("loading libraries")

function debug_message(message)
	#println(message)
end
include("interior_point_algorithm.jl")

println("internal libraries loaded")

using MAT
using JuMP
using Gurobi
using Ipopt

println("external libraries loaded")


function main() # this function is called at the bottom of the code
	
	max_iter             = 50;  # Total number of iterarions
    linear_feas_tol      = 1e-6;  # Threshold for feasibility of the linear constraints.
    comp_tol            = 1e-6;  # Threshold for complementry optimality condition s^Tz
    bkscale              = 0.98;  # Back scaling for line search.

   # Initialize configuration variable
    settings = class_settings(max_iter,linear_feas_tol,comp_tol,bkscale);

	srand(1234)
	# We are creating an instance of LP. In practice, we should read the problem data from input stream.
	problem_data = construct_instance_from_mat("Problems/QAP8.mat");
	
	xTrue = solve_with_JuMP(problem_data)
	
	# The main function that run interior point algorithm.
	variables = @time interior_point_algorithm(problem_data,settings);
	
	# do we have an optimal solution ?
	interpret_variables(variables)
	
	println(variables.tau)
	println(variables.kappa)
	
	xOur = variables.x/variables.tau;
	println("Our objective: " * string(get_objective_value(problem_data,xOur)))
	println("JuMP objective: " * string(get_objective_value(problem_data,xTrue)))
	
	println("Our feasiblity: " * string(calculate_feasibility(problem_data,xOur)))
	println("JuMP feasibility: " * string(calculate_feasibility(problem_data,xTrue)))
	
	
	#println(variables.x/variables.tau)
	
	println("Distance to JuMP solver solution: "* string(norm(xTrue - xOur)))
end

function get_objective_value(problem_data,x)
	return(problem_data.c' * x + 0.5*x'*problem_data.P*x)
end

function calculate_feasibility(problem_data,x)
	return norm([problem_data.A*x - problem_data.b, max(problem_data.G*x - problem_data.h,0)])
end

function interpret_variables(variables)
	if variables.tau > variables.kappa
		println("Optimal solution found")
	else:
		println("Problem is infeasible or unbounded")
	end
end

function solve_with_JuMP(problem_data::class_linear_program_input)
	model = Model(solver=GurobiSolver())
	#model = Model(solver=IpoptSolver())
	k = problem_data.k
	m = problem_data.m
	n = problem_data.n
	
	A = problem_data.A
	G = problem_data.G
	c = problem_data.c
	b = problem_data.b
	h = problem_data.h
	
	@defVar(model, x[1:k] )
	@setObjective(model, Min, sum{c[j]*x[j] + 20*0.5*x[j]*x[j], j=1:k} )
	
	#Afull = full(A)
	#indicies = M.rowval[a.colptr[col] : M.colptr[col+1]-1] 

	for i = 1:n
		@addConstraint(model, sum{A[i,j]*x[j], j=1:k} == b[i])
	end
	
	for i = 1:m
		@addConstraint(model, sum{G[i,j]*x[j], j=1:k} <= h[i])
	end
	
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

function construct_instance_from_mat(filename)
	println("Reading Matlab file " * filename )

	file = matopen(filename)
	A = sparse(read(file, "A"))
	n,k = size(A)
	m = k
	G = -speye(m)
	
	c = read(file,"c")
	b = read(file,"b")
	h = zeros(m,1)
	
	close(file)
	
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
	problem_data.P = 20*speye(k,k)
	problem_data.A = sparse(A)
	problem_data.G = sparse(G)
	problem_data.c = c[:]
	problem_data.h = h[:]
	problem_data.b = b[:]
	
	problem_data.m = m
	problem_data.k = k
	problem_data.n = n
	
	return(problem_data)
end

######################
#  run the program   #
######################
main()
