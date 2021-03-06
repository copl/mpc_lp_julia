using JuMP

include("../src/homogeneous_algorithm/ip_algorithm.jl")

EMPTY_ARRAY = spzeros(0,1)*[0.0]


function is_problem_successful(problem_name,status,correct_value,details=None)
	if status == correct_value
		println(problem_name, " sucessfully solved")
	else
		println("---------------------------------------")
		println(problem_name, " NOT SOLVED SUCESSFULLY.")
		println("Status = ", status, " instead of ", correct_value)
		
		if details != None
			println("Details:")
			println(details) # i.e. problem description
		end	
		println("---------------------------------------")
	end
	
	return (status == correct_value)
end

function standard_settings()
	settings = class_settings();
	
	settings.max_iter = 70;  # Total number of iterarions
	settings.min_alpha = 1e-4;
	
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
	settings.beta6 = 0.75
	
	settings.diagonal_modification = 1e-8
	
	return settings
end

function get_netlib_problem(file_name)
	file = matopen(file_name)
	A= sparse(read(file, "A"))
	b = read(file,"b")
	b = b[:]
	c = read(file,"c")
	c = c[:];
	
	return A,b,c
end

function solve_with_JuMP(A, b, c, solver=IpoptSolver(max_iter=300)) #Output=0
	model = Model(solver=solver)
	
	n, k = size(A)
	
	@defVar(model, x[1:k] >= 0)
	@setObjective( model, Min, sum{ c[j]*x[j] , j=1:k} )

	for i = 1:n
		@addConstraint(model, sum{A[i,j]*x[j], j=1:k} == b[i])
	end
		
	status = solve(model)
	
	#println("Objective value: ", getObjectiveValue(model))
	if status == :Optimal
		return 1
	elseif status == :Infeasible
		return 2
	elseif status == :Unbounded
		return 3
	else
		return 0
	end
end

function solve_net_lib_problem(A,b,c,settings)
	n,k = size(A)
	
	lp = class_non_linear_program();
	lp.set_linear_objective(c);
	lp.set_linear_constraints(spzeros(0,k),spzeros(0,1)*[1.0],A,b,k);
	
	vars = class_variables(lp);
	
	vars, status, iter = ip_algorithm(lp, settings, vars, false);
	
	return status, iter
end