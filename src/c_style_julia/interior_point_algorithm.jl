include("interior_point_data_structures.jl")
include("linear_system_solver.jl")

function interior_point_algorithm(problem_data::class_linear_program_input,	settings::class_settings)

	println("Iter\tgap \t      mu\talpha\ttau\tfeasibility")

	variables = class_linear_program_variables(problem_data) 
	state = class_algorithm_state()
	state.update_mu(variables,problem_data)
	
	
	# allocate memory for working variables	
	K_newton_matrix = class_K_newton_matrix(problem_data); # Data structure used to share info between two Linear Eq. Solver in affine direction and corrector step.
	rhs = class_linear_system_rhs(problem_data);
	direction = class_direction(problem_data); # Stores affine and corrector directions
	residuals = class_residuals(problem_data);
	
	early_termination = false 
	for itr = 1:settings.max_iter
		# santiago's solver needs this
		#println("Time to construct K-newton matrix")
		K_newton_matrix.update(variables);
		
		# compute the residuals
		residuals.compute_residuals(problem_data,variables);
		
		if termination_criterion_met(settings,state,residuals)
			early_termination = true
			print_status(state, direction, variables, residuals, itr)
			println("Termination criteron met")
			break
		end
		
		# compute rhs for affine direction
		rhs.compute_affine_rhs(residuals,variables)
		# compute affine direction using new affine rhs
		direction.compute_affine_direction(rhs,problem_data,variables,state,settings,K_newton_matrix); # TO DO - FIX ORDER
		
		# update corrector rhs using new affine direction
		rhs.compute_corrector_rhs(residuals,variables,state,direction,problem_data)
		# update corrector direction using new corrector rhs
	    direction.compute_corrector_direction(rhs,problem_data,variables,state,settings,K_newton_matrix);
		
		# take step in the corrector direction
		variables.take_step(direction)
	    
		# compute the gap between primal and dual after the step
		state.update_mu(variables,problem_data)
		state.update_gap(variables,problem_data)
		
		print_status(state, direction, variables, residuals, itr)
	end
	
	if ~early_termination
		println("MAXIMUM ITERATIONS REACHED")		
	end
	
	# TO DO RETURN MORE INFORMATION
	return(variables)
end

function print_status(state, direction, variables, residuals, itr)
	@printf("%i\t%2.1e\t%2.1e\t%2.1e\t%2.1e\t%2.1e\n", itr, state.gap, state.mu, direction.alpha, variables.tau, residuals.residuals_norm)
end

function termination_criterion_met(settings::class_settings,state::class_algorithm_state,residuals::class_residuals)
	# TO DO
	# store a bunch of norms
	
	#Evaluate termination criteria	
	#TODO: This part will have to be a more sophisticated test to detect 
	#unbounded and infeasible problems.
	
	if (residuals.residuals_norm < settings.linear_feas_tol && 
		state.mu < settings.comp_tol)
		 return true
	else
		return false
	end
end
