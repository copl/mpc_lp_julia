include("interior_point_data_structures.jl")
include("linear_system_solver.jl")


function interior_point_algorithm(problem_data::class_linear_program_input,	settings::class_settings)

	println("Iter \t gap     \t mu   \t alpha     \t tau       \t feasibility")

	variables = class_linear_program_variables(problem_data) 
	state = class_algorithm_state()
	state.update_mu(variables,problem_data)
	
	K_newton_matrix = class_K_newton_matrix(problem_data);
	rhs = class_linear_system_rhs(problem_data);
	direction = class_direction(problem_data);
	residuals = class_residuals(problem_data);
		
	for itr =1:settings.max_iter
		K_newton_matrix.update(variables);
		
		# compute the residuals
		residuals.compute_residuals(problem_data,variables);
		
		if termination_criterion_met(settings,state,residuals)
			break
		end
		
		# compute affine rhs
		rhs.compute_affine_rhs(residuals,variables)
		# compute affine direction using new affine rhs
		direction.compute_affine_direction(rhs,problem_data,variables,K_newton_matrix); # TO DO - FIX ORDER
		
		# update corrector rhs using new affine direction
		rhs.compute_corrector_rhs(residuals,variables,state,direction,problem_data)
		# update corrector direction using new corrector rhs
	    direction.compute_corrector_direction(direction,rhs,problem_data,
											  variables,
											  state,
											  settings,
											  K_newton_matrix);
		
		
		# take step in the corrector direction
		variables.take_step(direction)
	    
		# compute the gap after the step
		state.update_mu(variables,problem_data)
		state.update_gap(variables,problem_data)
		
		@printf("%3i\t%3.3e\t%3.3e\t%3.3e\t%3.3e\t%3.3e\n", itr, state.gap ,state.mu, direction.alpha, variables.tau, residuals.normed_squared)
	end
	
	# TO DO RETURN MORE INFORMATION
	return(variables)
end


function termination_criterion_met(settings::class_settings,state::class_algorithm_state,residuals::class_residuals)
	# TO DO
	# store a bunch of norms
	
	#Evaluate termination criteria	
	#TODO: This part will have to be a more sofisticated test to detect 
	#unbounded and infeasible problems.
	
	if (residuals.r1_norm < settings.linear_feas_tol && 
		residuals.r2_norm < settings.linear_feas_tol &&
		residuals.r3_norm < settings.linear_feas_tol && 
		state.mu < settings.comp_tol)
		 println("Ended");
		 return true
	else
		return false
	end
end
