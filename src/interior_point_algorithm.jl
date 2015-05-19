include("interior_point_data_structures.jl")

function interior_point_algorithm(problem_data::class_linear_program_input,	settings::class_settings, variables)
	try
		problem_data.validate()
		problem_data.update(variables);
		
		println("Iter\tgap \t      mu\talpha\ttau\tfeasibility")
		
		state = class_algorithm_state()
		
		# allocate memory for working variables	
		debug_message_high_level("allocate memory")
		K_newton_matrix = class_K_newton_matrix(problem_data,settings); # Data structure used to share info between two Linear Eq. Solver in affine direction and corrector step.
		rhs = class_linear_system_rhs(problem_data);
		direction = class_direction(problem_data); # Stores affine and corrector directions
		residuals = class_residuals(problem_data);
		
		early_termination = false 
		for itr = 1:settings.max_iter
			# compute the gap between primal and dual after the step
			problem_data.update(variables);
			state.update_mu(variables,problem_data)
			state.update_gap(variables,problem_data)
			
			# santiago's solver needs this
			debug_message_high_level("K newton matrix")
			K_newton_matrix.update(problem_data,variables);
			
			# compute the residuals
			debug_message_high_level("compute_residuals")
			residuals.compute_residuals(problem_data,variables);
			
			print_status(state, direction, variables, residuals, itr)
			
			debug_message_high_level("check termination criterion")
			if termination_criterion_met(settings,state,residuals)
				early_termination = true
				println("Termination criteron met")
				break
			end
			
			# compute rhs for affine direction
			debug_message_high_level("affine rhs")
			rhs.compute_affine_rhs(residuals,variables)
			# compute affine direction using new affine rhs
			debug_message_high_level("affine direction")
			direction.compute_affine_direction(rhs,problem_data,variables,state,settings,K_newton_matrix); # TO DO - FIX ORDER
			
			debug_message_high_level("corrector rhs")
			# update corrector rhs using new affine direction
			rhs.compute_corrector_rhs(residuals,variables,state,direction,problem_data)
			# update corrector direction using new corrector rhs
			debug_message_high_level("corrector direction")
			direction.compute_corrector_direction(rhs,problem_data,variables,state,settings,K_newton_matrix);
			
			debug_message_high_level("step")
			# take step in the corrector direction
			variables.take_step(direction)
			variables.check()			
		end
		
		if ~early_termination
			println("MAXIMUM ITERATIONS REACHED")
			status = -1;
		else
			status = evaluate_solution_status(variables, residuals, problem_data, settings)
		end
		
		# TO DO RETURN MORE INFORMATION
		return(variables,status)
	catch e
		println("*** ERROR in interior point algorithm ***")
		throw(e)
	end
end

function print_status(state, direction, variables, residuals, itr) 
	@printf("%i\t%2.1e\t%2.1e\t%2.1e\t%2.1e\t%2.1e\n", itr, state.gap, state.mu, direction.alpha, variables.tau, residuals.summary_norm)
end

function evaluate_solution_status(variables, residuals, problem_data, settings)
	# TO DO
	# update so it cannot make mistakes i.e. termination implies a descision
	
	# update criteron to more accurately decide termination
	# println(residuals.primal_feasibility)
	# println(residuals.dual_feasibility)
	
	status = [true,true]; 
	if (variables.tau/variables.kappa)*norm(problem_data.c,1) < settings.comp_tol*100
		if residuals.primal_feasibility < -settings.comp_tol
			status[1] = false;
		end
		if residuals.dual_feasibility < -settings.comp_tol
			status[2] = false;
		end
	elseif variables.kappa/variables.tau > settings.comp_tol
		println(variables.kappa/variables.tau)
		println("could not conclude if problem was feasible or infeasible. perhaps modify stopping conditions")
		Error("could not conclude if problem was feasible or infeasible. perhaps modify stopping conditions")
	end
	
	return status
end

function termination_criterion_met(settings::class_settings,state::class_algorithm_state,residuals::class_residuals)
	# TO DO
	# store a bunch of norms
	
	#Evaluate termination criteria	
	#TODO: This part will have to be a more sophisticated test to detect 
	#unbounded and infeasible problems.
	
	if (
		residuals.r1_norm < settings.linear_feas_tol && residuals.r2_norm < settings.linear_feas_tol && residuals.r3_norm < settings.linear_feas_tol
		&& residuals.r4_norm < settings.linear_feas_tol 
		||
		(residuals.kkt_is_infeasible < settings.comp_tol && state.mu < settings.comp_tol)
		)
		 return true
	else
		return false
	end
end
