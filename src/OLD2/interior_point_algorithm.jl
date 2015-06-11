include("interior_point_data_structures.jl")

function interior_point_algorithm(problem_data::class_linear_program_input,	settings::class_settings, variables::class_linear_program_variables)
	try
		status = {"primal_feasible"=>None,"dual_feasible"=>None};
		
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
			
			# santiago"s solver needs this
			debug_message_high_level("K newton matrix")
			K_newton_matrix.update(problem_data,variables);
			
			# compute the residuals
			debug_message_high_level("compute_residuals")
			residuals.compute_residuals(problem_data,variables);
			
			print_status(state, direction, variables, residuals, itr)
			
			debug_message_high_level("check termination criterion")
			
			terminate, status = evaluate_solution_status(variables, residuals, settings);
						
			if terminate # termination_criterion_met(settings,state,residuals)
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
			#variables.line_search(direction,residuals,problem_data)
			
			debug_message_high_level("corrector rhs")
			# update corrector rhs using new affine direction
			rhs.compute_corrector_rhs(residuals,variables,state,direction,problem_data)
			# update corrector direction using new corrector rhs
			debug_message_high_level("corrector direction")
			direction.compute_corrector_direction(rhs,problem_data,variables,state,settings,K_newton_matrix);
			
			debug_message_high_level("step")
			# take step in the corrector direction
			#variables.line_search(direction,residuals,problem_data)
			variables.take_step(direction)
			variables.check()
		end
		
		if ~early_termination
			println(variables.tau*norm(problem_data.grad,1))
			println("MAXIMUM ITERATIONS REACHED")
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


function evaluate_solution_status(variables::class_linear_program_variables, residuals::class_residuals, settings::class_settings)
	try
		# TO DO
		# update so it cannot make mistakes i.e. termination implies a descision
		
		# update criteron to more accurately decide termination

		status = {"primal_feasible"=>None,"dual_feasible"=>None}; # primal 
		if residuals.primal_feasibility < settings.primal_feas_tol
			status["primal_feasible"] = true
		elseif 0 < residuals.primal_infeasibility && residuals.primal_infeasibility < settings.primal_infeas_tol
			status["primal_feasible"] = false
		end
		
		if residuals.dual_feasibility < settings.dual_feas_tol
			status["dual_feasible"] = true
		elseif 0 < residuals.dual_infeasibility && residuals.dual_infeasibility < settings.dual_infeas_tol
			status["dual_feasible"] = false
		end
				
		if status["primal_feasible"] == true && status["dual_feasible"] == true
			if residuals.duality_gap < settings.duality_gap_tol
				return true, status
			else
				return false, status
			end
		elseif status["primal_feasible"] == false || status["dual_feasible"] == false
			return true, status
		else
			return false, status
		end
	catch e
		println("*** ERROR in evaluate_solution_status ***")
		throw(e)
	end
end
