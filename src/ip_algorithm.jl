# 	Problem input is in the form:
#	
# 	c(x,x_bar)
# 	a_i(x,x_bar) >= 0 		for	i = 1, ..., m_1
# 	a_i(x,x_bar) =  0		for	i = 1, ..., m_2
#	x >= 0
#
# 	this solver is based on:
# 	"A Computational Study of the Homogeneous Algorithm for Large-scale Convex Optimization"
#	Computational Optimization and Applications
#  	By Erling D. Anderson, Yinyu Ye, 1998
#


include("ip_core.jl")

println("loading ip_algorithm")

function ip_algorithm(nlp::class_non_linear_program,	settings::class_settings, variables::class_variables)
	try
		status = {"primal_feasible"=>None,"dual_feasible"=>None};
		
		validate(nlp,variables)
	
		println("It | alpha | gamma | tau  | kappa  |  mu  |  gap  | primal | dual | trial#")
		
		
		# allocate memory for working variables	
		debug_message_high_level("allocate memory")
		local_approx = class_local_approximation();
		direction = class_direction();
		

		number_of_merit_function_evaluations = 0;
		direction.alpha = 0.0;
		
		early_termination = false 
		for itr = 0:settings.max_iter
			try				
				orginal_merit_function_value = local_approx.update_approximation(nlp, variables, settings);
				
				print_status(direction, local_approx.state, itr, variables, number_of_merit_function_evaluations, local_approx.gamma) 
				
				debug_message_high_level("check termination criterion")
				
				status = settings.solution_status(local_approx.state);
				
				if status != 0
					early_termination = true
					println("Termination criteron met")
					settings.print_status(status)
					break
				end
				
				direction.update_factorization(local_approx);
				
				direction.compute_p_vector(local_approx,variables);
				
				debug_message_high_level("predictor")
				
				
				# compute predictor
				local_approx.gamma = 0.0;
				direction.compute_direction(local_approx, variables);
				predictor_merit_function_value, temp, number_of_merit_function_evaluations_predictor = direction.step_size_line_search(variables, nlp, local_approx, settings);
				
				# compute corrector
				v = predictor_merit_function_value/orginal_merit_function_value
				local_approx.gamma = v*min(v, 1.5*settings.beta6); # Mehrotra heuristic, see pg 257					
				
				direction.compute_direction(local_approx, variables);
				merit_function_value, variables, number_of_merit_function_evaluations_corrector = direction.step_size_line_search(variables, nlp, local_approx, settings);
				
				#println("alpha",direction.alpha)
				debug_message_high_level("step")
				# step to next point
				number_of_merit_function_evaluations = number_of_merit_function_evaluations_predictor + number_of_merit_function_evaluations_corrector
				#println(number_of_merit_function_evaluations_corrector)
				#variables.take_step(direction);
			catch e
				println("ERROR iteration ", itr)
				throw(e)
			end
		end
		
		if ~early_termination
			println("MAXIMUM ITERATIONS REACHED")
		end
		
		# TO DO RETURN MORE INFORMATION
		return(variables,status)
	catch e
		println("ERROR ip_algorithm")
		throw(e)
	end
end

function evaluate_mu(nlp,variables,direction)
	mu_x = (variables.x + direction.alpha*direction.dx)'*(variables.s + direction.alpha*direction.ds);
	mu_tau = (variables.tau + direction.alpha*direction.dtau)*(variables.kappa + direction.alpha*direction.dkappa);
	mu_y = (variables.y + direction.alpha*direction.dy)'*(variables.z + direction.alpha*direction.dz);
	return (mu_x + mu_tau + mu_y)[1]/(nlp.n_1 + nlp.m_1 + 1)
end

function print_status(direction, state, itr, variables, number_of_merit_function_evaluations, gamma) 
	# %%%%%%%%%%%%
	# %%%%%%%%%%%%
	# %%%%%%%%%%%%
	#if local_approx.r_gap >= 0
	#	gap_sign_string = " "
	#else
	#	gap_sign_string = "-"
	#end
	
	@printf("%s %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %i\n", rpad(string(itr),3), direction.alpha, gamma, variables.tau, variables.kappa, state.mu, state.r_gap_norm, state.r_primal_norm, state.r_dual_norm, number_of_merit_function_evaluations)
end

function evaluate_solution_status(state::class_state, settings::class_settings)
	try 
		terminate = (state.mu < settings.primal_feas_tol) && (state.r_gap_norm < settings.primal_feas_tol) && (state.r_primal_norm < settings.primal_feas_tol) && (state.r_dual_norm < settings.primal_feas_tol)
		
		return  terminate, {"primal_feasible"=>None,"dual_feasible"=>None}
	catch e
		println("ERROR evaluate_solution_status")
		throw (e)
	end
end