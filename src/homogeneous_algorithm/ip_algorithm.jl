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

function ip_algorithm(nlp::class_non_linear_program,	settings::class_settings, variables::class_variables, print_output=true)
	try
		validate(nlp,variables)
	
		if print_output
			println("It | alpha | gamma | tau  | kappa  |  mu  |  gap  | primal | dual | trial#")
		end
		
		# allocate memory for working variables	
		local_approx = class_local_approximation();
		direction = class_direction();
		

		number_of_merit_function_evaluations = 0;
		direction.alpha = 0.0;
		
		delta = 1e-4

		for itr = 0:settings.max_iter
			try	

				local_approx.sigma = 0.0;
				orginal_merit_function_value = local_approx.update_approximation(nlp, variables, settings);

				if print_output
					print_status(direction, local_approx.state, itr, variables, number_of_merit_function_evaluations, local_approx.gamma)
				end
								
				status = settings.solution_status(local_approx.state);
				
				# terminate if we have solved the problem or the maximum iterations are reached
				if status != 0	
				
					if print_output
						println("Termination criteron met")
						settings.print_status(status)
					end
					
					return(variables, status, itr)
					
				elseif itr == settings.max_iter
					if print_output
						println("MAXIMUM ITERATIONS REACHED")
					end
					
					return(variables, status,itr)
					
				end
				

				#direction.update_factorization(local_approx,nlp,variables,settings)
				#sigma = direction.factored_K_bar.calc_min_mod()

				
				#println("Omega = ",omega)				
				
				for j = 1:100
					nlp_convex = nlp.convexify(delta, variables.x_scaled()) ##################################
					local_approx.update_approximation(nlp_convex, variables, settings);
					
					if direction.update_factorization(local_approx,nlp_convex,variables,settings)
						#println("is inertia correct: ",direction.factored_K_bar.is_convex)

						direction.compute_p_vector(local_approx,variables);				
				
						# compute predictor
						predictor_merit_function_value = orginal_merit_function_value;
						number_of_merit_function_evaluations_predictor = 0;

						orginal_mu = local_approx.state.mu
						try
							local_approx.gamma = 0.0;
							direction.compute_direction(local_approx, variables);
							direction.compute_maximum_step_size(variables, settings)
							#@assert(direction.alpha > settings.min_alpha)
							predictor_merit_function_value, temp, number_of_merit_function_evaluations_predictor = direction.step_size_line_search(variables, nlp_convex, local_approx, settings);
						
						
						catch e
							#println(full(direction.K_bar))
							#println("x-scaled",variables.x_scaled())	
							#println("x/s",(variables.s)./variables.x)
							#println("x",(variables.x))
							#println("dx",direction.dx)					
							println("ERROR predictor")
							println(e)
							println("===============")
							#throw(e)
						end

						# compute corrector
						number_of_merit_function_evaluations_corrector = 0;
						step_alpha = NaN;
						try
							#new_mu = local_approx.state.mu
							#v = new_mu/orginal_mu
							v = predictor_merit_function_value/orginal_merit_function_value;
						
							local_approx.gamma = min(v,1.0)*min(v, settings.beta6); # Mehrotra heuristic, see pg 257	

							direction.compute_direction(local_approx, variables);
							direction.compute_maximum_step_size(variables, settings)
							#step_alpha = direction.alpha
							#@assert(direction.alpha > 10*settings.min_alpha)
										
							#variables.take_step(direction)
							#direction.alpha = 0.9*direction.alpha						
							merit_function_value, variables, number_of_merit_function_evaluations_corrector = direction.step_size_line_search(variables, nlp_convex, local_approx, settings);
							#variables.take_step(direction);
							#@assert(direction.alpha > settings.min_alpha)				
						catch e
							println("ERROR corrector")
							println("step Alpha = ", step_alpha)							
							println("Alpha = ", direction.alpha)						
							println("Gamma = ", local_approx.gamma)
							throw(e)
						end
				
				
					 
						# step to next point
						number_of_merit_function_evaluations = number_of_merit_function_evaluations_predictor + number_of_merit_function_evaluations_corrector

						delta = delta/(1.3) #max(delta/((direction.alpha/0.3)),1e-10)
						break
					else
						delta = delta*11.0;
					end
				end

				println("delta =", delta)
				@assert(delta < 1e10)
			catch e
				println("ERROR iteration ", itr + 1)
				#println(variables.x_scaled(),variables.y_scaled())
				throw(e)
			end
		end
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
	scale = variables.tau + variables.kappa
	@printf("%s %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %2.1e %i\n", rpad(string(itr),3), direction.alpha, gamma, variables.tau, variables.kappa, state.mu/scale, state.r_gap_norm/scale, state.r_primal_norm/scale, state.r_dual_norm/scale, number_of_merit_function_evaluations)

	println("merit = ", (state.mu/sqrt(length(variables.x) + length(variables.y) + 1) + state.r_norm)/scale)
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
