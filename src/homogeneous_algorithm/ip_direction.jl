println("loading ip_direction")
# TO DO:
#
# - complete line (2) i.e w + alpha*d_w \in F_w
#
# - test MOSEK on infeasible convex problems to check it works
#
# - direction of sufficent descent
#



type class_direction
	# this class is used to store and compute the search direction
	# it handles the linear algebra i.e. solves a KKT system to find the search direction
	
	update_factorization::Function
	factored_K_bar
	line_search::Function
	compute_p_vector::Function
	compute_direction::Function
	
	step_size_line_search::Function
	compute_maximum_step_size::Function
	compute_maximum_step_size_vector::Function
	compute_maximum_step_size_variable::Function
	line_search_merit_function::Function
	compute_alternative_point::Function
	is_valid_step::Function
	
	validate_direction_dimensions::Function
	
	p::Array{Float64,1} 			# solution to system (34)
	
	# variable direction
	dx::Array{Float64,1} 			# positive primal variable
	dx_bar::Array{Float64,1} 		# unbounded primal variable
	dz::Array{Float64,1} 			# primal slack variable
	
	dy::Array{Float64,1} 			# positive dual variable
	dy_bar::Array{Float64,1} 		# slack dual variable
	ds::Array{Float64,1}		 	# dual slack variable
	
	dtau::Float64
	dkappa::Float64
	
	get_point::Function
	
	alpha::Float64
	
	function class_direction()
		this = new();
		
		this.get_point = function()
			return [ this.dx; this.dx_bar; this.dx_bar; this.dy; this.dy_bar; this.ds; this.dz; this.dtau; this.dkappa ];
		end
		
		this.update_factorization = function(local_approximation::class_local_approximation)
			# function creates the K_bar matrix (37)
			# [ D_x	  J'  ]
			# [ J    -D_z ]
			# And factorizes it.
			#
			# return true if factorization corresponds to convex problem
			# otherwise returns false
			
			try
				K_bar = sparse(
						[[local_approximation.D_x local_approximation.J']
						[local_approximation.J -local_approximation.D_z]]
					);
				
								
				this.factored_K_bar = lufact(K_bar); # LDL!!!!
				
			catch e
				println("ERROR class_direction.update_factorization")
				throw(e)
			end
		end
		
		this.compute_p_vector = function(local_approximation,variables)
			try
				this.p = this.factored_K_bar \ [local_approximation.v2; local_approximation.v3];
				this.p[(length(variables.x_scaled())+1):length(this.p)] = -this.p[(length(variables.x_scaled())+1):length(this.p)]
				@assert length(this.p) - length(variables.x_scaled()) == length(variables.y_scaled())
			catch e
				println("ERROR class_direction.compute_p_vector")
				throw(e)
			end	
		end
		
		this.compute_direction = function(local_approximation::class_local_approximation,vars::class_variables)
			try
				gamma = local_approximation.gamma;
				eta = 1.0 - gamma;
				
				# given some eta, gamma this function computes a search direction
				mu = local_approximation.state.mu;
				weights_x = ones(length(vars.x)); # weights for the log barrier
				weights_y = ones(length(vars.y));
				
				# evaluate LHS of (29)
				xs = -(vars.x) .* vars.s + gamma * mu * weights_x;
				zy = -(vars.z) .* vars.y + gamma * mu * weights_y;
				tk = -vars.tau * vars.kappa + gamma * mu;
				
				t_dual = [xs ./ vars.x; zeros(length(vars.x_bar))];
				#t_primal = [zy ./ vars.z; zeros(length(vars.y_bar))];
				t_primal = [zy ./ vars.y; zeros(length(vars.y_bar))];
				
				r_hat = [eta*local_approximation.r_dual + t_dual; eta*local_approximation.r_primal + t_primal]
								
				q = this.factored_K_bar \ r_hat; # (35)
				
				# perform tranformations
				q[(length(vars.x_scaled())+1):length(q)] = -q[(length(vars.x_scaled())+1):length(q)];
				
				v1_v3 = [local_approximation.v1; -local_approximation.v3];
				this.dtau = (eta * local_approximation.r_gap + tk/vars.tau - v1_v3' * q)[1]/(local_approximation.D_g - v1_v3' * this.p)[1];
				direction_vector = q - this.p * (this.dtau); # ERROR IN YINYU ANDERSON 1
				
				k1 = length(vars.x);
				k2 = k1 + length(vars.x_bar);
				k3 = k2 + length(vars.y);
				k4 = k3 + length(vars.y_bar);
				
				this.dx = direction_vector[1:k1];
				this.dx_bar = direction_vector[(k1+1):k2];
				this.dy = direction_vector[(k2+1):k3];
				this.dy_bar = direction_vector[(k3+1):k4];
				
				# (29)
				this.ds = (xs - (vars.s) .* (this.dx))./(vars.x);
				this.dz = (zy - (vars.z) .* (this.dy))./(vars.y);
				this.dkappa = (tk - vars.kappa * this.dtau)/(vars.tau);
				
			catch e
				println("ERROR class_direction.compute_direction")
				throw(e)
			end
		end
		
		this.step_size_line_search = function(vars::class_variables, nlp::class_non_linear_program, local_approx::class_local_approximation, settings::class_settings)
			# Computes the maximum step size
			# Then does a line search on the merit function
			# returns the merit function value and number of evaluations required
			
			try 
				this.compute_maximum_step_size(vars, settings)
				return this.line_search_merit_function(local_approx, nlp, vars, settings)
			catch e
				println("ERROR class_direction.step_size_line_search")
				throw(e)
			end
		end
		
		this.compute_maximum_step_size = function(vars::class_variables, settings::class_settings)
			# compute the largest alpha, such that the bounds are strictly satisfied
			this.alpha = 1/settings.beta5;
			this.compute_maximum_step_size_vector( vars.x, this.dx );
			this.compute_maximum_step_size_vector( vars.y, this.dy );
			this.compute_maximum_step_size_vector( vars.s, this.ds );
			this.compute_maximum_step_size_vector( vars.z, this.dz );
			this.compute_maximum_step_size_variable( vars.tau, this.dtau );
			this.compute_maximum_step_size_variable( vars.kappa, this.dkappa );
			#this.alpha = this.alpha*settings.bkscale;
		end
		
		this.compute_maximum_step_size_vector = function(var::Array{Float64,1}, dvar::Array{Float64,1})
			@assert (length(var) == length(dvar))
			for i = 1:length(var)
				this.compute_maximum_step_size_variable(var[i], dvar[i])
			end
		end
		
		this.compute_maximum_step_size_variable = function(var::Float64, dvar::Float64)
			candidate_alpha = -var/dvar;
			
			if candidate_alpha > 0
				this.alpha = min(this.alpha, candidate_alpha);
			end
		end
				
		this.line_search_merit_function = function(local_approx::class_local_approximation, nlp::class_non_linear_program, vars::class_variables, settings::class_settings)# nlp::class_non_linear_programn)
			try 
				# TO DO:
				# - compute gradient of residuals
				# - do back-tracking line search
				# - make sure x-values match
								
				# compute the alpha that minimizes (or maximizes) the merit function				
				
				previous_merit_function_value = local_approx.update_approximation(nlp,vars,settings);				
				state_vars = deepcopy( local_approx.state );
				#merit_function_derivative = local_approx.calculate_merit_function_derivative( nlp, vars, settings )
				#merit_function_derivative = local_approx.finite_difference_merit_function_derivative( nlp, vars, settings )
				#expected_gain = dot_directions( merit_function_derivative, this );
				
				old_alpha = this.alpha;
				this.alpha = settings.min_alpha;
				temp_var = deepcopy(vars);
				temp_var.take_step(this);
				temp_merit_function = local_approx.update_approximation(nlp,temp_var,settings);
				expected_gain = (previous_merit_function_value - temp_merit_function)/this.alpha
				this.alpha = old_alpha
				
				
				
				if expected_gain <= 0
					println("ERROR Expected gain is negative!!!")
					println("Expected gain =  ", expected_gain)
					println("r-norm ",(state_vars.r_norm - local_approx.state.r_norm)/settings.min_alpha)
					
					@assert(false)
				end
								
				merit_function_value = Inf;
				new_vars = deepcopy(vars);
				this.alpha = this.alpha*settings.beta5;
				
				i = 0;
				max_it = ceil(log(settings.min_alpha)/log(settings.beta4));
				for i = 1:max_it
					w_l = deepcopy(vars);
					w_l.take_step(this);
					X_l = local_approx.update_approximation(nlp,w_l,settings);
					state_w_l = deepcopy(local_approx.state);
					if ~ this.is_valid_step(w_l, state_w_l, state_vars, settings) 
						X_l = Inf
					end
					
					
					w_n = this.compute_alternative_point(w_l,local_approx);
					X_n = local_approx.update_approximation(nlp,w_n,settings);
					state_w_n = deepcopy(local_approx.state);
					if ~ this.is_valid_step(w_n, state_w_n, state_vars, settings) 
						X_n = Inf
					end
					
					 #minimum([X_l,X_n]); #minimum(X_l,X_n)
					
					goal = previous_merit_function_value - this.alpha * settings.beta3 * expected_gain;
					if X_n < goal && X_l < X_n
						new_vars = w_n;
						merit_function_value = X_n
						break;
					elseif X_l < goal
						new_vars = w_l;
						merit_function_value = X_l;
						break;
					end
					
					if this.alpha < settings.min_alpha
						merit_function_value = previous_merit_function_value;
						break
					end
					
					this.alpha = this.alpha*settings.beta4;
					
				end
								
				local_approx.update_approximation(nlp,vars,settings); # restore local approximation
				return merit_function_value, new_vars, i
			catch e
				println("ERROR class_direction.line_search_merit_function")
				throw(e)
			end
		end
		
		this.compute_alternative_point = function(vars::class_variables,local_approx::class_local_approximation)
			# compute the point w_n
			w_n = deepcopy(vars)
			w_n.z = (1 - this.alpha*(1 - local_approx.gamma))*local_approx.r_primal[1:length(vars.z)] + vars.tau*local_approx.a[1:length(vars.z)]
			w_n.s = (1 - this.alpha*(1 - local_approx.gamma))*local_approx.r_dual[1:length(vars.s)] + vars.tau*local_approx.g[1:length(vars.s)]
			
			return(w_n)
			#r_dual_n = [s_n; zeros(length(vars.x_bar))] - vars.tau*local_approx.g
			#r_primal_n = [z_n; zeros(length(vars.y_bar))] - vars.tau*local_approx.a
			
			# residual_change = norm(r_dual_n,1) + norm(r_primal_n,1) - local_approx.r_dual_norm - local_approx.r_primal_norm
			# mu_change = (vars.x)'*(s_n - vars.s) + (vars.y)'*(z_n - vars.z);			
		end
		
		this.is_valid_step = function(vars::class_variables, state_w_new::class_state, state_vars::class_state, settings::class_settings)
			try
				if state_w_new.mu / state_vars.mu < settings.beta1 * state_w_new.r_norm / state_vars.r_norm
					return false
				end
				
				closeness_to_boundary = minimum([minimum([vars.x .* vars.s, [Inf]]),minimum([ vars.z .* vars.y, [Inf]]), vars.tau*vars.kappa])
					
				if closeness_to_boundary < settings.beta2 * state_w_new.mu
					return false;
				end
				
				return true;
			catch e
				println("ERROR class_direction.is_valid_step")
				throw(e)
			end
		end
		
		this.validate_direction_dimensions = function(vars::class_variables)
			try
				@assert(length(vars.x) == length(this.dx))
				@assert(length(vars.x_bar) == length(this.dx_bar))
				@assert(length(vars.y) == length(this.dy))
				@assert(length(vars.y_bar) == length(this.dy_bar))
				@assert(length(vars.s) == length(this.ds))
				@assert(length(vars.z) == length(this.dz))
			catch e
				println("ERROR class_direction.validate_direction_dimensions")
				throw(e)
			end
		end
		
		return(this)
	end
end

function dot_directions(d1::class_direction, d2::class_direction)
	# computes the dot product of two directions
	try
		value = 0
		
		value += (d1.dx)' * d2.dx
		value += (d1.dx_bar)' * d2.dx_bar
		value += (d1.dy)' * d2.dy
		value += (d1.dy_bar)' * d2.dy_bar
		value += (d1.ds)' * d2.ds
		value += (d1.dz)' * d2.dz
		value += (d1.dkappa) * d2.dkappa
		value += (d1.dtau) * d2.dtau
		
		return value[1]
	catch
		println("ERROR dot_directions")
		throw(e)
	end
end