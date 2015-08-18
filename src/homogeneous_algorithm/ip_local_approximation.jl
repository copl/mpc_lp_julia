println("loading ip_local_approximation")
# at each stage we create a quadratic approximation of the current problem we
# use this to take the netwon step and improve the performance

function diagonally_dominant(H)
	# finds the smallest diagonal that needs to be added to make the matrix diagonally dominant.
	# will replace in the future with a modified cholesky factorization
	
	n = size(H,1)
	@assert(size(H) == (n,n))
	#@assert()
	
	vec = zeros(n)
	for i = 1:n
		row_sum = 0;
		for j = 1:n
			if i != j
				row_sum += H[i,j];
			end
		end
		vec[i] = max(row_sum - H[i,i], 0);
	end
	
	return vec
end


type class_state
	r_dual_norm::Float64
	r_primal_norm::Float64
	r_gap_norm::Float64
	r_norm::Float64
	mu::Float64

	kkt_feasible::Float64
	
	relative_gap::Float64
	dual_feasibility::Float64
	primal_feasibility::Float64
	
	
	homogeneous_primal_sign::Int64
	homogeneous_dual_sign::Int64
	homogeneous_primal_feasibility::Float64
	homogeneous_dual_feasibility::Float64
	
	update_state::Function
	
	function class_state()
		this = new();
		
		this.update_state = function(local_approx::class_local_approximation, vars::class_variables, nlp::class_non_linear_program)
			try 
				# update class_state
				this.r_dual_norm = norm(local_approx.r_dual, GLOBAL_P_NORM)
				this.r_primal_norm = norm(local_approx.r_primal, GLOBAL_P_NORM)
				this.r_gap_norm = abs(local_approx.r_gap)
				this.r_norm = norm([ this.r_dual_norm, this.r_primal_norm, this.r_gap_norm ], GLOBAL_P_NORM)
				this.mu = (vars.x'*vars.s + vars.z'*vars.y + vars.tau*vars.kappa)[1]/(nlp.n_1 + nlp.m_1 + 1);
				
				this.kkt_feasible = vars.tau/vars.kappa
				
				# relative stopping criteron we may want to change this match erling and yinyu.
				this.relative_gap = this.mu/(vars.tau + vars.kappa) # / norm( [vars.x; vars.x_bar; vars.y; vars.y_bar; vars.s; vars.s; vars.tau; vars.kappa], GLOBAL_P_NORM );
				this.dual_feasibility = this.r_dual_norm/vars.tau #/  norm( [vars.y; vars.y_bar; vars.s ], GLOBAL_P_NORM );
				this.primal_feasibility = norm(local_approx.a, GLOBAL_P_NORM) #/ norm( [vars.x; vars.x_bar; vars.z ], GLOBAL_P_NORM );

				homogeneous_dual_objective = -(local_approx.v3' * [vars.y; vars.y_bar])[1];
				homogeneous_primal_objective = -(local_approx.c' * [vars.x; vars.x_bar])[1];
				
				this.homogeneous_primal_sign = sign(homogeneous_primal_objective)
				this.homogeneous_dual_sign = sign(homogeneous_dual_objective)
				
				this.homogeneous_primal_feasibility = norm(local_approx.J * [vars.x; vars.x_bar] - [vars.z; zeros(length(vars.y_bar))], GLOBAL_P_NORM ) / abs(homogeneous_primal_objective);
				
				this.homogeneous_dual_feasibility = norm((local_approx.J)' * [vars.y; vars.y_bar] + [vars.s; zeros(
				length(vars.x_bar))], GLOBAL_P_NORM ) / abs(homogeneous_dual_objective);			
				#println(this.primal_infeasibility)
			catch e
				println("ERROR class_state.update_state")
				throw (e)
			end
		end
		
		return this
	end
end

type class_local_approximation
	# makes a local approximation of constraints and objective at current point
	
	# deal with non-convexity
	convexify_hessian::Function
	reset_hessian_modifications::Function
	
	update_approximation::Function # update the local approximation at the new point

	theta::Function	
	calculate_merit_function::Function
	calculate_merit_function_expected_progress::Function
	calculate_merit_function_expected_progress_finite_difference::Function

	calculate_merit_function_derivative::Function
	calculate_residual_norm_derivative::Function
	finite_difference_merit_function_derivative::Function
	finite_difference_merit_function_derivative_vector::Function
	finite_difference_merit_function_derivative_point::Function
	
	calculate_potential_merit_function::Function
	calculate_potential_merit_function_derivative::Function
	
	gamma::Float64 # predictor-corrector parameter
	sigma::Float64 # parameters that ensures hessian is positive-definite
	


	c::Array{Float64,1}
	#double_c::SparseMatrixCSC{Float64,Int64} # objective hessian
	g::Array{Float64,1}
	J::SparseMatrixCSC{Float64,Int64}
	H::SparseMatrixCSC{Float64,Int64}
	convex_H::SparseMatrixCSC{Float64,Int64}
	a::Array{Float64,1}
	
	current_objective_value::Float64
	
	v1::Array{Float64,1}
	v2::Array{Float64,1}
	v3::Array{Float64,1}
	
	D_g::Float64
	D_x::SparseMatrixCSC{Float64,Int64}
	D_z::SparseMatrixCSC{Float64,Int64}
	
	r_dual::Array{Float64,1}
	r_primal::Array{Float64,1}
	r_gap::Float64
	
	state::class_state
	
	function class_local_approximation()
		this = new();
		this.gamma = 1.0;
		this.state = class_state();
		
		
		this.update_approximation = function(nlp::class_non_linear_program,vars::class_variables,settings::class_settings)
			try
				this.current_objective_value = nlp.objective_function(vars.x_scaled());
				
				this.c = nlp.objective_function_gradient(vars.x_scaled());
				#this.double_c = nlp.objective_function_hessian(vars.x_scaled())
				this.J = nlp.evaluate_constraint_gradients(vars.x_scaled());
				this.g = this.c - this.J'*vars.y_scaled();
				#this.H = this.double_c - nlp.evaluate_constraint_lagrangian_hessian(vars.x_scaled(),vars.y_scaled())
				this.H = nlp.objective_function_hessian(vars.x_scaled()) - nlp.evaluate_constraint_lagrangian_hessian(vars.x_scaled(),vars.y_scaled())
				#this.update_lambda()	
				
				(n,m) = size(this.H)
				@assert(n == m)			
				

				#this.convex_H = this.H #+ this.sigma * speye(size(this.H,1));
				#this.H = this.convex_H
				#this.H = this.convex_H
				#println(full(this.H))
				#Delta_H = this.convex_H * (vars.x_scaled());
				# - mod * [vars.x; vars.x_bar]

				# In (30)
				this.a  = nlp.evaluate_constraints(vars.x_scaled());
				this.v1 = -this.c - this.H * (vars.x_scaled()); # ERROR IN YINYU ANDERSON 1
				this.v2 = this.c - this.H * (vars.x_scaled());# + 2*mod*vars.x_scaled(); # ERROR IN YINYU ANDERSON 1
				this.v3 = this.a - this.J * (vars.x_scaled());

				# following (31)
				# 1e-8*ones are added incase the constraints are not linearly independent
				this.D_g = (vars.x_scaled()' * this.H * vars.x_scaled())[1] + vars.kappa/vars.tau;

				#this.D_g = max(this.D_g,10)

				this.D_x = this.H + spdiagm([(vars.s)./(vars.x); zeros(length(vars.x_bar))] + settings.diagonal_modification*ones(length(vars.x_scaled())));
				this.D_z = spdiagm(  [(vars.z)./(vars.y); zeros(length(vars.y_bar))]  + settings.diagonal_modification*ones(length(vars.y_scaled())) );
				

				# (27) the residuals of our problem
				this.r_dual = ([vars.s; zeros(length(vars.x_bar))] - vars.tau*(this.g));
				this.r_gap = (vars.kappa + vars.tau*((vars.x_scaled()'*this.g) + (vars.y_scaled()'*this.a)))[1];
				this.r_primal = ([vars.z; zeros(length(vars.y_bar))] - vars.tau*this.a);

				this.state.update_state(this,vars,nlp)

				return this.calculate_merit_function(nlp,vars,settings);
			catch e
				println("ERROR in class_local_approximation.update_approximation")
				throw(e)
			end
		end
		
		this.theta = function(nlp::class_non_linear_program)
			return sqrt(nlp.m_1 + nlp.n_1 + 1.0)
		end

		this.calculate_merit_function = function(nlp::class_non_linear_program,vars::class_variables,settings::class_settings)
			try
				return this.theta(nlp)*(this.state.mu) + this.state.r_norm;
			catch e
				println("ERROR class_local_approximation.calculate_merit_function")
				throw(e)
			end
		end

		this.calculate_merit_function_expected_progress = function(nlp::class_non_linear_program)
			return ( 1.0 - this.gamma ) * ( this.theta(nlp) * this.state.mu + this.state.r_norm )
		end
		
		this.calculate_merit_function_expected_progress_finite_difference = function(nlp::class_non_linear_program,vars,dir)
			old_alpha = dir.alpha;
			dir.alpha = settings.min_alpha;
			temp_var = deepcopy(vars);
			temp_var.take_step(dir);
			temp_merit_function = local_approx.update_approximation(nlp,temp_var,settings);
			expected_gain_est = (previous_merit_function_value - temp_merit_function)/this.alpha
			dir.alpha = old_alpha

			return expected_gain_est
		end

		this.calculate_merit_function_derivative = function(nlp::class_non_linear_program,vars::class_variables,settings::class_settings)
			try
				# compute the derivative of the residual norm
				merit_function_derivative = this.calculate_residual_norm_derivative(nlp,vars,settings);
				
				# add the derivative of duality gap theta*g(w)
				theta = 1.0 / sqrt(nlp.m_1 + nlp.n_1 + 1);
				merit_function_derivative.dx += theta * vars.s;
				merit_function_derivative.ds += theta * vars.x;
				
				merit_function_derivative.dy += theta * vars.z;
				merit_function_derivative.dz += theta * vars.y;
				
				merit_function_derivative.dtau += theta * vars.kappa;
				merit_function_derivative.dkappa += theta * vars.tau;
				
				return merit_function_derivative
			catch e
				println("ERROR class_local_approximation.calculate_merit_function_derivative")
				throw(e)
			end
		end
		
		this.calculate_residual_norm_derivative = function(nlp::class_non_linear_program,vars::class_variables,settings::class_settings)
			try
				# 
				sigma_D = GLOBAL_P_NORM * (this.r_dual) .^ (GLOBAL_P_NORM - 1);
				sigma_G = GLOBAL_P_NORM * (this.r_gap) ^ (GLOBAL_P_NORM - 1);
				sigma_P = GLOBAL_P_NORM * (this.r_primal) .^ (GLOBAL_P_NORM - 1);
				
				residual_norm_derivative = class_direction();
				
				#dx = -this.H * sigma_D + sigma_G * ( this.g + (this.J)' * vars.y_scaled() ) - (this.J)' * sigma_P; 
				dx = this.H * sigma_D - sigma_G * this.v1 + (this.J)' * sigma_P; 
				#dy = this.J * sigma_D + sigma_G * ( -this.J * vars.x_scaled() + this.a ); 
				dy = -this.J * sigma_D + this.v3 * sigma_G; 
				
				residual_norm_derivative.dkappa = sigma_G;
				
				# tau (the difficult one)
				#dg_dtau = (-1/vars.tau) * (this.H * vars.x_scaled() - this.J' * vars.y_scaled());
				
				#drD_dtau = -this.g - vars.tau * dg_dtau
				#drG_dtau = vars.x_scaled()' * dg_dtau - vars.y_scaled()' * this.J * vars.x_scaled();
				#drP_dtau = this.a - this.J * vars.x_scaled();
				#residual_norm_derivative.dtau = (drD_dtau' * sigma_D + drG_dtau' * sigma_G + drP_dtau' * sigma_P)[1];
				H_g = (vars.x_scaled()' * this.H * vars.x_scaled())[1];
				residual_norm_derivative.dtau = (this.v2' * sigma_D + H_g * sigma_G + this.v3' * sigma_P)[1];
				residual_norm_derivative.ds = -sigma_D[1:length(vars.x)];
				residual_norm_derivative.dz = -sigma_P[1:length(vars.y)];
				
				residual_norm_derivative.dx = dx[1:length(vars.x)];
				residual_norm_derivative.dx_bar = dx[(length(vars.x)+1):length(vars.x_scaled())];
				residual_norm_derivative.dy = dy[1:length(vars.y)];
				residual_norm_derivative.dy_bar = dy[(length(vars.y)+1):length(vars.y_scaled())];
				
				residual_norm_derivative.validate_direction_dimensions(vars)
				
				return(residual_norm_derivative)
			catch e
				println("ERROR class_local_approximation.calculate_norm_derivative")
				throw(e)
			end
		end	
		
		this.finite_difference_merit_function_derivative = function(nlp::class_non_linear_program, vars::class_variables, settings::class_settings)
			try
				merit_derivative = class_direction();
				orginal_value = this.calculate_merit_function(nlp,vars,settings);
				
				h = 1e-6
				merit_derivative.dx = this.finite_difference_merit_function_derivative_vector(h, nlp, vars, settings, orginal_value, vars.x)
				merit_derivative.dx_bar = this.finite_difference_merit_function_derivative_vector(h, nlp, vars, settings, orginal_value, vars.x_bar)
				merit_derivative.dy = this.finite_difference_merit_function_derivative_vector(h, nlp, vars, settings, orginal_value, vars.y)
				merit_derivative.dy_bar = this.finite_difference_merit_function_derivative_vector(h, nlp, vars, settings, orginal_value, vars.y_bar)
				merit_derivative.ds = this.finite_difference_merit_function_derivative_vector(h, nlp, vars, settings, orginal_value, vars.s)
				merit_derivative.dz = this.finite_difference_merit_function_derivative_vector(h, nlp, vars, settings, orginal_value, vars.z)
				
				merit_derivative.dtau = this.finite_difference_merit_function_derivative_point(h, nlp, vars, settings, orginal_value, vars.tau)
				merit_derivative.dkappa = this.finite_difference_merit_function_derivative_point(h, nlp, vars, settings, orginal_value, vars.kappa)
				
				return(merit_derivative)
			catch e
				println("ERROR class_local_approximation.finite_difference_merit_function_derivative")
				throw(e)
			end
		end
		
		this.finite_difference_merit_function_derivative_vector = function(h, nlp, vars, settings, orginal_value, vector)
			result = zeros(length(vector))

			for i in 1:length(vector)
				orginal_point = vector[i];
				vector[i] = orginal_point + h
				f1 = this.update_approximation(nlp, vars, settings);
				vector[i] = orginal_point - h
				f2 = this.update_approximation(nlp, vars, settings);
				result[i] = (f1 - f2)/(2.0*h)
				
				vector[i] = orginal_point
			end
				
			return result
		end
		
		this.finite_difference_merit_function_derivative_point = function(h, nlp, vars, settings, orginal_value, point)
			point = point + h
			f1 = this.update_approximation(nlp, vars, settings);
			point = point - 2.0*h
			f2 = this.update_approximation(nlp, vars, settings);
			
			result = (f1 - f2)/(2.0*h)
			
			point = point + h;
				
			return result
		end
			
		this.calculate_potential_merit_function = function(nlp::class_non_linear_program,vars::class_variables,settings::class_settings)
			try
				#theta = 1.0/sqrt(nlp.m_1 + nlp.n_1 + 1);
				#boundary_distance = sum(log(vars.x)) + sum(log(vars.s)) + sum(log(vars.z)) + sum(log(vars.y)) + sum(log(vars.tau)) + sum(log(vars.kappa))
				#q = (nlp.m_1 + nlp.n_1 + 1) + sqrt(nlp.m_1 + nlp.n_1 + 1)
				#potential_function = 2*log(vars.x'*vars.s + vars.z'*vars.y + vars.tau*vars.kappa)[1] - boundary_distance
				#merit_function_value = -(this.mu)*boundary_distance + log(this.r_dual_norm) + log(this.r_primal_norm) + log(abs(this.r_gap));
		
				
					
				
				#merit_function_value = sqrt(nlp.m_1 + nlp.n_1 + 1)*(this.state.mu) + this.state.r_norm;
				
				#return merit_function_value;
			catch e
				println("ERROR in class_local_approximation.calculate_merit_function")
				throw(e)
			end
		end
		
		
		return(this);
	end
end
