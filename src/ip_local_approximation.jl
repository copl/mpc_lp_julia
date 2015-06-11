println("loading ip_local_approximation")

type class_state
	r_dual_norm::Float64
	r_primal_norm::Float64
	r_gap_norm::Float64
	r_norm::Float64
	mu::Float64
	
	dual_feasibility::Float64
	primal_feasibility::Float64
	dual_infeasibility::Float64
	primal_infeasibility::Float64
	
	
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
				
				this.dual_feasibility = this.r_dual_norm / vars.tau
				this.primal_feasibility = this.r_primal_norm / vars.tau;
				
				homogenous_dual_objective = -(local_approx.v3' * [vars.y; vars.y_bar])[1];
				homogeneous_primal_objective = -(local_approx.c' * [vars.x; vars.x_bar])[1];
				
				this.dual_infeasibility = norm(local_approx.J * [vars.x; vars.x_bar] - [vars.z; zeros(length(vars.y_bar))], GLOBAL_P_NORM ) / homogeneous_primal_objective;
				this.primal_infeasibility = norm( (local_approx.J)' * [vars.y; vars.y_bar] + [vars.s; zeros(
				length(vars.x_bar))], GLOBAL_P_NORM ) / homogenous_dual_objective;
				
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
	
	update_approximation::Function # update the local approximation at the new point
	calculate_merit_function::Function
	calculate_merit_function_derivative::Function
	
	gamma::Float64
	
	c::Array{Float64,1}
	g::Array{Float64,1}
	J::SparseMatrixCSC{Float64,Int64}
	H::SparseMatrixCSC{Float64,Int64}
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
		this.state = class_state();
		
		this.update_approximation = function(nlp::class_non_linear_program,vars::class_variables,settings::class_settings)
			try
				this.current_objective_value = nlp.objective_function(vars.x_scaled);
				
				this.c = nlp.objective_function_gradient(vars.x_scaled);
				this.J = nlp.evaluate_constraint_gradients(vars.x_scaled);
				this.g = this.c - this.J'*vars.y_scaled;
				this.H = nlp.objective_function_hessian(vars.x_scaled) - nlp.evaluate_constraint_lagrangian_hessian(vars.x_scaled,vars.y_scaled);
			
				Delta_H = this.H * (vars.x_scaled);
				
				# In (30)
				this.a  = nlp.evaluate_constraints(vars.x_scaled);
				this.v1 = -this.c - Delta_H; # ERROR IN YINYU ANDERSON 1
				this.v2 = this.c - Delta_H; # ERROR IN YINYU ANDERSON 1
				this.v3 = this.a - this.J*(vars.x_scaled);
				
				# following (31)
				this.D_g = ((vars.x_scaled'*Delta_H)[1] + vars.kappa/vars.tau);
				this.D_x = this.H + spdiagm([(vars.s)./(vars.x); zeros(length(vars.x_bar))]);
				#this.D_z = spdiagm([(vars.y)./(vars.z); zeros(length(vars.y_bar))]);
				this.D_z = spdiagm([(vars.z)./(vars.y); zeros(length(vars.y_bar))]);
				
				# (27) the residuals of our problem
				this.r_dual = ([vars.s; zeros(length(vars.x_bar))] - vars.tau*(this.g));
				this.r_gap = (vars.kappa + vars.tau*((vars.x_scaled'*this.g) + (vars.y_scaled'*this.a)))[1];
				this.r_primal = ([vars.z; zeros(length(vars.y_bar))] - vars.tau*this.a);
				
				this.state.update_state(this,vars,nlp)
				
				return this.calculate_merit_function(nlp,vars,settings);
			catch e
				println("ERROR in class_local_approximation.update_approximation")
				throw(e)
			end
		end
		
		this.calculate_merit_function = function(nlp::class_non_linear_program,vars::class_variables,settings::class_settings)
			try
				#theta = 1.0/sqrt(nlp.m_1 + nlp.n_1 + 1);
				#boundary_distance = sum(log(vars.x)) + sum(log(vars.s)) + sum(log(vars.z)) + sum(log(vars.y)) + sum(log(vars.tau)) + sum(log(vars.kappa))
				#q = (nlp.m_1 + nlp.n_1 + 1) + sqrt(nlp.m_1 + nlp.n_1 + 1)
				#potential_function = 2*log(vars.x'*vars.s + vars.z'*vars.y + vars.tau*vars.kappa)[1] - boundary_distance
				#merit_function_value = -(this.mu)*boundary_distance + log(this.r_dual_norm) + log(this.r_primal_norm) + log(abs(this.r_gap));
		
				closest_point = minimum([minimum([vars.x, [Inf]]),minimum([vars.s, [Inf]]),minimum([vars.z, [Inf]]),minimum([vars.y, [Inf]]),vars.tau,vars.kappa])
					
				if closest_point < settings.beta2*this.state.mu
					merit_function_value = Inf;
				else
					merit_function_value = sqrt(nlp.m_1 + nlp.n_1 + 1)*(this.state.mu) + this.state.r_norm;
				end
				
				return merit_function_value;
			catch e
				println("ERROR in class_local_approximation.calculate_merit_function")
				throw(e)
			end
		end
		
		this.calculate_merit_function_derivative = function(nlp::class_non_linear_program,vars::class_variables,settings::class_settings)
			try
				sigma_D = GLOBAL_P_NORM * (this.r_dual) .^ GLOBAL_P_NORM;
				sigma_G = GLOBAL_P_NORM * (this.r_gap) .^ GLOBAL_P_NORM;
				sigma_P = GLOBAL_P_NORM * (this.r_primal) .^ GLOBAL_P_NORM;
				
				derivative_merit_function = class_direction();
				
				dx = -this.H * sigma_D + sigma_G * ( this.g + (this.J)' * vars.y_scaled ) - (this.J)' * sigma_P; 
				dy = this.J * sigma_D + sigma_G * ( this.J * vars.x_scaled + this.a ); 
				
				derivative_merit_function.dx = dx[1:length(vars.x)];
				derivative_merit_function.dx_bar = dx[(length(vars.x)+1):length(vars.x_scaled)];
				derivative_merit_function.dy = dy[1:length(vars.y)];
				derivative_merit_function.dy_bar = dy[(length(vars.y)+1):length(vars.y_scaled)];
				
				derivative_merit_function.ds = sigma_D[1:length(vars.x)];
				derivative_merit_function.dz = sigma_P[1:length(vars.y)];
				derivative_merit_function.dtau = -( this.g' * sigma_D + this.a ' * sigma_P )[1]; 
				derivative_merit_function.dkappa = sigma_G;
				
				derivative_merit_function.validate_direction_dimensions(vars)
				
				return(derivative_merit_function)
			catch e
				println("ERROR class_local_approximation.calculate_merit_function_derivative")
				throw(e)
			end
		end
		
		
		
		return(this);
	end
end