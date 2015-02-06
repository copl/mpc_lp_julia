#Settings for the algorithm
type class_settings
    
    max_iter::Int            
    linear_feas_tol::Real   #This is a relative tolerance w.r.t. 
                                 #some normalizing norms
    comp_tol::Real          #How small must s^Tz must be when we stop

    #Constant length of the 
    #maximum combined step to the boundary to use
    bkscale::Real
    
    #Configuration for solver
    linear_solver_settings 

    function class_settings(max_iter::Int,linear_feas_tol::Real,comp_tol::Real,
                            bkscale::Real)
        this = new()
        this.max_iter = max_iter
        this.linear_feas_tol = linear_feas_tol
        this.comp_tol        = comp_tol
        this.bkscale         = bkscale
        return this
    end
end

type class_linear_program_input
	A
	G
	c
	h
	b
	m
	n
	k
	
	function class_linear_program_input()
		this = new();
		return(this);
	end
end

type class_K_newton_matrix
	update::Function
	
	function class_K_newton_matrix(problem_data::class_linear_program_input)
		this = new();
		this.update = function(variables)
			
		end
		return(this);
	end
end

type class_linear_program_result
	###########
	# TO DO
	###########
	result # integer codes: optimal_solution_found, unbounded, infeasible
	variables
end

type class_algorithm_state
	# TO DO
	mu::Real
	sigma::Real
	gap::Real
	
	update_mu::Function
	update_gap::Function
	
	function class_algorithm_state()
		this = new();
		
		this.update_mu = function(variables::class_linear_program_variables,problem_data::class_linear_program_input)
			this.mu = ( ((variables.s)'*(variables.z) + (variables.tau)*(variables.kappa))/(problem_data.m+1))[1] 
		end
		
		this.update_gap = function(variables::class_linear_program_variables,problem_data::class_linear_program_input)
			this.gap = ((problem_data.c)'*(variables.x) + (problem_data.h)'*(variables.z) + (problem_data.b)'*(variables.y))[1];
		end
		
		return(this)
	end
end

type class_linear_program_variables 
	x
	s
	z
	y
	tau::Real
	kappa::Real
	
	take_step::Function
	
	function class_linear_program_variables(problem_data::class_linear_program_input)
		this = new();
		
		this.x = zeros(problem_data.k, 1);
		this.s = ones(problem_data.m, 1);
		this.z = ones(problem_data.m, 1);
		this.y = zeros(problem_data.n, 1);
		this.tau = 1;
		this.kappa = 1;
		
		this.take_step = function(direction::class_direction)
			alpha = direction.alpha; # get the step size
			
			this.x = this.x + alpha*direction.dx;
			this.s = this.s + alpha*direction.ds;
			this.y = this.y + alpha*direction.dy;
			this.z = this.z + alpha*direction.dz;
			this.tau = this.tau + alpha*direction.dtau;
			this.kappa = this.kappa + alpha*direction.dkappa;
		end
		
		return(this)
	end
end

type class_linear_system_rhs
	q1
	q2
	q3
	q4
	q5
	q6
	
	update_values::Function
	compute_affine_rhs::Function
	compute_corrector_rhs::Function
	
	function class_linear_system_rhs(problem_data::class_linear_program_input)
		this = new();
		
		this.update_values = function (q1,q2,q3,q4,q5,q6)
			this.q1 = q1
			this.q2 = q2
			this.q3 = q3
			this.q4 = q4
			this.q5 = q5
			this.q6 = q6
		end
		
		this.compute_affine_rhs = function(residuals::class_residuals,	variables::class_linear_program_variables)
			this.update_values(-residuals.r1, -residuals.r2, -residuals.r3, -residuals.r4, -(variables.z).*(variables.s), -(variables.tau)*(variables.kappa))
		end
		
		this.compute_corrector_rhs = function(residuals::class_residuals,variables::class_linear_program_variables,state::class_algorithm_state,affine_direction::class_direction,problem_data::class_linear_program_input)
			m = problem_data.m
			
			z = variables.z
			s = variables.s
			tau = variables.tau
			kappa = variables.kappa
			
			dx_a = affine_direction.dx
			dy_a = affine_direction.dy
			ds_a = affine_direction.ds
			dz_a = affine_direction.dz
			dtau_a = affine_direction.dtau
			dkappa_a = affine_direction.dkappa
			alpha = affine_direction.alpha
			
			mu = state.mu
			
			mu_a = ((s+alpha*ds_a)'*(z+alpha*dz_a) + (tau + alpha*dtau_a)*((kappa) + alpha*dkappa_a))/(m+1);
			sigma = ((mu_a/(mu))^3)[1]
			
			state.sigma = sigma
			
			this.update_values(-(1-sigma)*residuals.r1, -(1-sigma)*residuals.r2, -(1-sigma)*residuals.r3, -(1-sigma)*residuals.r4, -z.*s -ds_a.*dz_a + sigma*mu,  -tau*kappa-dtau_a*dkappa_a + sigma*mu)
		end
		
		return this
	end
end

type class_residuals
	r1
	r2
	r3
	r4
	
	r1_norm::Real
	r2_norm::Real
	r3_norm::Real
	normed_squared::Real
	
	
	update_values::Function
	compute_residuals::Function
	
	function class_residuals(problem_data::class_linear_program_input)
		this = new();
		
		this.update_values = function(r1,r2,r3,r4)
			this.r1 = r1
			this.r2 = r2
			this.r3 = r3
			this.r4 = r4
			
			this.r1_norm = Base.norm(r1)
			this.r2_norm = Base.norm(r2)
			this.r3_norm = Base.norm(r3)
			
			this.normed_squared = 0
		end
		
		this.compute_residuals = function(pd::class_linear_program_input, variables::class_linear_program_variables)
			if false:
				r = [zeros(1,pd.k) zeros(1,pd.n) (variables.s)'  variables.kappa]' -
						[ zeros(pd.k,pd.k)  pd.A'          pd.G'           pd.c;
						 -pd.A           zeros(pd.n,pd.n)  zeros(pd.n,pd.m)  pd.b;
						 -pd.G           zeros(pd.m,pd.n)  zeros(pd.m,pd.m)  pd.h;
						 -pd.c'          -pd.b'         -pd.h'          0]*
						 [variables.x; variables.y; variables.z; variables.tau];

				r1 = r[1:pd.k];
				r2 = r[(pd.k+1):(pd.k+pd.n)];
				r3 = r[(pd.k+pd.n+1):(pd.k+pd.n+pd.m)];
				r4 = r[pd.k+pd.n+pd.m+1];
			
			this.r1 = pd.A'*variables.y + pd.G'*variables.z + pd.c*variables.tau;
			
			#this.update_values(r1,r2,r3,r4)
		end
		
		return this
	end
end

type class_direction
	dx
	dy
	dz
	ds
	dtau::Real
	dkappa::Real
	alpha::Real
	
	update_values::Function
	compute_affine_direction::Function
	compute_corrector_direction::Function
	
	function class_direction(problem_data::class_linear_program_input)
		this = new();
		
		this.update_values = function(dx,dy,dz,dtau,ds,dkappa,alpha)
			this.dx = dx;
			this.dy = dy;
			this.dz = dz;
			this.dtau = dtau;
			this.ds = ds
			this.dkappa = dkappa;
			this.alpha = alpha;
		end
		
		this.compute_affine_direction = function(affine_rhs::class_linear_system_rhs,
													problem_data::class_linear_program_input,	
													variables::class_linear_program_variables,	
													K_newton_matrix::class_K_newton_matrix)
			
			dir = solveLinearEquation(problem_data, variables, affine_rhs, K_newton_matrix)
			
			m = problem_data.m
			n = problem_data.n
			k = problem_data.k
			
			x = variables.x
			z = variables.z
			s = variables.s
			y = variables.y
			tau = variables.tau
			kappa = variables.kappa
			
			dx_a = dir[1:k];
			dy_a = dir[(k+1):(k+n)];
			dz_a = dir[(k+n+1):(k+n+m)];
			dtau_a = dir[(k+n+m+1)];
			ds_a = ( -z.*s - dz_a.*s)./z;
			dkappa_a = (-(tau)*(kappa) - dtau_a*(kappa))/(tau)

			# Compute Step size a, and Centering Parameter s
			vv = vec([s./ds_a; z./dz_a; tau/dtau_a; kappa/dkappa_a]);
			alpha = 1;
			for i=1:length(vv)
				if (vv[i] < 0)
					alpha = minimum([alpha, -vv[i]]);
				end
			end
			
			this.update_values(dx_a,dy_a,dz_a,dtau_a,ds_a,dkappa_a,alpha)
		end
		
		this.compute_corrector_direction = function(
									 affine_direction::class_direction,
									 corrector_rhs::class_linear_system_rhs,
									 problem_data::class_linear_program_input,
                                     variables::class_linear_program_variables,
									 state::class_algorithm_state,
                                     settings::class_settings,
									 K_newton_matrix::class_K_newton_matrix)
			m = problem_data.m
			n = problem_data.n
			k = problem_data.k
			
			kappa = variables.kappa
			tau = variables.tau
			s = variables.s
			z = variables.z
			
			dx_a = affine_direction.dx
			dy_a = affine_direction.dy
			ds_a = affine_direction.ds
			dz_a = affine_direction.dz
			dtau_a = affine_direction.dtau
			dkappa_a = affine_direction.dkappa
			alpha = affine_direction.alpha
			
			mu = state.mu
			sigma = state.sigma
			
			dir = solveLinearEquation(problem_data, variables, corrector_rhs, K_newton_matrix);
			
			dx = dir[1:k];
			dy = dir[(k+1):(k+n)];
			dz = dir[(k+n+1):(k+n+m)];
			dtau = dir[(k+n+m+1)];
			ds = ( -z.*s -ds_a.*dz_a + sigma*mu - dz.*s)./z;
			dkappa = (-tau*kappa-dtau_a*dkappa_a + sigma*mu - dtau*kappa)/tau

			# Update
			vv = vec([s./ds; z./dz; tau/dtau; kappa/dkappa]);
			alpha = 1;
			for i=1:length(vv)
				if (vv[i] < 0)
					alpha = minimum([alpha, -vv[i]]);
				end
			end
			alpha = alpha*settings.bkscale
			
			this.update_values(dx,dy,dz,dtau,ds,dkappa,alpha)
		end
		
		return this
	end
end