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
	P::SparseMatrixCSC{Float64,Int64} # matrix
	A::SparseMatrixCSC{Float64,Int64} # matrix
	G::SparseMatrixCSC{Float64,Int64} # matrix
	c::Array{Float64,1} # vector
	h::Array{Float64,1} # vector
	b::Array{Float64,1} # vector
	m::Int
	n::Int
	k::Int
	
	function class_linear_program_input()
		this = new();
		return(this);
	end
end

type class_K_newton_matrix
	update::Function
	update2::Function
	solve::Function
	solve2::Function
	
	n::Real
	m::Real
	k::Real
	
	Q::SparseMatrixCSC{Float64,Int64}
	#Q2::SparseMatrixCSC{Float64,Int64}
	F #::UmfpackLU{Float64,Int64} #::CholmodFactor
	#F2 #::CholmodFactor{Float64,Int64}
	#F2_defined::Bool
	
	function class_K_newton_matrix(problem_data::class_linear_program_input,settings::class_settings)
		this = new();
		n = problem_data.n;
		m = problem_data.m;
		k = problem_data.k;
		
		this.n = n;
		this.m = m;
		this.k = k;
		
		this.Q = sparse([ problem_data.P  problem_data.A'          problem_data.G'    sparse(problem_data.c);
				 -problem_data.A           spzeros(n,n)  spzeros(n, m)      sparse(problem_data.b);
				 -problem_data.G           spzeros(m,n)  speye(m)  sparse(problem_data.h);
				sparse(-problem_data.c')         sparse(-problem_data.b')        sparse(-problem_data.h')       sparse([1.0])
				]);
		
				
		this.update = function(variables::class_linear_program_variables)
			n = this.n
			m = this.m
			
			offset = k + n;
			for i = 1:m
				this.Q[offset + i,offset + i] = variables.s[i]/variables.z[i];
			end
			
			this.Q[offset + m + 1, offset + m + 1] = variables.kappa/variables.tau;
			
			this.F = lufact(this.Q);
			#@time lufact(this.Q2);
			
			#this.update2(variables);
		end
		
		this.solve = function(problem_data::class_linear_program_input,variables::class_linear_program_variables,rhs::class_linear_system_rhs)
			q7 = rhs.q3 - rhs.q5./variables.z;
			q8 = rhs.q4 - rhs.q6./variables.tau;
						
			result = -(this.F\[rhs.q1; rhs.q2; q7; q8])
			
			return(result)
		end
		
		this.solve2 = function(problem_data::class_linear_program_input,variables::class_linear_program_variables,rhs::class_linear_system_rhs)
			real_rhs = [-problem_data.c; problem_data.b problem_data.h]
			
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
	x::Array{Float64,1} # vector
	s::Array{Float64,1} # vector
	z::Array{Float64,1} # vector
	y::Array{Float64,1} # vector
	tau::Real
	kappa::Real
	
	take_step::Function
	
	function class_linear_program_variables(problem_data::class_linear_program_input)
		this = new();
		
		this.x = zeros(problem_data.k);
		this.s = ones(problem_data.m);
		this.z = ones(problem_data.m);
		this.y = zeros(problem_data.n);
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
	q1::Array{Float64,1} # vector
	q2::Array{Float64,1} # vector
	q3::Array{Float64,1} # vector
	q4::Array{Float64,1} # vector
	q5::Array{Float64,1} # vector
	q6::Real # vector
	
	update_values::Function
	compute_affine_rhs::Function
	compute_corrector_rhs::Function
	
	function class_linear_system_rhs(problem_data::class_linear_program_input)
		this = new();
		
		this.compute_affine_rhs = function(residuals::class_residuals,	variables::class_linear_program_variables)
			this.q1 = -residuals.r1
			this.q2 = -residuals.r2
			this.q3 = -residuals.r3
			this.q4 = -residuals.r4
			this.q5 = -(variables.z).*(variables.s)
			this.q6 = -(variables.tau)*(variables.kappa)
		
			debug_message("q1")
			debug_message(this.q1)
			debug_message("q2")
			debug_message(this.q2)
			debug_message("q3")
			debug_message(this.q3)
			debug_message("q4")
			debug_message(this.q4)
			debug_message("q5")
			debug_message(this.q5)
			debug_message("q6")
			debug_message(this.q6)
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
			
			mu_a = (((s+alpha*ds_a)'*(z+alpha*dz_a))[1] + (tau + alpha*dtau_a)*((kappa) + alpha*dkappa_a))/(m+1);
			sigma = ((mu_a/(mu))^3)
			
			state.sigma = sigma
			
			this.q1 = -(1-sigma)*residuals.r1
			this.q2 = -(1-sigma)*residuals.r2
			this.q3 = -(1-sigma)*residuals.r3 
			this.q4 = -(1-sigma)*residuals.r4
			this.q5 = -z.*s -ds_a.*dz_a + sigma*mu
			this.q6 = -tau*kappa-dtau_a*dkappa_a + sigma*mu
		end
		
		return this
	end
end

type class_residuals
	r1 #::Array{Float64,1} # vector
	r2 #::Array{Float64,1} # vector
	r3 #::Array{Float64,1} # vector
	r4 #::Array{Float64,1} # vector
	
	r1_norm::Real
	r2_norm::Real
	r3_norm::Real
	residuals_norm::Real
	
	
	update_norms::Function
	compute_residuals::Function
	
	function class_residuals(problem_data::class_linear_program_input)
		this = new();
		
		this.update_norms = function()
			this.r1_norm = Base.norm(this.r1)
			this.r2_norm = Base.norm(this.r2)
			this.r3_norm = Base.norm(this.r3)
			
			this.residuals_norm = sqrt((this.r1_norm)^2 + (this.r2_norm)^2 + (this.r3_norm)^2)
		end
		
		this.compute_residuals = function(pd::class_linear_program_input, variables::class_linear_program_variables)
			this.r1 = -pd.A'*variables.y - pd.G'*variables.z - pd.c*variables.tau;
			this.r2 = pd.A*variables.x - pd.b*variables.tau;
			this.r3 = variables.s + pd.G*variables.x - variables.tau*pd.h;
			this.r4 = variables.kappa + pd.c'*variables.x + pd.b'*variables.y + pd.h'*variables.z;
			
			debug_message("r1")
			debug_message(this.r1)
			debug_message("r2")
			debug_message(this.r2)
			debug_message("r3")
			debug_message(this.r3)
			debug_message("r4")
			debug_message(this.r4)
			
			this.update_norms()
		end
		
		return this
	end
end

type class_direction
	dx::Array{Float64,1} # vector
	dy::Array{Float64,1} # vector
	dz::Array{Float64,1} # vector
	ds::Array{Float64,1} # vector
	dtau::Real
	dkappa::Real
	alpha::Real
	
	update_values::Function
	compute_affine_direction::Function
	compute_corrector_direction::Function
	
	compute_alpha::Function
	compute_min_ratio_alpha::Function
	
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
									 state::class_algorithm_state,
                                     settings::class_settings,
									 K_newton_matrix::class_K_newton_matrix)
			
			dir = K_newton_matrix.solve(problem_data, variables, affine_rhs)
			
			m = problem_data.m
			n = problem_data.n
			k = problem_data.k
			
			x = variables.x
			z = variables.z
			s = variables.s
			y = variables.y
			tau = variables.tau
			kappa = variables.kappa
			
			this.dx = dir[1:k];
			this.dy = dir[(k+1):(k+n)];
			this.dz = dir[(k+n+1):(k+n+m)];
			this.dtau = dir[(k+n+m+1)];
			this.ds = ( -z.*s - this.dz.*s)./z;
			this.dkappa = (-(tau)*(kappa) - this.dtau*(kappa))/(tau)

			this.compute_alpha(variables,settings)
		end
		
		this.compute_corrector_direction = function(
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
			
			 # this requires clever manipulation in C to prevent use
			dx_a = this.dx
			dy_a = this.dy
			ds_a = this.ds
			dz_a = this.dz
			dtau_a = this.dtau
			dkappa_a = this.dkappa
			alpha = this.alpha
			
			mu = state.mu
			sigma = state.sigma
			
			dir = K_newton_matrix.solve(problem_data, variables, corrector_rhs);
			
			this.dx = dir[1:k];
			this.dy = dir[(k+1):(k+n)];
			this.dz = dir[(k+n+1):(k+n+m)];
			this.dtau = dir[(k+n+m+1)];
			this.ds = ( -z.*s -ds_a.*dz_a + sigma*mu - (this.dz).*s)./z;
			this.dkappa = (-tau*kappa-dtau_a*dkappa_a + sigma*mu - (this.dtau)*kappa)/tau

			# Update
			this.compute_alpha(variables,settings)
		end
		
		this.compute_alpha = function(variables::class_linear_program_variables,settings::class_settings)
			this.alpha = 99999;
			this.compute_min_ratio_alpha(variables.s,this.ds)
			this.compute_min_ratio_alpha(variables.z,this.dz)
			this.compute_min_ratio_alpha([variables.kappa],[this.dkappa])
			this.compute_min_ratio_alpha([variables.tau],[this.dtau])
			this.alpha = this.alpha*settings.bkscale
			this.alpha = min(1.0,this.alpha)
		end
		
		this.compute_min_ratio_alpha = function(var,dvar)
			for i=1:length(var)
				candidate_alpha = -var[i]/dvar[i]
				if candidate_alpha > 0
					this.alpha = minimum([this.alpha, candidate_alpha]) 
				end
			end
		end
		return this
	end
end