function sparse_condition_number(M)
	return sqrt(abs(eigs(M'*M; nev=1, which=:LM, ritzvec=false)[1][1] / eigs(M'*M; nev=1, which=:SM, ritzvec=false)[1][1]));
end

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
	grad::Array{Float64,1}
	L::SparseMatrixCSC{Float64,Int64}
	
	P::SparseMatrixCSC{Float64,Int64} # matrix
	A::SparseMatrixCSC{Float64,Int64} # matrix
	G::SparseMatrixCSC{Float64,Int64} # matrix
	c::Array{Float64,1} # vector
	h::Array{Float64,1} # vector
	b::Array{Float64,1} # vector
	m::Int
	n::Int
	k::Int
	
	barrier_weights::Array{Float64,1}
	
	currentObjectiveValue::Real
	
	evaluate_objective::Function
	
	black_box_objective::Bool
	objective_function::Function
	objective_function_gradient::Function
	objective_function_hessian::Function
	
	black_box_constraints::Bool
	evaluate_inequality_constraints::Function
	evaluate_inequality_constraint_gradients::Function
	evaluate_inequality_constraint_lagrangian::Function
	
	update::Function
	add_black_box_objective::Function
	add_black_box_constraints::Function
	
	validate::Function
	
	add_inequality_constraint::Function
	
	function class_linear_program_input()
		this = new();
		
		this.black_box_objective = false;
		this.black_box_constraints = false;
		
		this.evaluate_objective = function(x)
			if this.black_box_objective
				return this.objective_function(x)
			else
				return this.c'*x + 0.5*x'*this.P*x
			end
		end
		
		this.add_inequality_constraint = function(constraint,rhs,barrier_weight)
			if ~this.black_box_constraints
				this.m = this.m + 1;
				this.G = [this.G; constraint]; 
				this.h = [this.h; rhs];
				this.barrier_weights = [this.barrier_weights, barrier_weight];
			else
				Error("black box constraints and adding constraints is dangerous!")
			end
		end
		
		this.add_black_box_objective = function(objective_function::Function, objective_function_gradient::Function, objective_function_hessian::Function)
			this.black_box_objective = true;
			
			this.objective_function = objective_function;
			this.objective_function_gradient = objective_function_gradient;
			this.objective_function_hessian = objective_function_hessian;
		end
		
		this.add_black_box_constraints = function(evaluate_inequality_constraints::Function, evaluate_inequality_constraint_gradients::Function, evaluate_inequality_constraint_lagrangian::Function)
			this.black_box_constraints = true;
			
			this.evaluate_inequality_constraints = evaluate_inequality_constraints;
			this.evaluate_inequality_constraint_gradients = evaluate_inequality_constraint_gradients;
			this.evaluate_inequality_constraint_lagrangian = evaluate_inequality_constraint_lagrangian;
		end
		
		this.update = function(variables::class_linear_program_variables)
			try
				xNormalized = variables.x/variables.tau;
				zNormalized = variables.z/variables.tau;
				yNormalized = variables.y/variables.tau;
				
				if this.black_box_objective # update the objective
					this.currentObjectiveValue = this.objective_function(xNormalized);
					this.P = this.objective_function_hessian(xNormalized);
					this.grad = this.objective_function_gradient(xNormalized);
					# - 0.5*(this.P)*(xNormalized);
				else
					this.grad = this.c + 0.5*(this.P)*(xNormalized);
					this.currentObjectiveValue = ((this.c')*xNormalized +  0.5*(xNormalized')*(this.P)*(xNormalized))[1];
				end
				
				this.L = this.P; # default
				if this.black_box_constraints # update the constraints
					this.G = this.evaluate_inequality_constraint_gradients(xNormalized);
					this.h = this.G*xNormalized - this.evaluate_inequality_constraints(xNormalized);
					
					this.L = this.L + this.evaluate_inequality_constraint_lagrangian(xNormalized,zNormalized);
				end
								 
			catch e
				println("*** ERROR in class_linear_program_input.update ***")
				throw(e)
			end
		end
		
		this.validate = function()
			# check problem data has been correctly created
			
			# check dimensions			
			try
				if (size(this.P) != (this.k,this.k))
					error("size(this.P) != (this.k,this.k)")
				end
				
				if size(this.A) != (this.n,this.k)
					error("size(this.A) != (this.n,this.k)")
				end
				
				if size(this.G) != (this.m,this.k)
					error("size(this.G) != (this.m,this.k)")
				end
				
				if size(this.c) != (this.k,)
					error("size(this.c) != (this.k,)")
				end
				
				if size(this.b) != (this.n,)
					error("size(this.b) != (this.n,)")
				end
				
				if size(this.h) != (this.m,)
					error("size(this.h) != (this.m,)")
				end
				
				if size(this.barrier_weights) != (this.m,)
					error("size(this.barrier_weights) != (this.m,)")
				end
			catch e
				println("*** Error while validating problem data ***")
				throw(e)
			end
			# TO DO: check derivatives
		end
		
		return(this);
	end
end

type class_K_newton_matrix
	update::Function
	block_update::Function
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
		
		#println(cond(full(problem_data.A*problem_data.A')))
		
		#problem_data.P
		#rho = 1.0
		#problem_data.P = speye(k);
		
				
		#this.K = sparse([ 	
		#					problem_data.P   	problem_data.A'     problem_data.G'		spzeros(k,1);
		#					problem_data.A      spzeros(n,n)  		spzeros(n, m) 		spzeros(n,1);
		#					problem_data.G      spzeros(m,n)  		-speye(m)			spzeros(m,1);
		#					spzeros(1,k)        spzeros(1,n)       	spzeros(1,m)		sparse([1.0])
		#				]);
		
		#Q2 = sparse([ problem_data.G'*speye(m)*problem_data.G   problem_data.A'          ;
		#		 problem_data.A           spzeros(n,n)
		#		]);
				
		#Q3 = problem_data.A*problem_data.A' + 10^(-6.0)*speye(n,n);
		
		#println("rank etc")
		#println(n)
		#println(rank(full(problem_data.A)))
		#this.F = @time lufact(Q2);
		
		#this.F = @time cholfact(Q3);
		
		#this.solve_for = function(Factor,u,v,rhs)
		#	temp = Factor \ rhs;
		#	result = temp - Factor \ temp
		#end
		
		#this.update_Q = function(variables::class_linear_program_variables)
		#	this.Q[]
		#end
		
		this.block_update = function(variables::class_linear_program_variables,block_length)
			n = this.n
			m = this.m
			k = this.k
			
			diag_vars = -vec((variables.s) ./ (variables.z));
			#println(variables.s)
			#println(variables.z)
			#println(diag_vars)
			offset = k + n;
			sorted_indicies = (sortperm(diag_vars))
			#sorted_indicies = [1,2,3,4];
			
			number_of_blocks = ceil(length(sorted_indicies)/block_length);
			block_start = 1
			block_end = block_length;
			for block = 1:number_of_blocks
				block_size = 1 + block_end - block_start;
				values = diag_vars[sorted_indicies[block_start:block_end]]
				average_value_of_block = mean(values)
				#println(values)
				#println(average_value_of_block)
				#println(offset + (block_start:block_end))
				#println(size(this.Q))
				indicies = sorted_indicies[block_start:block_end]
				this.Q[offset + indicies, offset + indicies] = (average_value_of_block)*ones(block_size, block_size)
				#this.Q[offset + indicies, offset + indicies] = average_value_of_block # -speye(block_size) + 
				block_start = block_end + 1; #- speye(block_size)
				block_end = min((block_end + block_length),length(sorted_indicies));
				#println(block_end)
				#println(length(sorted_indicies))
			end
		end
		
		this.update = function(problem_data::class_linear_program_input,variables::class_linear_program_variables)
			n = this.n
			m = this.m
			k = this.k
			
			# update diagonal
			
			#sortperm
			 #  + (0.1 - minimum(problem_data.P))
			this.Q = sparse([ problem_data.L   problem_data.A'          problem_data.G'    sparse(problem_data.c);
				 problem_data.A           spzeros(n,n)  spzeros(n, m)      -sparse(problem_data.b);
				 problem_data.G           spzeros(m,n)  -spdiagm((variables.s)./variables.z)  -sparse(problem_data.h);
				sparse(-problem_data.c')         sparse(-problem_data.b')        sparse(-problem_data.h')      variables.kappa/variables.tau
				]);
			
			#this.block_update(variables,2)
			
			# update gradient
			#last_entry = k + n + m + 1
			#grad = problem_data.P*(variables.x/variables.tau) + problem_data.c;
			
			#for i = 1:k
			#	this.Q[last_entry, i] = -grad[i] 
			#	this.Q[i, last_entry] = grad[i]
			#end
			
			M = this.Q;
			
			
			
			#condition_number = cond(full(M))
			#condition_number = sparse_condition_number(M);
			#println("Condition number " * string(condition_number))
			
			this.F = lufact(this.Q);
			#@time lufact(this.Q2);
			
			#this.update2(variables);
		end
		
		this.solve = function(problem_data::class_linear_program_input,variables::class_linear_program_variables,rhs::class_linear_system_rhs)
			q7 = rhs.q3 - rhs.q5./variables.z;
			q8 = rhs.q4 - rhs.q6./variables.tau;
						
			result = -(this.F\[rhs.q1; -rhs.q2; -q7; q8])
			
			return(result)
		end
		
		#this.solve2 = function(problem_data::class_linear_program_input,variables::class_linear_program_variables,rhs::class_linear_system_rhs)
		#	real_rhs = [-problem_data.c; problem_data.b problem_data.h]
		#	
		#end
		
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
			#problem_data.barrier_weights.*
			this.mu = ( ((variables.s)'*(variables.z) + (variables.tau)*(variables.kappa))/(problem_data.m+1))[1] 
		end
		
		this.update_gap = function(variables::class_linear_program_variables,problem_data::class_linear_program_input)
			#grad = problem_data.P*(variables.x/variables.tau) + problem_data.c;
			this.gap = ((problem_data.grad)'*(variables.x) + (problem_data.h)'*(variables.z) + (problem_data.b)'*(variables.y))[1];
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
	check::Function
	
	function class_linear_program_variables(problem_data::class_linear_program_input)
		this = new();
		
		this.s = ones(problem_data.m);
		this.z = ones(problem_data.m);# + rand(problem_data.m);
		this.x = 20*rand(problem_data.k); #zeros(problem_data.k)*
		#this.s = rand(problem_data.m) + 1;
		#this.z = rand(problem_data.m) + 1;
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
		
		this.check = function()
			try
				for val in this.s
					assert(val > 0)
				end
				
				for val in this.z
					assert(val > 0)
				end
			catch e
				println("one variable is less than or equal to zero")
				throw(e)
			end
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
			try
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
			catch e
				println("*** ERROR in compute_affine_rhs ***")
				throw(e)
			end
		end
		
		this.compute_corrector_rhs = function(residuals::class_residuals,variables::class_linear_program_variables,state::class_algorithm_state,affine_direction::class_direction,problem_data::class_linear_program_input)
			try
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
				this.q5 = -(z).*s -ds_a.*dz_a + sigma*mu*problem_data.barrier_weights #problem_data.barrier_weights.*
				this.q6 = -tau*kappa-dtau_a*dkappa_a + sigma*mu
				
			catch e
				println("*** ERROR in  compute_corrector_rhs ***")
				throw(e)
			end
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
	r4_norm::Real	
	summary_norm::Real
	kkt_is_infeasible::Real
	
	primal_feasibility::Real
	dual_feasibility::Real
	
	update_norms::Function
	compute_residuals::Function
	
	function class_residuals(problem_data::class_linear_program_input)
		this = new();
		
		this.compute_residuals = function(pd::class_linear_program_input, variables::class_linear_program_variables)
			try
				#grad = pd.P*(variables.x/variables.tau) + problem_data.c;
				# this.r1 =  grad*variables.tau
				
				this.r1 = -pd.A'*variables.y - pd.G'*variables.z - pd.grad*variables.tau;
				this.r2 = pd.A*variables.x - pd.b*variables.tau;
				this.r3 = variables.s + (pd.G*variables.x) - variables.tau*pd.h; 
				this.r4 = variables.kappa + pd.grad'*variables.x + pd.b'*variables.y + pd.h'*variables.z;
				
				# 5.4 pg 395 homogenous; A homogeneous algorithm for MCP Anderson and Ye
				this.r1_norm = norm(this.r1,1)/(norm(variables.y,1) + norm(variables.z,1) + variables.tau);
				this.r2_norm = norm(this.r2,1)/(norm(variables.x,1) + variables.tau);
				this.r3_norm = norm(this.r3,1)/(norm(variables.s,1) + norm(variables.x,1) + variables.tau);
				this.r4_norm = norm(this.r4,1)/(variables.kappa + norm(variables.x,1) + norm(variables.z,1));
				
				this.kkt_is_infeasible = variables.tau/variables.kappa;
				
				this.summary_norm = this.r1_norm + this.r2_norm + this.r3_norm + this.r4_norm;
				
				this.primal_feasibility = (pd.b'*variables.y + pd.h'*variables.z)[1];
				this.dual_feasibility = (pd.grad'*variables.x)[1];
				
				#this.primal_feasibility = pd.A
				#this.dual_feasibility = this.r1_norm/(variables.tau + norm(variables.s,1))
				#this.relative_gap = grad'*variables.x + pd.b'*variables.y + pd.h'*variables.z				
				
				debug_message("r1")
				debug_message(this.r1)
				debug_message("r2")
				debug_message(this.r2)
				debug_message("r3")
				debug_message(this.r3)
				debug_message("r4")
				debug_message(this.r4)
			catch (e)
				println("*** ERROR in compute_residuals")
				throw(e)
			end
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
	line_search_alpha::Function
	
	function class_direction(problem_data::class_linear_program_input)
		this = new();
		this.alpha = -1;
		
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
			try
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
			catch e
				println("*** ERROR in compute_affine_direction ***")
				throw(e)
			end
		end
		
		this.compute_corrector_direction = function(
									 corrector_rhs::class_linear_system_rhs,
									 problem_data::class_linear_program_input,
                                     variables::class_linear_program_variables,
									 state::class_algorithm_state,
                                     settings::class_settings,
									 K_newton_matrix::class_K_newton_matrix)
									 
			try
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
				this.ds = ( -z.*s -ds_a.*dz_a + sigma*mu*problem_data.barrier_weights - (this.dz).*s )./z;
				this.dkappa = ( -tau*kappa-dtau_a*dkappa_a + sigma*mu - (this.dtau)*kappa )/tau

				this.compute_alpha(variables,settings)
			catch e
				println("*** ERROR in  compute_corrector_direction ***")
				throw(e)
			end
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
		
		this.line_search_alpha = function(x,problem_data::class_linear_program_input)
			step_size = this.alpha;
			
			best_obj = 999999
			best_step = nothing;
			for i = 1:10
				obj = problem_data.objective_function(x + step_size*this.dx)
				
				if obj < best_obj
					best_obj = obj;
					best_step = step_size;
				end
				
				step_size = step_size*0.9;
			end
			
			this.alpha = best_step;
		end
		
		return this
	end
end