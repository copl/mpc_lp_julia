function sparse_condition_number(M)
	return sqrt(abs(eigs(M'*M; nev=1, which=:LM, ritzvec=false)[1][1] / eigs(M'*M; nev=1, which=:SM, ritzvec=false)[1][1]));
end

function ldl(A) # temporary until julia gets hooked up to correct libraries	
	(n,m) =size(A)
	assert(n==m)
	d=zeros(n)
	L=spzeros(n,n)
	
	for j=1:n
		  L[j,j]=1;
		  S=A[j,j];
		  for k=1:j-1
				 S=S-d[k]*L[j,k]^2;
		  end
		  
		  d[j]=S
		  for i=(j+1):n
			  L[j,i]=0;
			  S=A[i,j];
			  for k=1:j-1
					 S=S-d[k]*L[i,k]*L[j,k];
			  end
			  L[i,j]=S/d[j]
		  end
	end

	D = spdiagm(d)
	
	return(L,D)
end


#Settings for the algorithm
type class_settings
    
    max_iter::Int   
	
	primal_feas_tol::Real
	dual_feas_tol::Real
	duality_gap_tol::Real
	primal_infeas_tol::Real
	dual_infeas_tol::Real
	
    #Constant length of the 
    #maximum combined step to the boundary to use
    bkscale::Real

    function class_settings()
        this = new()
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
				println("black box constraints and adding constraints is dangerous!")
				Error()
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
					#println("z norm " * string(norm(zNormalized,1)))
					#println("x/tau norm " * string(norm(xNormalized,1)))
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
				@assert isdefined(this,:P)
				@assert isdefined(this,:A)
				@assert isdefined(this,:G)
				@assert isdefined(this,:c)
				@assert isdefined(this,:b)
				@assert isdefined(this,:h)
				@assert isdefined(this,:k)
				@assert isdefined(this,:n)
				@assert isdefined(this,:m)
				@assert isdefined(this,:barrier_weights)
				
				@assert (size(this.P) == (this.k,this.k))
				
				@assert size(this.A) == (this.n,this.k)
					
				@assert size(this.G) == (this.m,this.k)
				
				@assert size(this.c) == (this.k,)
				
				@assert size(this.b) == (this.n,)
				
				@assert size(this.h) == (this.m,)
					
				@assert size(this.barrier_weights) == (this.m,)
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
	K_bar::SparseMatrixCSC{Float64,Int64}
	calc_inertia::Function
	
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
			xNorm = variables.x/variables.tau
			
			D_g = variables.kappa/variables.tau #+  0.5*xNorm'*problem_data.L*xNorm;
			
			@time this.Q = sparse([ problem_data.L   problem_data.A'          problem_data.G'    sparse(problem_data.c);
				 problem_data.A           spzeros(n,n)  spzeros(n, m)      -sparse(problem_data.b);
				 problem_data.G           spzeros(m,n)  -spdiagm((variables.s)./variables.z)  -sparse(problem_data.h);
				sparse(-problem_data.c')         sparse(-problem_data.b')        sparse(-problem_data.h')      sparse([variables.kappa/variables.tau])
				]);
				
					
			
			#condition_number = cond(full(M))
			#condition_number = sparse_condition_number(M);
			#println("Condition number " * string(condition_number))
			
			@time this.F = lufact(this.Q);
			
			if false
				this.K_bar = 	sparse([ problem_data.L   problem_data.A'  problem_data.G';
					 problem_data.A           spzeros(n,n)  spzeros(n, m)      ;
					 problem_data.G           spzeros(m,n)  -spdiagm((variables.s)./variables.z) ;
					]);	
					
				L, D = ldl(this.K_bar);
				println("nz",length(nonzeros(L)))
				inertia = this.calc_inertia(diag(D));
				#expected_inertia = (n+m,l+m,0) pg 5. (3.6) Knitro: An Integrated Package for Nonlinear Optimization
				# k + n, 
				println((k,m,n))
				expected_inertia = (k,n+m,0)
				println(inertia)
				println(expected_inertia)
			end
		end
				
		this.calc_inertia = function(perturbed_eigens)
			# TO DO
			#print("lu factor correct ")
			#print(norm(F[:L]*diagm(1./diag(F[:L]))-F[:U]'*diagm(1./diag(F[:U]))))
			#print("\n")
			
			pos_eigs = 0
			neg_eigs = 0
			zero_eigs = 0
			for val in perturbed_eigens
				if val > 0
					pos_eigs = pos_eigs + 1
				elseif val < 0
					neg_eigs = neg_eigs + 1
				else
					zero_eigs = zero_eigs + 1
				end
			end
			
			return pos_eigs, neg_eigs, zero_eigs
		end
		
		this.solve = function(problem_data::class_linear_program_input,variables::class_linear_program_variables,rhs::class_linear_system_rhs)
			q7 = rhs.q3 - rhs.q5./variables.z;
			q8 = rhs.q4 - rhs.q6./variables.tau;
						
			result = -(this.F\[rhs.q1; -rhs.q2; -q7; q8])
			
			return(result)
		end
		
		this.solve2 = function(eta,gamma,problem_data::class_linear_program_input,variables::class_linear_program_variables,residuals)
			r_hat = [residuals.r_dual; 
			-residuals.r_primal + diag(s)gamma*mu*problem_data.barrier_weights];
			
			
			#
			#
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
	line_search::Function
	check::Function
	
	function class_linear_program_variables(problem_data::class_linear_program_input)
		this = new();
		
		this.s = ones(problem_data.m);
		this.z = ones(problem_data.m);# + rand(problem_data.m);
		this.x = ones(problem_data.k); #zeros(problem_data.k)*
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
		
		this.line_search = function(direction,residuals,problem_data)
			variables = this;
			
			orginal_obj = residuals.summary_norm;
			
			variables.take_step(direction)
			residuals.compute_residuals(problem_data,variables)
			best_obj = residuals.summary_norm;
			last_step = direction.alpha;
			best_step = last_step;
			next_step = last_step;
			
			if true
				for i = 1:5
					try
						next_step = last_step*0.6
						direction.alpha = next_step - last_step
						variables.take_step(direction)
						variables.check()
						problem_data.update(variables)
						residuals.compute_residuals(problem_data,variables)
						
						if residuals.summary_norm < best_obj
							best_obj = residuals.summary_norm;
							best_step = next_step;
						end
						last_step = next_step;
					catch e
						println("error " * string(i))
						throw(e)
					end
				end
			end
			
			direction.alpha = -next_step;
			variables.take_step(direction)
			direction.alpha = best_step;
			
			if orginal_obj <= best_obj
				println("error, residuals norm did not decrease")
			end
			
			println("best step" * string(best_step))
			println("best value" * string(best_obj))
			
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
				rho = mu_a/mu;
				#sigma = ((mu_a/(mu))^3)
				sigma = max(0,(min(1,rho)))^3
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
	r1::Array{Float64,1} # vector
	r2::Array{Float64,1} # vector
	r3::Array{Float64,1} # vector
	r4 #::Array{Float64,1} # real
	
	r_primal::Array{Float64,1}
	r_dual::Array{Float64,1}
	r_gap::Real
	
	
	r1_norm::Real
	r2_norm::Real
	r3_norm::Real
	r4_norm::Real
	
	
	homogeneous_dual_objective_value::Real
	homogeneous_primal_objective_value::Real
	duality_gap::Real
	
	dual_feasibility::Real
	primal_infeasibility::Real
	primal_feasibility::Real
	dual_infeasibility::Real
	
	summary_norm::Real
	
	update_norms::Function
	compute_residuals::Function
	
	function class_residuals(problem_data::class_linear_program_input)
		this = new();
		
		this.compute_residuals = function(pd::class_linear_program_input, variables::class_linear_program_variables)
			try
				# grad = pd.P*(variables.x/variables.tau) + problem_data.c;
				# this.r1 =  grad*variables.tau
				
				this.r1 = -pd.A'*variables.y - pd.G'*variables.z - pd.grad*variables.tau;
				this.r2 = pd.A*variables.x - pd.b*variables.tau;
				this.r3 = variables.s + (pd.G*variables.x) - variables.tau*pd.h; 
				this.r4 = variables.kappa + pd.grad'*variables.x + pd.b'*variables.y + pd.h'*variables.z;
				
				this.r_primal = [this.r2; this.r3];
				this.r_dual = this.r1;
				this.r_gap = this.r4[1];
				
				# 5.4 pg 395 homogenous; A homogeneous algorithm for MCP Anderson and Ye
				# I added kappa to denominator r1, r2 and r3
				normalization_factor = 1.0/(variables.tau + variables.kappa)
				this.r1_norm = normalization_factor*norm(this.r1,1)#/(norm(variables.y,1) + norm(variables.z,1));
				this.r2_norm = normalization_factor*norm(this.r2,1)#/(norm(variables.x,1));
				this.r3_norm = normalization_factor*norm(this.r3,1)#/(norm(variables.s,1) + norm(variables.x,1));
				this.r4_norm = normalization_factor*norm(this.r4,1)#/(norm(variables.x,1) + norm(variables.y,1) + norm(variables.z,1));
				
				# this.r1_norm = this.r1_norm
				
				xNormalized = variables.x/variables.tau
				yNormalized = variables.y/variables.tau
				zNormalized = variables.z/variables.tau
				sNormalized = variables.s/variables.tau
				
				this.homogeneous_dual_objective_value = -(pd.b'*yNormalized + pd.h'*zNormalized)[1];
				this.homogeneous_primal_objective_value = -(pd.grad'*xNormalized)[1];
				this.duality_gap = (sNormalized'*zNormalized)[1]
				
				this.dual_feasibility = norm(this.r1,1)/variables.tau;
				this.primal_feasibility = norm(this.r2/variables.tau,1) + norm(this.r3/variables.tau,1);
				
				this.primal_infeasibility = norm(pd.A'*yNormalized + pd.G'*zNormalized,1)/this.homogeneous_dual_objective_value;
				this.dual_infeasibility = (norm(pd.A*xNormalized,1) + norm(sNormalized + (pd.G*xNormalized),1))/this.homogeneous_primal_objective_value;
				
				#println(this.dual_feasibility)
				#println(this.primal_feasibility)
				#println(this.dual_feasibility)
				#println(scaled_duality_gap2)
				#println("duality gap: " * string(this.duality_gap))
				#println("primal feasibility: " * string(this.primal_feasibility))
				#println("dual feasibility: " * string(this.dual_feasibility))
				#println("primal infeasibility: " * string(this.primal_infeasibility))
				#println("duality infeasibility: " * string(this.dual_infeasibility))
				#println("tau: " * string(variables.tau))
				#println("kappa: " * string(variables.kappa))
				#println(homogeneous_dual_objective_value)
				#println(homogeneous_primal_objective_value)
				
				this.summary_norm = this.r1_norm + this.r2_norm + this.r3_norm + this.r4_norm;
			
				#println(this.r1_norm," ",this.r2_norm," ",this.r3_norm," ",this.r4_norm)
			
				debug_message("r1")
				debug_message(this.r1)
				debug_message("r2")
				debug_message(this.r2)
				debug_message("r3")
				debug_message(this.r3)
				debug_message("r4")
				debug_message(this.r4)
				
			catch e
				println("*** ERROR in compute_residuals ***")
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
		
		return this
	end
end