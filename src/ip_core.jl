# TO DO:
# - stopping criteron and infeasibility detection
# - get inequality stuff going**
# - get split stuff going
# - test feasible non-linear constraints more throughly
# - get infeasible non-linear constraints working - more test examples

include("ip_local_approximation.jl")
include("ip_direction.jl")


println("loading ip_core")

GLOBAL_P_NORM = 2;

type class_settings
    
    max_iter::Int64  
	max_iter_line_search::Int64
	
	duality_gap_tol::Float64
	primal_feas_tol::Float64
	dual_feas_tol::Float64
	
	primal_infeas_tol::Float64
	dual_infeas_tol::Float64
	
	# see page 260
	beta1::Float64
	beta2::Float64
	beta3::Float64
	beta4::Float64
	beta5::Float64
	beta6::Float64
	
	solution_status::Function
	print_status::Function
	
    #Constant length of the 
    #maximum combined step to the boundary to use
    #bkscale::Float64

    function class_settings()
        this = new();
		
		this.solution_status = function(state::class_state)
			# criteron on resiudals ??????
			
			# status:
			# 0 = unknown
			# 1 = optimal
			# 2 = primal infeasible
			# 3 = dual infeasible
			
			if state.relative_gap < this.duality_gap_tol && state.primal_feasibility < this.primal_feas_tol && state.dual_feasibility < this.dual_feas_tol
				return 1 # optimal solution found
			end
			
			if state.homogeneous_dual_sign == 1 && state.homogeneous_dual_feasibility < this.primal_infeas_tol
				return 2 # primal infeasible
			end
			
			if state.homogeneous_primal_sign == 1 && state.homogeneous_primal_feasibility < this.dual_infeas_tol
				return 3 # dual infeasible
			end
			
			# status unknown
			return 0
		end
		
		this.print_status = function(status::Int64)
			if status == 0
				println("Unknown solutions status")
			elseif status == 1
				println("Optimal solution found")
			elseif status == 2
				println("Primal infeasible")
			elseif	status == 3
				println("Dual infeasible")
			else
				println("ERROR class_state.class_settings")
				ERROR("")
			end
		end
		
        return this;
    end
end

type class_non_linear_program
	###########################################################################
	# user defined 
	###########################################################################
	
	# variable dimensions
	n_1::Int64 				# number of primal variables with x >= 0
	n_2::Int64 				# number of unbounded primal variables
	n::Int64 				# total number of primal variables
	
	m_1::Int64 				# number of dual variables y >= 0 (same as the number of inequality constraints)
	m_2::Int64 				# number of unbounded dual variables
	m::Int64				# total number of dual variables
	
	# the objective and constraints of the non-linear program
	
	# objective_function(x_scaled)
	# c(x,x_bar)
	objective_function::Function						
	
	# objective_function_gradient(x_scaled)
	# D(x,x_bar) c(x,x_bar)
	objective_function_gradient::Function				
	
	# objective_function_hessian(x_scaled)
	# D(x,x_bar)^2 c(x,x_bar)
	objective_function_hessian::Function				
	
	# evaluate_constraints(x_scaled)
	# evaluates a(x,x_bar)
	evaluate_constraints::Function						
	
	# evaluate_constraint_gradients(x_scaled)
	# evaluates D(x,x_bar) a(x,x_bar)
	evaluate_constraint_gradients::Function				
	
	# evaluate_constraint_lagrangian_hessian(x_scaled,y_scaled)
	# evaluates D(x,x_bar)^2 (y,y_bar)' * a(x,x_bar)
	evaluate_constraint_lagrangian_hessian::Function	
	
	# weights on log barrier function. standard is all ones.
	barrier_weights::Array{Float64,1}
	
	###########################################################################
	# NOT user defined
	###########################################################################
	
	# for simple setup i.e. a quadratic program
	set_linear_objective::Function
	set_quadratic_objective::Function
	set_linear_constraints::Function
	
	# validate function: calls validate_dimensions and validate_derivatives
	validate::Function
	
	# check the dimensions are correct of the nlp, variables and direction vector
	validate_dimensions::Function
	
	# check the derivatives are correct
	validate_derivatives::Function
	
	# check variables satisfy bounds
	check_bounds_are_satisfied::Function
	
	function class_non_linear_program()
        this = new()
		
		this.set_linear_objective = function(c::Array{Float64,1})
			this.set_quadratic_objective(c,spzeros(length(c),length(c)))
		end
		
		this.set_quadratic_objective = function(c::Array{Float64,1},Q::SparseMatrixCSC{Float64,Int64})
			# min c'*x_scaled + 0.5*x_scaled'*Q*x_scaled
			this.objective_function = function(x_scaled)
				return ((c ' * x_scaled) + 0.5*(x_scaled)'*Q*(x_scaled))[1]
			end
			this.objective_function_gradient = function(x_scaled)
				return c + Q*(x_scaled)
			end
			this.objective_function_hessian = function(x_scaled)
				return Q
			end
		end
		
		this.set_linear_constraints = function(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, A_bar::SparseMatrixCSC{Float64,Int64}, b_bar::Array{Float64,1}, n_1::Int64 )
			# creates a problem with the following constraints:
			# A*x_scaled - b >= 0
			# A_bar*x_scaled - b_bar == 0
			
			# check dimensions
			try
				this.n_1 = n_1;
				this.m_1, this.n = size(A)
				this.m_2, this.n = size(A_bar)
				this.m = this.m_1 + this.m_2
				
				this.n_2 = this.n - this.n_1;
				
				# check everything is ok.
				@assert(this.n_1 <= this.n)
				@assert(this.n_2 <= this.n)
				@assert(this.m_1 <= this.m)
				@assert(this.m_2 <= this.m)
				
				@assert(this.m_1 + this.m_2 == this.m)
				@assert(this.n_1 + this.n_2 == this.n)
				@assert (size(A) == (this.m_1,this.n))
				@assert (size(b) == (this.m_1,))
				@assert (size(b_bar) == (this.m_2,))
			
				this.evaluate_constraints = function(x_scaled)
					return [A*x_scaled - b; A_bar*x_scaled - b_bar];
				end
				
				this.evaluate_constraint_gradients = function(x_scaled)
					return [A; A_bar];
				end
				
				this.evaluate_constraint_lagrangian_hessian = function(x_scaled,y_scaled)
					return spzeros(length(x_scaled),length(x_scaled))
				end
				
			catch e
				println("ERROR class_non_linear_program.set_linear_constraints")
				throw(e)
			end
		end
		
        return this
    end
end

type class_variables 
	# see page 251
	x::Array{Float64,1} 		# bounded primal variables
	x_bar::Array{Float64,1} 	# unbounded primal variables
	y::Array{Float64,1} 		# dual variables
	y_bar::Array{Float64,1} 	# unbounded dual variables
	z::Array{Float64,1} 		# primal slack
	s::Array{Float64,1} 		# dual slack
	tau::Float64				# tau (large tau => feasible)
	kappa::Float64				# kappa (large kappa => infeasible)
	
	x_scaled::Array{Float64,1}
	y_scaled::Array{Float64,1}
	
	take_step::Function			# 
	check::Function
	
	function class_variables(nlp::class_non_linear_program)
		class_variables(nlp.n_1,nlp.n_2,nlp.m_1,nlp.m_2)
	end
	
	function class_variables(n_1::Int64,n_2::Int64,m_1::Int64,m_2::Int64)
		this = new();
		
		this.x = ones(n_1);
		this.x_bar = zeros(n_2);
		this.y = ones(m_1);
		this.y_bar = zeros(m_2);
		this.s = ones(n_1);
		this.z = ones(m_1);
		this.tau = 1.0;
		this.kappa = 1.0;
		
		this.x_scaled = [this.x; this.x_bar]/this.tau;
		this.y_scaled = [this.y; this.y_bar]/this.tau;
		
		this.take_step = function(direction::class_direction)
			try
				alpha = direction.alpha; # get the step size
				this.x = this.x + alpha * direction.dx;
				this.x_bar = this.x_bar + alpha * direction.dx_bar;
				this.y = this.y + alpha * direction.dy;
				this.y_bar = this.y_bar + alpha * direction.dy_bar;
				this.s = this.s + alpha * direction.ds;
				this.z = this.z + alpha * direction.dz;
				
				this.tau = this.tau + alpha * direction.dtau;
				this.kappa = this.kappa + alpha * direction.dkappa;
				
				this.x_scaled = [this.x; this.x_bar]/this.tau;
				this.y_scaled = [this.y; this.y_bar]/this.tau;
				
				this.check()
			catch e
				println("error inside class_variables.take_step")
				throw(e)
			end
		end
		
		this.check = function()
			try
				for x_i in this.x
					@assert(x_i > 0)
				end
				
				for z_i in this.z
					@assert(z_i > 0)
				end
				
				for y_i in this.y
					@assert(y_i > 0)
				end
				
				for s_i in this.s
					@assert(s_i > 0)
				end
				
				@assert(this.tau > 0)
				@assert(this.kappa > 0)
				
				@assert(norm(this.x_scaled - [this.x; this.x_bar]/this.tau,1) < 10^(-8.0));
				@assert(norm(this.y_scaled - [this.y; this.y_bar]/this.tau,1) < 10^(-8.0));
				
			catch e
				println("one variable is less than or equal to zero")
				throw(e)
			end
		end
		
		#this.copy = function()
		#	copy_vars = class_variables()
		#	return copy_vars;
		#end
		
		return(this)
	end
end



function randomly_generate_variables(nlp)
	vars = class_variables(nlp)
	vars.x = rand(vars.n_1 + vars.n_2)
end

function validate(nlp,vars)
	# derivative checker
	
	validate_dimensions(nlp,vars) 
	vars.check()
end

function validate_dimensions(nlp::class_non_linear_program,vars::class_variables) 
	# this function checks that dimensions of the non-linear program and variables match
	
	try
		@assert(nlp.n_1 == length(vars.x))
		@assert(nlp.n_2 == length(vars.x_bar))
		@assert(nlp.m_1 == length(vars.y))
		@assert(nlp.m_2 == length(vars.y_bar))
		@assert(nlp.n == nlp.n_1 + nlp.n_2)
		@assert(nlp.m == nlp.m_1 + nlp.m_2)
		@assert(nlp.n == length(vars.x_scaled))
		@assert(nlp.m == length(vars.y_scaled))
		
		@assert(nlp.m_1 == length(vars.z))
		@assert(nlp.n_1 == length(vars.s))
		
		
		try
			obj = nlp.objective_function(vars.x_scaled);
			obj::Float64
		catch e
			println("ERROR objective_function")
			throw(e)
		end
		
		try
			c = nlp.objective_function_gradient(vars.x_scaled);
			c::Array{Float64,1}
			@assert((nlp.n,) == size(c))
		catch e
			println("ERROR objective_function_gradient")
			throw(e)
		end
		
		try
			H = nlp.objective_function_hessian(vars.x_scaled)
			H::SparseMatrixCSC{Float64,Int64}
			@assert((nlp.n,nlp.n) == size(H))
		catch e
			println("ERROR objective_function_hessian")
			throw(e)
		end
		
		try
			a = nlp.evaluate_constraints(vars.x_scaled);
			a::Array{Float64,1}
			@assert((nlp.m,) == size(a))
		catch e
			println("ERROR evaluate_constraints")
			throw(e)
		end
		
		try
			J = nlp.evaluate_constraint_gradients(vars.x_scaled);
			J::SparseMatrixCSC{Float64,Int64}
			@assert((nlp.m,nlp.n) == size(J))
		catch e
			println("ERROR evaluate_constraint_gradients")
			throw(e)
		end
		
		try
			L = nlp.evaluate_constraint_lagrangian_hessian(vars.x_scaled,vars.y_scaled);
			L::SparseMatrixCSC{Float64,Int64}
			@assert((nlp.n,nlp.n) == size(L))
		catch e
			println("ERROR evaluate_constraint_lagrangian_hessian")
			throw(e)
		end
		
		
		
	catch e
		println("ERROR validate_dimensions")
		throw(e)
	end
end



function validate_derivatives(nlp::class_non_linear_program,vars::class_variables) 
	# this function uses finite differences to check that the hessian etc, are being correctly computed
end


