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

#Settings for the linear solver
type linear_solver_settings
#Empty for now
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
	# get_m
	# get_n
	# get_k

	function class_linear_program_input()
		this = new ();

		return this
	end
end

type class_linear_program_result
	# empty for now
	result # integer codes: optimal_solution_found, unbounded, infeasible
	x
end

type class_linear_program_state
	x
	s
	z
	y
	tau
	kappa

	function class_linear_program_state(problem_data)
	  this = new();
	  this.x = zeros(problem_data.k, 1);
	  this.s = ones(problem_data.m, 1);
	  this.z = ones(problem_data.m, 1);
	  this.y = zeros(problem_data.n, 1);
	  this.tau = 1;
	  this.kappa = 1;

	  return this
	end
end

type class_linear_system_rhs
	q1
	q2
	q3
	q4
	q5
	q6
	
	function class_linear_system_rhs(q1,q2,q3,q4,q5,q6)
		this = new()
		this.q1 = q1
		this.q2 = q2
		this.q3 = q3
		this.q4 = q4
		this.q5 = q5
		this.q6 = q6
		
		return this
	end
end

type class_residuals
	r1
	r2
	r3
	r4
	function class_residuals(r1,r2,r3,r4)
		this = new()
		this.r1 = r1
		this.r2 = r2
		this.r3 = r3
		this.r4 = r4
		
		return this
	end
end

type class_direction
	dx
	dy
	dz
	dtau
	ds
	dkappa
	alpha
	
	function class_direction(dx,dy,dz,dtau,ds,dkappa,alpha)
		this = new();
		this.dx = dx;
		this.dy = dy;
		this.dz = dz;
		this.dtau = dtau;
		this.ds = ds
		this.dkappa = dkappa;
		this.alpha = alpha;
		
		return this
	end
end