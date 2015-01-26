include("interior_point_data_structures.jl")
include("linear_system_solver.jl")


function interior_point_algorithm(problem_data::class_linear_program_input,	settings::class_settings)

	println("Iter \t gap     \t mu   \t alpha     \t tau       \t feasibility")

	state = class_linear_program_state(problem_data) # intialize state of linear program solver with default values
	
	#initial value of mu
    mu = ((state.s)'*(state.z) + (state.tau)*(state.kappa))/(problem_data.m+1)
    mu = mu[1];
	
	for itr =1:settings.max_iter
		residuals = compute_residuals(problem_data,state);
		
		#Evaluate termination criteria	
        #TODO: This part will have to be a more sofisticated test to detect 
        #unbounded and infeasible problems.
		if (Base.norm(residuals.r1) < settings.linear_feas_tol && 
            Base.norm(residuals.r2) < settings.linear_feas_tol &&
            Base.norm(residuals.r3) < settings.linear_feas_tol && 
            mu < settings.comp_tol)
			 println("Ended");
			 break;
		end
		
			affine_direction = compute_affine_direction(problem_data,state,residuals);
		
	    corrector_direction = compute_corrector_direction(problem_data,
                                                          state,
                                                          residuals,
                                                          affine_direction,
                                                          settings,
                                                          mu);
		
		alpha = corrector_direction.alpha; # get the step size
        
		#Take the step
		state.x = state.x + alpha*corrector_direction.dx;
		state.s = state.s + alpha*corrector_direction.ds;
		state.y = state.y + alpha*corrector_direction.dy;
		state.z = state.z + alpha*corrector_direction.dz;
		state.tau = state.tau + alpha*corrector_direction.dtau;
		state.kappa = state.kappa + alpha*corrector_direction.dkappa;
	    
		
		#Compute the gap after the step
		gap = ((problem_data.c)'*(state.x) + (problem_data.h)'*(state.z) + (problem_data.b)'*(state.y))[1];
		mu = ((state.s)'*(state.z) + (state.tau)*(state.kappa))/(problem_data.m+1);
		mu = mu[1];
		
		@printf("%3i\t%3.3e\t%3.3e\t%3.3e\t%3.3e\t%3.3e\n",itr,gap,mu,alpha, state.tau, mean(abs(norm_squared(residuals))))
	end
end

function norm_squared(residual::class_residuals)
	r = [residual.r1, residual.r2, residual.r3, residual.r4]
	return r.*r
end

function compute_residuals(pd::class_linear_program_input, state::class_linear_program_state)

	r = [zeros(1,pd.k) zeros(1,pd.n) (state.s)'  state.kappa]' -
			[ zeros(pd.k,pd.k)  pd.A'          pd.G'           pd.c;
			 -pd.A           zeros(pd.n,pd.n)  zeros(pd.n,pd.m)  pd.b;
			 -pd.G           zeros(pd.m,pd.n)  zeros(pd.m,pd.m)  pd.h;
			 -pd.c'          -pd.b'         -pd.h'          0]*
			 [state.x; state.y; state.z; state.tau];

	r1 = r[1:pd.k];
	r2 = r[(pd.k+1):(pd.k+pd.n)];
	r3 = r[(pd.k+pd.n+1):(pd.k+pd.n+pd.m)];
	r4 = r[pd.k+pd.n+pd.m+1];
	
	return(class_residuals(r1,r2,r3,r4))
end

function compute_affine_direction(problem_data::class_linear_program_input,	state::class_linear_program_state,	residuals::class_residuals)

	affine_rhs = class_linear_system_rhs(-residuals.r1, -residuals.r2, -residuals.r3, -residuals.r4, -(state.z).*(state.s), -(state.tau)*(state.kappa))
	dir = solveLinearEquation(problem_data, state, affine_rhs)
	A = problem_data.A
	G = problem_data.G
	h = problem_data.h
	b = problem_data.b
	c = problem_data.c
	m = problem_data.m
	n = problem_data.n
	k = problem_data.k
	
	x = state.x
	z = state.z
	s = state.s
	y = state.y
	tau = state.tau
	kappa = state.kappa
	
	dx_a = dir[1:k];
	dy_a = dir[(k+1):(k+n)];
	dz_a = dir[(k+n+1):(k+n+m)];
	dtau_a = dir[(k+n+m+1)];
	ds_a = ( -z.*s - dz_a.*s)./z;
	dkappa_a = (-(tau)*(kappa) - dtau_a*(kappa))/(tau)

	# Compute Step size α, and Centering Parameter σ
	 vv = vec([s./ds_a; z./dz_a; tau/dtau_a; kappa/dkappa_a]);
    alpha = 1;
    for i=1:length(vv)
        if (vv[i] < 0)
            alpha = minimum([alpha, -vv[i]]);
        end
    end
	
	return(class_direction(dx_a,dy_a,dz_a,dtau_a,ds_a,dkappa_a,alpha))
end


function compute_corrector_direction(problem_data::class_linear_program_input,
                                     state::class_linear_program_state,
                                     residuals::class_residuals,
                                     affine_direction::class_direction,
                                     settings::class_settings,
                                     mu::Real)
	A = problem_data.A
	G = problem_data.G
	h = problem_data.h
	b = problem_data.b
	c = problem_data.c
	m = problem_data.m
	n = problem_data.n
	k = problem_data.k
	
	kappa = state.kappa
	tau = state.tau
	s = state.s
	z = state.z
	
	dx_a = affine_direction.dx
	dy_a = affine_direction.dy
	ds_a = affine_direction.ds
	dz_a = affine_direction.dz
	dtau_a = affine_direction.dtau
	dkappa_a = affine_direction.dkappa
	alpha = affine_direction.alpha
	
	mu_a = ((s+alpha*ds_a)'*(z+alpha*dz_a) + (tau + alpha*dtau_a)*((kappa) + alpha*dkappa_a))/(m+1);
	sigma = ((mu_a/(mu))^3)[1]
	
	corrector_rhs = class_linear_system_rhs(-(1-sigma)*residuals.r1, -(1-sigma)*residuals.r2, -(1-sigma)*residuals.r3, -(1-sigma)*residuals.r4, -z.*s -ds_a.*dz_a + sigma*mu,  -tau*kappa-dtau_a*dkappa_a + sigma*mu)
	dir = solveLinearEquation(problem_data, state, corrector_rhs);
	
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
	
	return(class_direction(dx,dy,dz,dtau,ds,dkappa,alpha))
end
