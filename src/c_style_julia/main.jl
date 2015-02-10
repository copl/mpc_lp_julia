include("interior_point_algorithm.jl")

function main() # this function is called at the bottom of the code
	max_iter             = 20;  Total number of iterarions
    linear_feas_tol      = 1e-8;  # Threshold for feasibility of the linear constrains.
    comp_tol            = 1e-8;  # Threshold for complementry optimality condition s^Tz
    bkscale              = 0.95;  # Back scaling for line search.

   # Initialize configuration variable
    settings = class_settings(max_iter,linear_feas_tol,comp_tol,bkscale);

	srand(1234)
 # We are creating an instance of LP. In practice, we should read the problem data from input stream.
	problem_data = construct_instance1()
        # The main function that run interior point algorithm.
	interior_point_algorithm(problem_data,settings)
end

# Constructing a random LP
function construct_instance1()
	n = 5;
	k = 10;
	m = 10;
	
	problem_data = class_linear_program_input()
	x0 = rand(k,1)

	A = rand(n, k);
	G = -diagm(ones(m));
	c = rand(k,1)
	h = zeros(m);
	b = A*x0;


	problem_data.A = A
	problem_data.G = G
	problem_data.c = c
	problem_data.h = h
	problem_data.b = b
	problem_data.m = m
	problem_data.k = k
	problem_data.n = n
	
	return(problem_data)
end

######################
#  run the program   #
######################
main()
