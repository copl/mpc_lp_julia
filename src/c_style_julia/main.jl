include("interior_point_algorithm.jl")

function main() # this function is called at the bottom of the code
	max_iter = 20;
	tol = 1e-8;

	srand(1234)
	problem_data = construct_instance1()
	interior_point_algorithm(problem_data,max_iter,tol)
end

function construct_instance1()
	n = 5;
	k = 10;
	m = 10;
	
	problem_data = linear_program_input()
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
	
	return problem_data
end
######################
#  run the program   #
######################
main()
