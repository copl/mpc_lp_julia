# Functions
# A,G, c, h, b, x, s, z, y, tau, kappa, q1, q2, q3, q4, q5, q6
function solveLinearEquation(problem_data::class_linear_program_input,variables::class_linear_program_variables,rhs::class_linear_system_rhs,K_newton_matrix::class_K_newton_matrix)
	A = problem_data.A
	G = problem_data.G
	h = problem_data.h
	b = problem_data.b
	c = problem_data.c
	m = problem_data.m
	n = problem_data.n
	k = problem_data.k
	
	kappa = variables.kappa
	tau = variables.tau
	s = variables.s
	z = variables.z
	
	q1 = rhs.q1
	q2 = rhs.q2
	q3 = rhs.q3
	q4 = rhs.q4
	q5 = rhs.q5
	q6 = rhs.q6
	
	
	Q =  [ zeros(k,k)  A'          G'               c;
     	      -A           zeros(n,n)  zeros(n, m)      b;
	      -G           zeros(m,n)  diagm(vec(s./z))  h;
	      -c'          -b'         -h'              kappa/tau];

	q7 = q3 - q5./z;
        q8 = q4 - q6./tau;

        -Q\[q1; q2; q7; q8];
end