# Functions
# A,G, c, h, b, x, s, z, y, tau, kappa, q1, q2, q3, q4, q5, q6
# Linear Equation Solver for steps (1) Affine direction (2) Corrector direction.
# Input to the function
#       problem_data::class_linear_program_input: Input data for the LP that includes A, G, h, b, c, m, n, k
#       variables::class_linear_program_variables: Variables (primal and dual) of LP which includes x, y, s, z, tau, kappa
#       rhs::class_linear_system_rhs: Right hand side for the Linear Equation. Right hand side is differnet for
#               Affine direction and  Corrector direction
#       K_newton_matrix::class_K_newton_matrix: This data structure is used to share information when we find
#               Affine direction and Corrector direction. In the algorithm we first find Affine direction, and next Corrector direction.
#               Since, the factorization calculated in  Affine direction can be used in finding Corrector direction,
#               data is shared between these two function calls.
#       
# Output: The function modified LP variables in variables::class_linear_program_variables

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
	
	q7 = rhs.q3 - rhs.q5./z;
    q8 = rhs.q4 - rhs.q6./tau;
	
	new_rhs = -[q1; q2; q7; q8];
	
	result = K_newton_matrix.F\new_rhs
	
    return(result)
end