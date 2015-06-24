function proximal_step(x_hat, sigma)
	# solves sigma*|x-x_hat|_2^2 - log(x)
	# 2*sigma*(x-x_hat) - 1/x = 0
	# 2*sigma*x^2 - 2*sigma*x*x_hat - 1 = 0
	# x = (sigma*x_hat + sqrt(sigma^2*x_hat^2 + 2*sigma))/(2*sigma)
	
	return (sigma*x_hat + sqrt(sigma^2*x_hat^2 + 2*sigma))/(2*sigma)
end

function evaluate_objective(A, b, x_hat, pos_indicies, mu)
	return norm(A*x_hat - b,1) - mu * sum( log(x_hat[pos_indicies] ))
end

function proximal_gradient_step(A,b,x_hat,pos_indicies,mu,t_max,beta)
	# takes a proximal gradient step on |A[x, x_bar] - b|_2^2 - mu*log(x)W
	
	gradient = A'*sign(A*x_hat - b)
	true_gradient = gradient
	true_gradient[pos_indicies] = gradient[pos_indicies] - mu./x_hat[pos_indicies]
	
	
	intial_value = evaluate_objective(A, b, x_hat, pos_indicies, mu)
	println("Residuals:",intial_value)
	
	x = zeros(length(x_hat))
	
	t = t_max
	for j = 1:60
		#x = x_hat - t*true_gradient
		x = x_hat - t*gradient
		
		for i = pos_indicies
			x[i] = proximal_step(x[i], (1/(2*t*mu)))
			x[i] = minimum([maximum([x[i],0.01]),100])
		end
		
		if evaluate_objective(A, b, x, pos_indicies, mu) < intial_value - 1e-2 * t * norm(true_gradient,2)^2
			println("improvement")
			break
		end
		
		t *= beta
	end
	println("step size = ",t)
	
	return x
end


	
function run_proximal_solver(A_bar,b_bar,c_bar)
	mu = 1.0
	t_max = 1e-2;
	beta = 0.7;

	n,k = size(A_bar)

	KKT_system = sparse([[A_bar spzeros(n,k) spzeros(n,n)]
							 [spzeros(k,k) speye(k,k) A_bar' ]
							 [c_bar' spzeros(1,k) -b_bar']]);
	KKT_rhs = [b_bar; c_bar; 0];

	point = ones(2*k + n);
	pos_indicies = [1:2*k]


	for i = 1:20
		#mu = (point[1:k]' * point[(k+1):(2*k)])[1]/(k)
		#println("mu ",mu)
		point = proximal_gradient_step(KKT_system,KKT_rhs,point,pos_indicies,mu,t_max,beta)
	end
end

	
function run_proximal_solver_homogenous(A_bar,b_bar,c_bar)
	mu = 1e-1
	t_max = 1e-2;
	beta = 0.7;

	n,k = size(A_bar)
							# x					tau				kappa				s					y
	KKT_system = sparse([[A_bar			 -b_bar 		spzeros(n,1)		spzeros(n,k) 		spzeros(n,n)]
							 [spzeros(k,k)	 c_bar 		spzeros(k,1)		speye(k,k) 		-A_bar' ]
							 [-c_bar' 			 1				1					spzeros(1,k) 		b_bar']]);

	KKT_rhs = zeros(k + n + 1);

	point = [ones(2*(k+1)); zeros(n)]
	pos_indicies = [1:2*(k+1)]


	for i = 1:100
		mu = (point[1:k]' * point[(k+1):(2*k)])[1]/(k)
		println("mu ",mu)
		point = proximal_gradient_step(KKT_system,KKT_rhs,point,pos_indicies,mu,t_max,beta)
	end
	
	return point
end

using MAT

file = matopen("problems/blend.mat")

A = sparse(read(file, "A"))
n,k = size(A)	

b = read(file,"b")
b = b[:];

c = read(file,"c")
c = c[:];

println("A matrix is size " * string((n,k)) * " and has " * string(nnz(A))  * " non-zeros")
point = run_proximal_solver_homogenous(A,b,c)

include("ip_algorithm.jl")
include("testing_tools.jl")


lp = class_non_linear_program();
lp.set_linear_objective(c);
lp.set_linear_constraints(spzeros(0,k),spzeros(0,1)*[1.0],A,b,k);

vars = class_variables(lp);

if true
	vars.x = point[1:k]
	vars.tau = point[k+1]
	vars.kappa = point[k+2]
	vars.s = point[(k+3):(2*k+2)]
	vars.y_bar = point[(2*k+3):(2*k+2 + n)]
	vars.x_scaled = vars.x/vars.tau
	vars.y_scaled = vars.y_bar/vars.tau
end


vars, status, iter = ip_algorithm(lp, standard_settings(), vars, true);



