max_iter = 20;
tol = 1e-8;

srand(1234)
# Problem 1
n = 5;
k = 10;
m = 10;

x0 = rand(k,1)
A = rand(n, k);
G = -diagm(ones(m));
c = rand(k,1)
h = zeros(m);
b = A*x0;


# Functions
function solveLinearEquation(A,G, c, h, b, x, s, 
	z, y, tau, kappa, q1, q2, q3, q4, q5, q6)
  
	Q =  [ zeros(k,k)  A'          G'               c;
     	      -A           zeros(n,n)  zeros(n, m)      b;
	      -G           zeros(m,n)  diagm(vec(s./z))  h;
	      -c'          -b'         -h'              kappa/tau];	

	q7 = q3 - q5./z;
        q8 = q4 - q6./tau;
       
        -Q\[q1; q2; q7; q8];
end

# Initialization

x = zeros(k, 1);
s = ones(m, 1);
z = ones(m, 1);
y = zeros(n, 1);
tau = 1;
kappa = 1;

println("Iter\tgap     \tmu      \tsigma   \talpha     \ttau       \tfeasibility")

for itr =1:max_iter
    # Computing Residuals

    r = [zeros(1,k) zeros(1,n) s' kappa]' - 
        [ zeros(k,k)  A'          G'           c;
         -A           zeros(n,n)  zeros(n, m)  b;
         -G           zeros(m,n)  zeros(m, m)  h;
         -c'          -b'         -h'          0]*
         [x; y; z; tau];
    r1 = r[1:k];
    r2 = r[(k+1):(k+n)];
    r3 = r[(k+n+1):(k+n+m)];
    r4 = r[k+n+m+1];

	#println(r1,r2,r3,r4)
	#break
	
    if (Base.norm(r1) < tol && Base.norm(r2) < tol && Base.norm(r3) < tol)
         println("Ended");
         break;
    end

    mu = (s'*z + tau*kappa)/(m+1)
    mu = mu[1];
    dir = solveLinearEquation(A,G, c, h, b, x, s,
        z, y, tau, kappa, -r1, -r2, -r3, -r4, -z.*s, -tau*kappa)

    # affine direction
    dx_a = dir[1:k];
    dy_a = dir[(k+1):(k+n)];
    dz_a = dir[(k+n+1):(k+n+m)];
    dtau_a = dir[(k+n+m+1)];
    ds_a = ( -z.*s - dz_a.*s)./z;
    dkappa_a = (-tau*kappa - dtau_a*kappa)/tau

    # Compute Step size α, and Centering Parameter σ
    vv = vec([s./ds_a; z./dz_a; tau/dtau_a; kappa/dkappa_a]);
    alpha = 1;
    for i=1:length(vv)
        if (vv[i] < 0)
            alpha = minimum([alpha, -vv[i]]);
        end
    end
	
	
    mu_a = ((s+alpha*ds_a)'*(z+alpha*dz_a) + (tau + alpha*dtau_a)*(kappa + alpha*dkappa_a))/(m+1);
    sigma = ((mu_a/mu)^3)[1]

	#println(mu_a)
	
	#println(dx_a,dy_a,dz_a)
	#break
	
    # Corrector direction


    dir = solveLinearEquation(A,G, c, h, b, x, s,
         z, y, tau, kappa, -(1-sigma)*r1, -(1-sigma)*r2, -(1-sigma)*r3,
	 -(1-sigma)*r4, -z.*s -ds_a.*dz_a + sigma*mu, 
	 -tau*kappa-dtau_a*dkappa_a + sigma*mu);

	#println(dir) 
	 
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
    alpha = alpha*.99
	
	#println(dx)
	#println(alpha)
	#break
	
    x = x + alpha*dx;
    s = s + alpha*ds;
    y = y + alpha*dy;
    z = z + alpha*dz;
    tau = tau + alpha*dtau;
    kappa = kappa + alpha*dkappa;

    gap = (c'*x + h'*z + b'*y)[1];
    @printf("%3i\t%3.3e\t%3.3e\t%3.3e\t%3.3e\t%3.3e\t%3.3e\n",itr,gap,mu,sigma,alpha, tau, mean(abs(r.*r)))
    #println("$itr\t$gap\t$mu\t$sigma\t$alpha")

end
