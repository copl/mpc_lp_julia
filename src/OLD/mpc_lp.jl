using Base.SparseMatrix

#Generate a random feasible problem 
A = sprand(10,100,0.5);
b = A*rand(100,1);
c = A'*randn(10,1) + rand(100,1);


#Basic mehrotra predictor corrector in julia
#-------------------------------------------

#Set algorithm parameters 
max_iter = 20
abs_pre  = 1.e-5
feas_res = 1.e-7
feas_tol = 1.e-7

#Define the function that evaluates the gradient and hessian
#Evaluates the hessian at the N-T scaling point and the gradient at x
evalHessianAndGradient(x,s) = (spdiagm((s./x,),0),-1./x)


#Determine the problem size 
(m,n) = size(A)
nu    = n; #For the lp case

#Intialize the vectors
x     = ones(n);
s     = ones(n);
y     = zeros(m);
t     = 1;
k     = 1;
mu    = (s'*x+t*k)/(nu+1);
mu    = mu[1]

#Allocate space for the work vectors
pr    = zeros(m);
dr    = zeros(n);
gr    = 0;
iter  = 0;

for iter =1:max_iter
    #Evaluate the Hessian and gradient 
    (H,g) = evalHessianAndGradient(x,s)
    #Evaluate the residuals 
    dr =  A'*y-t*c+s;
    pr = -(A*x-t*b);
    gr = c'*x-b'*y;
    gr = gr[1] + k;
     
    #Assemble the matrix 
    K = [H   -A'          c;
         A  spzeros(m,m) -b;
         -c'    b'       k./t]
    #Affine rhs
    rhs = [dr-s;pr;gr-k]
    
    #Brute force solving using backslash
    dir = K\rhs
    dx  = dir[1:n];
    dy  = dir[n+1:n+m];
    dt  = dir[n+m+1];
    ds  = -s-H*dx;
    dk  = -k-k/t*dt;
    
    # Find the largest steps to the boundary x+adx = 0 for a >0 iff dx<0. In that case a = x./dx
    a = 1.0;
    for i=1:n
        if(dx[i]<0)
            a = min(-x[i]/dx[i],a)
        end
        if(ds[i]<0)
            a = min(-s[i]/ds[i],a)
        end
        if(dt <0)
            a = min(-t/dt,a)
        end
        if(dk <0)
            a = min(-k/dk,a)
        end
    end

    sigma = (1-a)^3;

    #Find the combined direction 
    rhs = [dr-s-sigma*mu*g-(1-sigma)^2*(dx.*ds)./x;pr;gr-k+sigma*mu/t-(1-sigma)^2*(dt*dk)/t]
    dir = K\rhs
    dx  = dir[1:n];
    dy  = dir[n+1:n+m];
    dt  = dir[n+m+1];
    ds  = -s-sigma*mu*g-H*dx-(1-sigma)^2*(dx.*ds)./x;
    dk  = -k+sigma*mu*1./t-k/t*dt-(1-sigma)^2*(dt*dk)/t;
    
    #Find the largest step to the boundary
    a = 1.0;
    for i=1:n
        if(dx[i]<0)
            a = min(-x[i]/dx[i],a)
        end
        if(ds[i]<0)
            a = min(-s[i]/ds[i],a)
        end
        if(dt <0)
            a = min(-t/dt,a)
        end
        if(dk <0)
            a = min(-k/dk,a)
        end
    end
    #Backtrack away from the boundary 
    a = a*0.95;
    #Update the variables
    x = x+a*dx;
    y = y+a*dy;
    s = s+a*ds;
    t = t+a*dt;
    k = k+a*dk;
    mu = x'*s
    mu = mu[1]+t*k;
    mu = mu/(nu+1)
    
    #Calculate the residual norm 
    npr = Base.norm(pr)
    ndr = Base.norm(dr)
    ngr = Base.norm(gr)
    #println("Iter: $iter, npr: $npr, ndr: $ndr, ngr: $ngr, sigma: $sigma, alpha: $a")

    @printf("Iter: %3i, npr: %3.3e, ndr: %3.3e, ngr: %3.3e, sigma: %3.3e, alpha: %3.3e\n",iter,npr,ndr,ngr,sigma,a)
    #Evaluate the termination criteria
    if(npr<feas_tol&&ndr<feas_tol&&ngr<feas_tol) 
        println("Required precision achieved, terminating.")
        break;
    end
 end


