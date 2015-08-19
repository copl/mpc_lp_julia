using MUMPS
using MAT

# Load a UF problem
file="/home/ubuntu/COPL/mpc_lp_julia/src/OLD2/example_problems/ex3sta1/ex3sta1.mat"
vars = matread(file)
A = vars["Problem"]["A"]
b = vars["Problem"]["b"]
c = vars["Problem"]["aux"]["c"]
m=size(A)[2]
n=size(A)[1]
k=m
G=sparse(1:m,1:m,-1)

# Create Q from eq 17 of algorithm
#    m     n     m
# m   [0    A'    G'
# n   A     0     0
# m   G     0     D]
Q = spzeros(m+n+m, m+n+m)
for i = 1:n
    for j = 1:m
        Q[j, m+i]=A[i,j]
        Q[m+i, j] = A[i,j]
    end
end
for i = 1:m
    for j = 1:m
        Q[j, m+n+i] = G[i,j]
        Q[m+n+i, j] = G[i,j]
    end
end

for i=1:m
    Q[m+n+i,m+n+i] = rand(1)[1]
end

# create right hand side
h=zeros(m)
rhs = [c;b;h]
x = solveMUMPS(Q,rhs);
        #y = -Q\[q1; q2; q7; q8];
        A = sparse(-Q)
        rhs = [q1; q2; q7; q8];
        x = solveMUMPS(A,rhs);
 
