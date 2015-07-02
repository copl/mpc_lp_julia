#eliminate all duplicate rows
using JuMP
using MAT


function check()
  include("testing_tools.jl")
  println("loaded libraries")


  A,b,c = get_netlib_problem("small_problems/AFIRO.mat")
  # min c'*x
  # A x = b
  # x >= 0

  f = open("values.txt", "w")
  #=for i = 1:size(A,2)
    writecsv(f, A[:,i])
    write(f, "\n")
  end=#
  writecsv(f,A)
  writecsv(f,c)
  close(f)
  S1 = LPsolve(A,b,c)
  S2 = trivialcheck(A,b,c)
  if isequal(S1,S2) == false
    if (S1-S2)/S2 < 0.001
      return "ok"
    end
    return false
  else
    return true
  end

end
function checks(A,b,c)
  #all checks go here

end

function LPsolve(A,b,c)
  tic()
  m = Model()

  @defVar(m, x[1:size(A,2)] >= 0)

  for i = 1:size(A,1)
    @addConstraint(m, sum{A[i,j]*x[j], j = 1:size(A,2)} == b[i])
  end
  @setObjective(m, Min, sum{c[i]*x[i], i = 1:size(A,2)})

  print(m)
  status = solve(m)
  println("Status: ", status)
  println("Objective value: ", getObjectiveValue(m))
  toc()
  return getObjectiveValue(m)
end
function check(A,b,c,S)

  m = Model()
  for j = 1:size(A,2)
    cnt = 0
    for i = 1:size(A,1)
      if A[i,j] != 0
        cnt += 1
      end
    end
    println(j, " ", cnt, " ", c[j])
  end
  if size(A,1) == 0 || size(A,2) == 0
    return S
  end

  @defVar(m, x[1:size(A,2)] >= 0)


  for i = 1:size(A,1)
    @addConstraint(m, sum{A[i,j]*x[j], j = 1:size(A,2)} == b[i])
  end

  @setObjective(m, Min, S + sum{c[i]*x[i], i = 1:size(A,2)})

  print(m)
  status = solve(m)
  println("Status: ", status)
  println("Objective value: ", getObjectiveValue(m))

  return getObjectiveValue(m)
end
function trivialcheck(A, b, c, l, u) #if 0 or 1 values in row/column

  tic()
  row = falses(size(A,1))
  col = falses(size(A,2))
  A1 = Array(Float64, size(A,1), size(A,2))
  b1 = Array(Float64, size(A,1))
  c1 = Array(Float64, size(c))
  S = 0

  convert(Float64, A[1,1])
  convert(Float64, b[1])
  convert(Float64, c[1])

  for i = 1:size(A,2)
    c1[i] = c[i]
  end

  for i = 1:size(A,1)
    b1[i] = b[i]
    for j = 1:size(A,2)
      A1[i,j] = A[i,j]
    end
  end

  for i = 1:size(A,1)
	  cnt = 0
    k = 0
	  for j = 1:size(A,2)
      A1[i,j] = A[i,j]
		  if A[i,j] != 0
		  	cnt += 1
        if k == 0
          k = j
        end
      end
	  end
  	if cnt == 0
		  if b[i] == 0
		  	row[i] = true #remove row
	  	else
			  println("infeasible")
			  return 2
		  end
  	elseif cnt == 1 #row singleton
  		row[i] = true
      col[k] = true
      x = b[i]/A[i,k]
      S += c[k]*x
      #c1[k] = 0
      #A1[i,k] = 0
		  #substitute into other equations
  		for j = 1:size(A,1)
	  		if j != i
		  		b1[j] -= x*A[j,k]
			  end
		  end
	  end
  end


  for i = 1:size(A,2)
    if col[i] == false
      cnt = 0
      k = 0 #marker
      for j = 1:size(A,1)
        if A[j,i] != 0
          cnt += 1
          if k == 0
            k = j
          end
        end
      end
      if cnt == 0 #empty column
        col[i] = true
        if c[i] >= 0
          c[i] = 0 #equivalent to x[i] = 0
        else
          println("unbounded")
          return 3
        end
      elseif cnt == 1 #column singleton - substitution

        row[k] = true
        col[i] = true
        for j = 1:size(A,2)
          if col[j] == false
            c1[j] -= c[i]*A[k,j]/A[k,i] #check
          end
        end
        S += c[i]*b[k]/A[k,i]
      end
    end
  end

 #= for i = 1:size(A,1)
    if row[i] == true
      println(i)
    end
  end
  for i = 1:size(A,2)
    if col[i] == true
      println(i)
    end
  end=#

  n = size(A,1)-1

  while n > 1
    if row[n] == true
        A1 = A1[[1:n-1,n+1:end],:]
        b1 = b1[[1:n-1,n+1:end]]
    end
    n -= 1
  end

  if size(b1) == 0
    return S
  elseif size(b1) != 1
    if row[end] == true
      A1 = A1[1:end-1,:]
      b1 = b1[1:end-1]
    end
    if row[1] == true
      A1 = A1[2:end,:]
      b1 = b1[2:end]
    end
  end

  n = size(A,2)-1

  while n > 1
    if col[n] == true
        A1 = A1[:,[1:n-1,n+1:end]]
        c1 = c1[[1:n-1,n+1:end]]
    end
    n -= 1
  end

  if size(c1) == 0
    return S
  elseif size(c1) != 1
    if col[end] == true
      A1 = A1[:,1:end-1]
      c1 = c1[1:end-1]
    end
    if col[1] == true
      A1 = A1[:,2:end]
      c1 = c1[2:end]
    end
  end

  println("The new matrix A:")
	println(A1)
  println("Matrix B:")
  println(b1)
  println("Matrix C:")
  println(c1)

  print("Old size: ")
  print(size(A,1), size(A,2))
  print("New size: ")
  print(size(A1,1), size(A1,2))
  for j = 1:size(A,2)
    cnt = 0
    for i = 1:size(A,1)
      if A[i,j] != 0
        cnt += 1
      end
    end
      println(j, " ", cnt, " ", col[j])

  end
  sol = check(A1,b1,c1,S)
  toc()
  return sol
end

function rowcheck(A, b)
row = falses(size(b))


#normalization
for i = 1:size(A,1)
	if A[i,1] != 0
		b[i] /= A[i,1]
		for j = 2:size(A,2)
			A[i,j] /= A[i,1]
		end
		A[i,1] = 1
	end
end

#sort here

for i = 1:size(A,1)
	for j = i+1:size(A,1)
		k=1
		while (A[i, k] == A[j, k] && k <= size(A,2))
			if k == size(A,2)
				if b[i] == b[j]
					row[j] = true
				else
					println("infeasible")
					return

				end
			end
			k += 1

		end

	end
end
k = 1
A1 = zeros(size(A,1), size(A,2))
b1 = zeros(size(b))
for i = 1:size(A,1)

	if row[i] != true
		b1[k] = b[i]
		for j = 1:size(A,2)
			A1[k,j] = A[i,j]
		end
		k += 1
	end
end
k -= 1

for i = 1:k
	#=for j = 1:size(A,2)
		print(A1[i,j])
		print(" ")
	end=#
	print(A1[i,:])
	println()
end
println()
for i = 1:k
	print(b1[i])
	print(" ")
end

end


function columncheck(A,c)

col = falses(size(c))

for i = 1:size(A,2)
	for j = i+1:size(A,2)
		k=1
		while (A[k, i] == A[k, j] && k <= size(A,1))
			if k == size(A,1)
				if c[i]+c[j]<0 && c[i]
					println("unbounded")
					return
				elseif c[i] > c[j]
					col[i] = true
				elseif c[j] > c[i]
					col[j] = true
				end
				break

			end
			k += 1

		end
	end
end

k = 1
A1 = zeros(size(A,1), size(A,2))
c1 = zeros(size(c))
for i = 1:size(A,2)

	if col[i] != 1
		c1[k] = c[i]
		for j = 1:size(A,1)
			A1[j,k] = A[j,i]
		end
		k += 1
	end
end
k -= 1

for i = 1:size(A,1)
	for j = 1:k
		print(A1[i,j])
		print(" ")
	end
	println()
end
println()
for i = 1:k
	print(c1[i])
	print(" ")
end

end

function subset(A, b, c, n)
x = zeros(size(A,2))

for i = 1:size(A,1)
	cnt = 0
	for j = 1:size(A,2)
		if A[i,j] != 0
			cnt += 1
		end
	end
	if cnt <= n
		for j = 1:size(A,1)
			B = true
			if j != i
				for k = 1:size(A,2)
					if A[i,k] != 0 && A[j,k] == 0
						B = false
					end
				end
			end
			if B == true
				for k = 1:size(A,2)
					if A[i,k] != 0
            x[k] = A[j,k]/A[i,k] #optimize?
          end
				end
			end
		end
	end

end


end
