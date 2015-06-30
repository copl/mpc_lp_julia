using JuMP

function solve(A,b,c)
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
  for i = 1:size(A,2)
    println(getValue(x[i]))
  end
  toc()
end