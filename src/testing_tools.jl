
function is_problem_successful(problem_name,status,correct_value,details=None)
	if status == correct_value
		println(problem_name, " sucessfully solved")
	else
		println("---------------------------------------")
		println(problem_name, " NOT SOLVED SUCESSFULLY.")
		println("Status = ", status, " instead of ", correct_value)
		
		if details != None
			println("Details:")
			println(details) # i.e. problem description
		end	
		println("---------------------------------------")
	end
	
	return (status == correct_value)
end