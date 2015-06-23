EMPTY_ARRAY = spzeros(0,1)*[0.0]


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

function standard_settings()
	settings = class_settings();
	
	settings.max_iter = 70;  # Total number of iterarions
	settings.max_iter_line_search = 40;
	
	settings.primal_feas_tol = 1e-8
	settings.dual_feas_tol = 1e-8
	settings.duality_gap_tol = 1e-10
	settings.primal_infeas_tol = 1e-8
	settings.dual_infeas_tol = 1e-8
	
	settings.beta1 = 10.0^(-4)
	settings.beta2 = 10.0^(-8)
	settings.beta3 = 10.0^(-4)
	settings.beta4 = 0.9
	settings.beta5 = 0.999
	settings.beta6 = 0.5
	
	return settings
end

function get_netlib_problem(file_name)
	file = matopen(file_name)
	A= sparse(read(file, "A"))
	b = read(file,"b")
	b = b[:]
	c = read(file,"c")
	c = c[:];
	
	return A,b,c
end