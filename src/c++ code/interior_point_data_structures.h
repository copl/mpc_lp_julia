#ifndef INTERIOR_POINT_DATA_STRUCTURES
#define INTERIOR_POINT_DATA_STRUCTURES

#include "matrix_data_structures.h"

//class_settings

class class_settings {
	int max_iter;
	float linear_feas_tol;	//This is a relative tolerance w.r.t. 
                            //some normalizing norms
	float comp_tol; // How small must s^Tz must be when we stop
	
	//Constant length of the 
    //maximum combined step to the boundary to use
    float bkscale;
	
	//Configuration for solver
	//linear_solver_settings 
	public:
		class_settings(int max_iter,	float linear_feas_tol,	float comp_tol,	float bkscale);
};

class_settings::class_settings (int input_max_iter, float input_linear_feas_tol, float input_comp_tol, float input_bkscale) {
	max_iter 			= input_max_iter;
	linear_feas_tol 	= input_linear_feas_tol;
	comp_tol 			= input_comp_tol;
	bkscale 			=input_bkscale;
}
//--------End class_settings--------

//class_direction
class class_direction {
	class_vector dx;
	class_vector dy;
	class_vector dz;
	class_vector ds;
	float dtau;
	float dkappa;
	float alpha;
public:
	class_direction();
	void update_values(class_vector dx, class_vector dy, class_vector dz, class_vector ds, float dtau, float dkappa, float alpha);
	void compute_affine_direction(	
		class_linear_system_rhs affine_rhs,
		class_linear_program_input problem_data,
		class_linear_program_variables variables,
		class_K_newton_matrix K_newton_matrix
		);
	void compute_corrector_direction(
		class_linear_system_rhs corrector_rhs,
		class_linear_program_input problem_data,
		class_linear_program_variables variables,
		class_algorithm_state state,
		class_settings settings,
		class_K_newton_matrix K_newton_matrix
		);
	void compute_alpha(
		class_algorithm_state state,
		class_settings settings
		);
	void compute_min_ratio_alpha (
		class_vector var;
		class_vector dvar;
		);
	float get_alpha();
	float get_dtau(); 
	float get_dkappa(); 
	class_vector get_dx(); 
	class_vector get_dy(); 
	class_vector get_dz(); 
	class_vector get_ds(); 

};
//--------End class_direction--------

//class_residuals

class class_residuals {
	class_vector r1;
	class_vector r2;
	class_vector r3;
	class_vector r4;

	float r1_norm;
	float r2_norm;
	float r3_norm;
	float normed_squared;
public:
	class_residuals();
	void update_values();
	void compute_residuals();
};
//--------End class_residuals--------

//class_linear_program_input
class class_linear_program_input {
	class_matrix A;
	class_matrix G;
	
	class_vector c;
	class_vector h;
	class_vector b;
	
	int m;
	int n;
	int k;
public:
	class_linear_program_input();
};
//--------End class_linear_program_input--------

//class_K_newton_matrix
class class_K_newton_matrix {
	public:
		class_K_newton_matrix();

};

//--------End class_K_newton_matrix--------

//class_linear_program_result
class class_linear_program_result {

};
//--------End class_linear_program_result--------

//class_linear_program_variables
class class_linear_program_variables {
	class_vector x;
	class_vector s;
	class_vector z;
	class_vector y;
	float tau;
	float kappa;
public:
	class_linear_program_variables(class_linear_program_input problem_data);
	void take_step(class_direction direction);
};
//--------End class_linear_program_variables--------

//class_algorithm_state
class class_algorithm_state {
private:
	float mu;
	float sigma;
	float gap;

public:
	class_algorithm_state();
	void update_mu (class_linear_program_variables variables, class_linear_program_input progblem_data); //TODO
	void update_gap (class_linear_program_variables variables, class_linear_program_input problem_data); //TODO
};


//--------End class_algorithm_state--------



//class_linear_system_rhs

class class_linear_system_rhs {
	class_vector q1;
	class_vector q2;
	class_vector q3;
	class_vector q4;
	class_vector q5;
	class_vector q6;
	void update_values(class_vector q1, class_vector q2, class_vector q3, class_vector q4, class_vector q5, class_vector q6); //Not called in algorithm
public:
	class_linear_system_rhs(class_linear_program_input problem_data);
	
	void compute_affine_rhs(class_residuals residuals, class_linear_program_variables variables);
	void compute_corrector_rhs(class_residuals residuals, class_linear_program_variables variables);

};

//--------End class_linear_system_rhs--------




#endif
