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
	void update_values();
	void compute_affine_direction();
	void compute_corrector_direction();
	float get_alpha() {return alpha;}
	float get_dtau() {return dtau;}
	float get_dkappa() {return dkappa;}
	class_vector get_dx() {return dx;}
	class_vector get_dy() {return dy;}
	class_vector get_dz() {return dz;}
	class_vector get_ds() {return ds;}

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
private:
	class_matrix A;
	class_matrix G;
	
	class_vector c;
	class_vector h;
	class_vector b;
	
	int m;
	int n;
	int k;
public:
	class_linear_program_input() {return; }
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
class_linear_program_variables::class_linear_program_variables(class_linear_program_input problem_data){

	tau = 1;
	kappa = 1;
}

void class_linear_program_variables::take_step(class_direction direction){
	float alpha = direction.get_alpha();
	//TODO: Implement class_vector.multiply and class_vector.add
	// x = x.add(direction.get_dx.multiply(alpha));
	// s = s.add(direction.get_ds.multiply(alpha));
	// z = z.add(direction.get_dz.multiply(alpha));
	// y = y.add(direction.get_dy.multiply(alpha));
	tau = tau + alpha * direction.get_dtau();
	kappa = kappa + alpha * direction.get_dkappa();

}
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

class_algorithm_state::class_algorithm_state() {
	
}
void class_algorithm_state::update_mu(class_linear_program_variables variables, class_linear_program_input problem_data){
	
}

void class_algorithm_state::update_gap(class_linear_program_variables variables, class_linear_program_input problem_data){

}
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
class_linear_system_rhs::class_linear_system_rhs(class_linear_program_input problem_data){

}
void class_linear_system_rhs::update_values(class_vector q1, class_vector q2, class_vector q3, class_vector q4, class_vector q5, class_vector q6){

}

//--------End class_linear_system_rhs--------




#endif
