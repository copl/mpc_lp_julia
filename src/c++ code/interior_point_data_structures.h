#ifndef INTERIOR_POINT_DATA_STRUCTURES
#define INTERIOR_POINT_DATA_STRUCTURES

#include "matrix_data_structures.h"


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

void class_settings::class_settings (int max_iter, float linear_feas_tol, float comp_tol, float bkscale) {
	this->max_iter = max_iter;
	this->linear_feas_tol = linear_feas_tol;
	this->comp_tol = comp_tol;
	this-> bkscale=bkscale;
}
class class_linear_program_input {
	public:
		class_matrix A;
		class_matrix G;
		
		class_vector c;
		class_vector h;
		class_vector b;
		
		int m;
		int n;
		int k;
		class_linear_program_input() {};
};

class class_K_newton_matrix {

};
class class_linear_program_result {

};

class class_algorithm_state {
	float mu;
	float sigma;
	float gap;

public:
	class_linear_program_state();
	update_mu (); //TODO
	update_gap (); //TODO
};
class class_linear_program_variables {
	class_vector x;
	class_vector s;
	class_vector z;
	class_vector y;
	float tau;
	float kappa;
public:
	class_linear_program_variables();
	take_step();
};

class class_linear_system_rhs {
	class_vector q1;
	class_vector q1;
	class_vector q3;
	class_vector q4;
	class_vector q5;
	class_vector q6;
public:
	class_linear_system_rhs();
	update_values();
	compute_affine_rhs();
	compute_corrector_rhs();

};

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
	update_values();
	compute_residuals();
};

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
	update_values();
	compute_affine_direction();
	compute_corrector_direction();
};

#endif
