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

#endif
