#include "definitions.h"

// Problem parameters.
double H_INIT = -15.0;
double H_ELEVATION = 10.0;
double K_S_vals[4] = {350.2, 712.8, 1.68, 18.64};
double ALPHA_vals[4] = {0.01, 1.0, 0.01, 0.01};
double N_vals[4] = {2.5, 2.0, 1.23, 2.5};
double M_vals[4] = {0.864, 0.626, 0.187, 0.864};

double THETA_R_vals[4] = {0.064, 0.0, 0.089, 0.064};
double THETA_S_vals[4] = {0.14, 0.43, 0.43, 0.24};
double STORATIVITY_vals[4] = {0.1, 0.1, 0.1, 0.1};


const int MATERIAL_COUNT = 4;
const int CONSTITUTIVE_TABLE_METHOD = 2;
						  
const int NUM_OF_INTERVALS = 16;
const double INTERVALS_4_APPROX[16] = {-1.0, -2.0, -3.0, -4.0, -5.0, -8.0, -10.0, -12.0,
				     -15.0, -20.0, -30.0, -50.0, -75.0, -100.0,-300.0, -1000.0};

double TABLE_LIMIT = -1000.0;
const double TABLE_PRECISION = 0.1;
const double LOW_LIMIT=-1.0;
const int NUM_OF_INSIDE_PTS = 0;

bool POLYNOMIALS_READY = false;
bool POLYNOMIALS_ALLOCATED = false;