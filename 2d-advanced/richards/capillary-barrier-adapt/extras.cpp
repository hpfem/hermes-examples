#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace std;

//Debugging matrix printer.
bool printmatrix(double** A, int n, int m){
  for (int i=0; i<n; i++){
    for (int j=0; j<m; j++){
      printf(" %lf ", A[i][j]) ;
    }
    printf(" \n");
  }
  printf("----------------------------------\n");
  return true;
}

//Debugging vector printer.
bool printvector(double* vect, int n){
  for (int i=0; i<n; i++){
    printf(" %lf ", vect[i]);
  }
  printf("\n");
  printf("----------------------------------\n");
  return true;
}

// Creates a table of precalculated constitutive functions.
bool get_constitutive_tables(int method, ConstitutiveRelationsGenuchtenWithLayer* constitutive, int material_count)
{
  info("Creating tables of constitutive functions (complicated real exponent relations).");

  // Table values dimension.
  int bound = int(-constitutive->table_limit/constitutive->table_precision)+1;
  
  // Allocating arrays. 
  constitutive->k_table = new double*[material_count];
  for (int i=0; i<material_count; i++) {
    constitutive->k_table[i] = new double[bound];
  }
  
  constitutive->dKdh_table = new double*[material_count] ;
  for (int i=0; i<material_count; i++) {
    constitutive->dKdh_table[i] = new double[bound];
  }
  
  constitutive->dKdh_table = new double*[material_count] ;
  for (int i=0; i<material_count; i++) {
    constitutive->dKdh_table[i] = new double[bound];
  }

  constitutive->c_table = new double*[material_count] ;
  for (int i=0; i<material_count; i++) {
    constitutive->c_table[i] = new double[bound];
  }
  
  //If Newton method (method==1) selected constitutive function derivations are required.
  if (method==1){
    constitutive->dCdh_table = new double*[material_count] ;
    for (int i=0; i<material_count; i++) {
      constitutive->dCdh_table[i] = new double[bound];
    }
    
    constitutive->ddKdhh_table = new double*[material_count] ;
    for (int i=0; i<material_count; i++) {
      constitutive->ddKdhh_table[i] = new double[bound];
    }
  }
  
  // Calculate and save K(h).
  info("Calculating and saving K(h).");
  for (int j=0; j<material_count; j++) {
    for (int i=0; i< bound; i++) {
      constitutive->k_table[j][i] = constitutive->K(-constitutive->table_precision*i, j);
    }
  }
  // Calculate and save dKdh(h).
  info("Calculating and saving dKdh(h).");
  for (int j=0; j<material_count; j++) {
    for (int i=0; i< bound; i++) {
      constitutive->dKdh_table[j][i] = constitutive->dKdh(-constitutive->table_precision*i, j);
    }
  }
  // Calculate and save C(h).
  info("Calculating and saving C(h).");
  for (int j=0; j<material_count; j++) {
    for (int i=0; i< bound; i++) {
      constitutive->c_table[j][i] = constitutive->C(-constitutive->table_precision*i, j);
    }
  }
  
  
  //If Newton method (method==1) selected constitutive function derivations are required.
  if (method==1){
    // Calculate and save ddKdhh(h).
    info("Calculating and saving ddKdhh(h).");
    for (int j=0; j<material_count; j++) {
      for (int i=0; i< bound; i++) {
	constitutive->ddKdhh_table[j][i] = constitutive->ddKdhh(-constitutive->table_precision*i, j);
      }
    }
    // Calculate and save dCdh(h).
    info("Calculating and saving dCdh(h).");
    for (int j=0; j<material_count; j++) {
      for (int i=0; i< bound; i++) {
	constitutive->dCdh_table[j][i] = constitutive->dCdh(-constitutive->table_precision*i, j);
      }
    }
  }	
      
  return true;
}

// Simple Gaussian elimination for full matrices called from init_polynomials().
bool gem_full(double** A, double* b, double* X, int n){
  int i,j,k;
  
  double** aa;
  double dotproduct, tmp;
  aa = new double*[n];
  
  for (i=0; i<n; i++){
    aa[i] = new double[n+1];
  }
   
  for (i=0;i<n; i++){
    for (j=0; j<n; j++){
      aa[i][j] = A[i][j] ;
    }
    aa[i][n] = b[i];
  }

  for (j=0; j<(n-1); j++){
    for (i=j+1; i<n; i++){
    tmp = aa[i][j]/aa[j][j];

      for (k=0; k<(n+1); k++){
	aa[i][k] = aa[i][k] - tmp*aa[j][k] ;
      }
    }
  }

  for (i=n-1; i>-1; i--){
    dotproduct=0.0;
    for (j=i+1; j<n; j++){
      dotproduct = dotproduct + aa[i][j]*X[j] ;
    }
    X[i] = (aa[i][n]-dotproduct)/aa[i][i] ;
  }
    
  delete []aa;
  return true;
}

// Initialize polynomial approximation of constitutive relations close to full saturation for constitutive->constitutive_table_method=1.
// For constitutive->constitutive_table_method=2 all constitutive functions are approximated by polynomials, K(h) function by quintic spline, C(h) 
// function by cubic splines. Discretization is managed by variable int num_of_intervals and double* intervals_4_approx.
// ------------------------------------------------------------------------------
// For constitutive->constitutive_table_method=1 this function requires folowing arguments:
// n - degree of polynomials
// low_limit - start point of the polynomial approximation
// points - array of points inside the interval bounded by <low_limit, 0> to improve the accuracy, at least one is recommended. 
// An approriate amount of points related to the polynomial degree should be selected.
// n_inside_point - number of inside points
// layer - material to be considered.
//------------------------------------------------------------------------------
// For constitutive->constitutive_table_method=2, all parameters are obtained from global definitions.
bool init_polynomials(int n, double low_limit, double *points, int n_inside_points, int layer, ConstitutiveRelationsGenuchtenWithLayer* constitutive, int material_count, int num_of_intervals, double* intervals_4_approx)
{
  double** Aside;
  double* Bside;
  double* X;
  switch (constitutive->constitutive_table_method) 
    {
      // no approximation 
      case 0 :
	break ;
      // polynomial approximation only for the the K(h) function surroundings close to zero
      case 1 :
	
	if (constitutive->polynomials_allocated == false){

	  constitutive->polynomials = new double**[material_count];
	  
	  for (int i=0; i<material_count; i++){
	    constitutive->polynomials[i] = new double*[3] ;
	  }
	  
	  for (int i=0; i<material_count; i++){
	    for (int j=0; j<3; j++){
	      constitutive->polynomials[i][j] = new double[n_inside_points+6] ;
	    }
	  }
	  constitutive->polynomials_allocated = true;
	}

	Aside = new double*[n+n_inside_points] ;
	Bside = new double[n+n_inside_points] ;
	for (int i=0; i<n; i++){
	  Aside[i] = new double[n+n_inside_points] ;
	}

	// Evaluate the first three rows of the matrix (zero, first and second derivative at point low_limit).
	for (int i=0; i<n; i++){
	  for (int j=0; j<n; j++){
	    Aside[i][j] = 0.0;
	  }
	}
	for (int i=0; i<n; i++){
	  Aside[3][i] = pow(low_limit,i) ;
	  Aside[4][i] = i*pow(low_limit, i-1) ;
	  Aside[5][i] = i*(i-1)*pow(low_limit, i-2) ;
	}
	Bside[3] = constitutive->K(low_limit, layer) ;
	Bside[4] = constitutive->dKdh(low_limit, layer) ;
	Bside[5] = constitutive->ddKdhh(low_limit, layer) ; 
	
	// Evaluate the second three rows of the matrix (zero, first and second derivative at point zero).
	Aside[0][0] = 1.0 ;

	// For the both first and second derivative it does not really matter what value is placed there.
	Aside[1][1] = 1.0 ;
	Aside[2][2] = 2.0 ;
      
	Bside[0] = constitutive->K(0.0, layer) ;
	Bside[1] = 0.0;
	Bside[2] = 0.0;
      
	for (int i=6; i<(6+n_inside_points); i++){
	  for (int j=0; j<n; j++) {
	    Aside[i][j]=pow(points[i-6],j) ;
	  }
	  printf("poradi, %i %lf %lf \n", i, constitutive->K(points[i-6], layer), points[i-6]);

	  Bside[i] = constitutive->K(points[i-6], layer) ;
	  printf("layer, %i \n", layer);
	}

	gem_full(Aside, Bside, constitutive->polynomials[layer][0], (n_inside_points+6));
	
	for (int i=1; i<3; i++){
	  for (int j=0; j< (n_inside_points+5); j++){
	    constitutive->polynomials[layer][i][j] = (j+1)*constitutive->polynomials[layer][i-1][j+1];
	  }
	  constitutive->polynomials[layer][i][n_inside_points+6-i] = 0.0;
	}
	
	
	delete [] Aside;
	delete [] Bside;
	break ;
      // polynomial approximation for all functions at interval (table_limit, 0)
      case 2 :
	
	int pts = 0;
	
	if (constitutive->polynomials_allocated == false) {
	  // K(h) function is approximated by quintic spline.
	  constitutive->k_pols = new double***[num_of_intervals];
	  //C(h) function is approximated by cubic spline.
	  constitutive->c_pols = new double***[num_of_intervals];
	  
	  for (int i=0; i<num_of_intervals; i++){
	    constitutive->k_pols[i] = new double**[material_count];
	    constitutive->c_pols[i] = new double**[material_count];
	    for (int j=0; j<material_count; j++){
	      constitutive->k_pols[i][j] = new double*[3];
	      constitutive->c_pols[i][j] = new double*[2];
	      for (int k=0; k<3; k++) {
		constitutive->k_pols[i][j][k] = new double[6 + pts];
		if (k<2)
		  constitutive->c_pols[i][j][k] = new double[5];
	      }
	    }
	  }
	  
// 	  //allocate constitutive->pol_search_help array -- an index array with locations for particular pressure head functions
	  constitutive->pol_search_help = new int[int(-constitutive->table_limit)+1];
	  
	  for (int i=0; i<int(-constitutive->table_limit); i++){
	    for (int j=0; j<num_of_intervals; j++){
	      if (j < 1) {
		if (-i > intervals_4_approx[j]) {
		  constitutive->pol_search_help[i] = j ;
		  break ;
		}
	      }
	      else {
		if (-i > intervals_4_approx[j] && -i <= intervals_4_approx[j-1]) {
		  constitutive->pol_search_help[i] = j ;
		  break ;
		 }
	       }
	     }
	   }
	  
	  constitutive->polynomials_allocated = true;
	}
	//create matrix
	Aside = new double*[6 + pts];
	for (int i=0; i< (6 + pts); i++){
	  Aside[i] = new double[6+pts];
	}
	
	Bside = new double[6 + pts];


        for (int i=0; i<num_of_intervals; i++) {

	  if (i < 1){
	    for (int j=0; j<3; j++){
	      for (int k=0; k<(6 + pts); k++){
		Aside[j][k] = 0.0;
	      }
	    }
	   // Evaluate the second three rows of the matrix (zero, first and second derivative at point zero).
	    Aside[0][0] = 1.0 ;

	    // For the both first and second derivative it does not really matter what value is placed there.
	    Aside[1][1] = 1.0 ;
	    Aside[2][2] = 2.0 ;
	  
	  }
	  else {
	   for (int j=0; j<(6+pts); j++){
	     Aside[0][j] = pow(intervals_4_approx[i-1],j);
	     Aside[1][j] = j*pow(intervals_4_approx[i-1], j-1);
	     Aside[2][j] = j*(j-1)*pow(intervals_4_approx[i-1], j-2);
	   }
	  }
	  for (int j=0; j<(6+pts); j++){	   
	    Aside[3][j] = pow(intervals_4_approx[i],j);
	    Aside[4][j] = j*pow(intervals_4_approx[i], j-1);
	    Aside[5][j] =  j*(j-1)*pow(intervals_4_approx[i], j-2);
	    switch (pts) {
	      case 0 :
		break;
	      case 1 : 
		if (i > 0) {
		  Aside[6][j] =  pow((intervals_4_approx[i]+intervals_4_approx[i-1])/2, j);
		}
		else {
		  Aside[6][j] =  pow((intervals_4_approx[i])/2, j) ;
		}
		break;
	      default :
		printf("too many of inside points in polynomial approximation; not implemented!!! (check extras.cpp) \n");
		exit(1);
	    }
	  }
	  //Evaluate K(h) function.
	  if (i<1){
	    Bside[0] = constitutive->K(0.0, layer);
	    Bside[1] = 0.0;
	    Bside[2] = 0.0;
	  }
	  else {
	    Bside[0] = constitutive->K(intervals_4_approx[i-1], layer);
	    Bside[1] = constitutive->dKdh(intervals_4_approx[i-1], layer);
	    Bside[2] = constitutive->ddKdhh(intervals_4_approx[i-1], layer);
	  }
	 
	  Bside[3] = constitutive->K(intervals_4_approx[i], layer);
	  Bside[4] = constitutive->dKdh(intervals_4_approx[i], layer);
	  Bside[5] = constitutive->ddKdhh(intervals_4_approx[i], layer);  
	  
	  switch (pts) {
	      case 0 :
		break;
	      case 1 :
		if (i > 0) {
		  Bside[6] = constitutive->K((intervals_4_approx[i]+intervals_4_approx[i-1])/2, layer);
		}
		else {
		  Bside[6] = constitutive->K((intervals_4_approx[i])/2, layer);
		}
		break;
	  }



	  gem_full(Aside, Bside, constitutive->k_pols[i][layer][0], (6+pts));

	  for (int j=1; j<3; j++){
	    for (int k=0; k<5; k++){
	      constitutive->k_pols[i][layer][j][k] = (k+1)*constitutive->k_pols[i][layer][j-1][k+1];
	    }
	    constitutive->k_pols[i][layer][j][6-j] = 0.0;
	  }
	  
	  //Evaluate C(h) functions.
	  if (i<1){
	    Bside[0] = constitutive->C(0.0, layer);
	    Bside[1] = constitutive->dCdh(0.0, layer);
	  }
	  else {
	    Bside[0] = constitutive->C(intervals_4_approx[i-1], layer);
	    Bside[1] = constitutive->dCdh(intervals_4_approx[i-1], layer);
	  }
	 
	  //The first two lines of the matrix Aside stays the same.
	  for (int j=0; j<4; j++){	   
	    Aside[2][j] = pow(intervals_4_approx[i],j);
	    Aside[3][j] = j*pow(intervals_4_approx[i], j-1);
	  }

	  
	  Bside[2] = constitutive->C(intervals_4_approx[i], layer);
	  Bside[3] = constitutive->dCdh(intervals_4_approx[i], layer);
	 
	  
	  gem_full(Aside, Bside, constitutive->c_pols[i][layer][0], 4);
	  
	  for (int k=0; k<5; k++){
	    constitutive->c_pols[i][layer][1][k] = (k+1)*constitutive->c_pols[i][layer][0][k+1];
	  }
	  
	  constitutive->c_pols[i][layer][1][5] = 0.0;
	
    }
    delete [] Aside;
    delete [] Bside;
    break ;
  }
      
  return true;
}







 





