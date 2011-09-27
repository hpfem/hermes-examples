#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Global variables for forms.
double K_S, ALPHA, THETA_R, THETA_S, N, M;

// Problem parameters.
const double TAU = 5e-3;                          // Time step.
const double STARTUP_TIME = 1.1e-2;               // Start-up time for time-dependent Dirichlet boundary condition.
const double T_FINAL = 5.0;                       // Time interval length.
double TIME = 0;                                  // Global time variable initialized with first time step.
double H_INIT = -9.5;                             // Initial pressure head.
double H_ELEVATION = 5.2;

double K_S_1 = 0.108;
double K_S_3 = 0.0048;
double K_S_2 = 0.0168;
double K_S_4 = 1.061;

double ALPHA_1 = 0.01;
double ALPHA_3 = 0.005;
double ALPHA_2 = 0.01;
double ALPHA_4 = 0.05;

double THETA_R_1 = 0.1020;
double THETA_R_2 = 0.09849;
double THETA_R_3 = 0.08590;
double THETA_R_4 = 0.08590;

double THETA_S_1 = 0.4570;
double THETA_S_2 = 0.4510;
double THETA_S_3 = 0.4650;
double THETA_S_4 = 0.5650;

double N_1 = 1.982;
double N_2 = 1.632; 
double N_3 = 5.0;
double N_4 = 5.0;

double M_1 = 0.49546;
double M_2 = 0.38726;
double M_3 = 0.8;
double M_4 = 0.8;