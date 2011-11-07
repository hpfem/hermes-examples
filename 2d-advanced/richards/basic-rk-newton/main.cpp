#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"



//  This example solves the Tracy problem with arbitrary Runge-Kutta 
//  methods in time. 
//
//  PDE: C(h)dh/dt - div(K(h)grad(h)) - (dK/dh)*(dh/dy) = 0
//  where K(h) = K_S*exp(alpha*h)                          for h < 0,
//        K(h) = K_S                                       for h >= 0,
//        C(h) = alpha*(theta_s - theta_r)*exp(alpha*h)    for h < 0,
//        C(h) = alpha*(theta_s - theta_r)                 for h >= 0.
//
//  Domain: square (0, 100)^2.
//
//  BC: Dirichlet, given by the initial condition.
//  IC: Flat in all elements except the top layer, within this 
//      layer the solution rises linearly to match the Dirichlet condition.
//
//  NOTE: The pressure head 'h' is between -1000 and 0. For convenience, we
//        increase it by an offset H_OFFSET = 1000. In this way we can start
//        from a zero coefficient vector.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_GLOB_REF_NUM = 3;                  
// Number of initial refinements towards boundary.
const int INIT_REF_NUM_BDY = 5;                   
// Initial polynomial degree.
const int P_INIT = 2;                             
// Time step.
double time_step = 5e-4;                          
// Time interval length.
const double T_FINAL = 0.4;                       
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Problem parameters.
double K_S = 20.464;
double ALPHA = 0.001;
double THETA_R = 0;
double THETA_S = 0.45;

double M, N, STORATIVITY;

// Constitutive relations.
enum CONSTITUTIVE_RELATIONS {
    CONSTITUTIVE_GENUCHTEN,    // Van Genuchten.
    CONSTITUTIVE_GARDNER       // Gardner.
};
// Use van Genuchten's constitutive relations, or Gardner's.
CONSTITUTIVE_RELATIONS constitutive_relations_type = CONSTITUTIVE_GARDNER;

// Newton's method.
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-5;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 500;                  
const double DAMPING_COEFF = 1.0;

// Choose one of the following time-integration methods, or define your own Butcher's table. The last number 
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods: 
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2, 
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4, 
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_5_4.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
//   Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded, 
//   Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded, 
//   Implicit_DIRK_ISMAIL_7_45_embedded. 
ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Top", INIT_REF_NUM_BDY);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential(Hermes::vector<std::string>("Bottom", "Right", "Top", "Left"));
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d.", ndof);

  // Zero initial solutions. This is why we use H_OFFSET.
  ZeroSolution h_time_prev(&mesh), h_time_new(&mesh);

  // Initialize views.
  ScalarView view("Initial condition", new WinGeom(0, 0, 600, 500));
  view.fix_scale_width(80);

  // Visualize the initial condition.
  view.show(&h_time_prev);

  // Initialize the constitutive relations.
  ConstitutiveRelations* constitutive_relations;
  if(constitutive_relations_type == CONSTITUTIVE_GENUCHTEN)
    constitutive_relations = new ConstitutiveRelationsGenuchten(ALPHA, M, N, THETA_S, THETA_R, K_S, STORATIVITY);
  else
    constitutive_relations = new ConstitutiveRelationsGardner(ALPHA, THETA_S, THETA_R, K_S);

  // Initialize the weak formulation.
  CustomWeakFormRichardsRK wf(constitutive_relations);

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(&wf, &space, &bt, matrix_solver);

  // Time stepping:
  double current_time = 0;
  int ts = 1;
  do 
  {
    info("---- Time step %d, time %3.5f s", ts, current_time);

    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, time step = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool freeze_jacobian = false;
    bool block_diagonal_jacobian = false;
    bool verbose = true;
    double damping_coeff = 1.0;
    double max_allowed_residual_norm = 1e10;

    try
    {
      runge_kutta.rk_time_step_newton(current_time, time_step, &h_time_prev, 
          &h_time_new, freeze_jacobian, block_diagonal_jacobian, verbose,
          NEWTON_TOL, NEWTON_MAX_ITER, damping_coeff, max_allowed_residual_norm);
    }
    catch(Exceptions::Exception& e)
    {
      e.printMsg();
      error("Runge-Kutta time step failed");
    }

    // Copy solution for the new time step.
    h_time_prev.copy(&h_time_new);

    // Increase current time.
    current_time += time_step;

    // Visualize the solution.
    char title[100];
    sprintf(title, "Time %3.2f s", current_time);
    view.set_title(title);
    view.show(&h_time_prev);

    // Increase time step counter.
    ts++;
  }
  while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

