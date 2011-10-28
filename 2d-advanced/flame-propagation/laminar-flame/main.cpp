#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This example is a very simple flame propagation model (laminar flame,
//  zero flow velocity), and its purpose is to show how the Newton's method
//  is applied to a time-dependent two-equation system.
//
//  PDEs:
//
//  dT/dt - laplace T = omega(T,Y),
//  dY/dt - 1/Le * laplace Y = - omega(T,Y).
//
//  Domain: rectangle with cooled rods.
//
//  BC:  T = 1, Y = 0 on the inlet,
//       dT/dn = - kappa T on cooled rods,
//       dT/dn = 0, dY/dn = 0 elsewhere.
//
//  Time-stepping: a second order BDF formula.

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;                       
// Initial polynomial degree.
const int P_INIT = 2;                             
// Time step.
const double time_step = 0.01;                        
// Time interval length.
const double T_FINAL = 60.0;                      

// Newton's method.
const double DAMPING_COEFF = 1.0;
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-4;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 50;                   
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

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
ButcherTableType butcher_table_type = Implicit_RK_1;

// Problem constants.
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

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
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> left_t("Left", 1.0);
  EssentialBCs<double> bcs_t(&left_t);

  DefaultEssentialBCConst<double> left_c("Left", 0.0);
  EssentialBCs<double> bcs_c(&left_c);

  // Create H1 spaces with default shapesets.
  H1Space<double>* t_space = new H1Space<double>(&mesh, &bcs_t, P_INIT);
  H1Space<double>* c_space = new H1Space<double>(&mesh, &bcs_c, P_INIT);
  Hermes::vector<Space<double> *> spaces = Hermes::vector<Space<double> *>(t_space, c_space);
  int ndof = Space<double>::get_num_dofs(spaces);
  info("ndof = %d.", ndof);

  // Define initial conditions.
  InitialSolutionTemperature t_prev_time(&mesh, x1);
  InitialSolutionConcentration c_prev_time(&mesh, x1, Le);
  Hermes::vector<MeshFunction<double>*> meshfns_prev_time = Hermes::vector<MeshFunction<double>*>(&t_prev_time, &c_prev_time);
  Hermes::vector<Solution<double>*> slns_prev_time = Hermes::vector<Solution<double>*>(&t_prev_time, &c_prev_time);
  Solution<double> t_new_time(&mesh);
  Solution<double> c_new_time(&mesh);
  Hermes::vector<Solution<double>*> slns_new_time = Hermes::vector<Solution<double>*>(&t_new_time, &c_new_time);

  // Filters for the reaction rate omega and its derivatives.
  CustomFilter omega(slns_prev_time, Le, alpha, beta, kappa, x1);
  CustomFilterDt omega_dt(slns_prev_time, Le, alpha, beta, kappa, x1);
  CustomFilterDc omega_dc(slns_prev_time, Le, alpha, beta, kappa, x1);

  // Initialize visualization.
  ScalarView rview("Reaction rate", new WinGeom(0, 0, 800, 230));

  // Initialize weak formulation.
  CustomWeakForm wf(Le, alpha, beta, kappa, x1, &omega, &omega_dt, &omega_dc);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, spaces);

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(&dp, &bt, matrix_solver);

  // Time stepping:
  double current_time = 0.0;
  int ts = 1;
  bool jacobian_changed = true;
  do 
  {
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
      runge_kutta.rk_time_step_newton(current_time, time_step, slns_prev_time, 
                                  slns_new_time, freeze_jacobian, block_diagonal_jacobian, verbose,
                                  NEWTON_TOL, NEWTON_MAX_ITER, damping_coeff,
                                  max_allowed_residual_norm);
    }
    catch(Exceptions::Exception& e)
    {
      e.printMsg();
      error("Runge-Kutta time step failed");
    }

    // Saving solutions for the next time step.
    t_prev_time.copy(&t_new_time);
    c_prev_time.copy(&c_new_time);

    // Reinit filters.
    omega.reinit();
    omega_dc.reinit();
    omega_dt.reinit();

    // Visualization.
    rview.set_min_max_range(0.0,2.0);
    rview.show(&omega);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
