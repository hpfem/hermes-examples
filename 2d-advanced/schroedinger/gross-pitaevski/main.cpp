#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This example uses the Newton's method to solve a nonlinear complex-valued
//  time-dependent PDE (the Gross-Pitaevski equation describing the behavior
//  of Einstein-Bose quantum gases). For time-discretization one can use either
//  the first-order implicit Euler method or the second-order Crank-Nicolson
//  method.
//
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates.
//
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi.
//
//  Domain: square (-1, 1)^2.
//
//  BC:  homogeneous Dirichlet everywhere on the boundary.

// Number of initial uniform refinements.
const int INIT_REF_NUM = 3;                       
// Initial polynomial degree.
const int P_INIT = 4;                             
// Time step.
double time_step = 0.005;                         
// Time interval length.
const double T_FINAL = 2;                         
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-5;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                  
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
ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

// Problem constants
// Planck constant 6.626068e-34.
const double h = 1;                               
// Mass of boson.
const double m = 1;                               
// Coupling constant.
const double g = 1;                               
// Frequency.
const double omega = 1;                           

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
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Convert initial condition into a Solution<std::complex<double> >.
  CustomInitialCondition psi_time_prev(&mesh);
  Solution<std::complex<double> > psi_time_new(&mesh);

  // Initialize the weak formulation.
  double current_time = 0;

  CustomWeakFormGPRK wf(h, m, g, omega);
  
  // Initialize boundary conditions.
  DefaultEssentialBCConst<std::complex<double> > bc_essential("Bdy", 0.0);
  EssentialBCs<std::complex<double> > bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<std::complex<double> > space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);
 
  // Initialize the FE problem.
  DiscreteProblem<std::complex<double> > dp(&wf, &space);

  // Initialize views.
  ScalarView sview_real("Solution - real part", new WinGeom(0, 0, 600, 500));
  ScalarView sview_imag("Solution - imaginary part", new WinGeom(610, 0, 600, 500));
  sview_real.fix_scale_width(80);
  sview_imag.fix_scale_width(80);

  // Initialize Runge-Kutta time stepping.
  RungeKutta<std::complex<double> > runge_kutta(&wf, &space, &bt, matrix_solver);
  
  // Time stepping:
  int ts = 1;
  int nstep = (int)(T_FINAL/time_step + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, time step = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool jacobian_changed = false;
    bool verbose = true;
    
    try
    {
      runge_kutta.rk_time_step_newton(current_time, time_step, &psi_time_prev, 
                                  &psi_time_new, !jacobian_changed, false, verbose);
    }
    catch(Exceptions::Exception& e)
    {
      e.printMsg();
      error("Runge-Kutta time step failed");
    }

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Solution - real part, Time %3.2f s", current_time);
    sview_real.set_title(title);
    sprintf(title, "Solution - imaginary part, Time %3.2f s", current_time);
    sview_imag.set_title(title);
    RealFilter real(&psi_time_new);
    ImagFilter imag(&psi_time_new);
    sview_real.show(&real);
    sview_imag.show(&imag);

    // Copy solution for the new time step.
    psi_time_prev.copy(&psi_time_new);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;
  }

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
