#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// This example is a simple test case for the Debye-Maxwell model solved in terms of 
// E, H and P. Here E is electric field (vector), H magnetic field (scalar), and P 
// electric polarization (vector). The example comes with a known exact solution. 
// Time discretization is performed using an arbitrary Runge-Kutta method.
//
// PDE system:
//
//      \frac{partial H}{\partial t} + \frac{1}{\mu_0} \curl E = 0,
//      \frac{\partial E}{\partial t} - \frac{1}{\epsilon_0 \epsilon_\infty} \curl H
//          + \frac{\epsilon_q - 1}{\tau}E - \frac{1}{\tau} P = 0,
//      \frac{\partial P}{\partial t} - \frac{(\epsilon_q - 1)\epsilon_0 \epsilon_\infty}{\tau} E
//          + \frac{1}{\tau} P = 0.
//
// Domain: Square (0, 1) x (0, 1).
//
// BC:  Perfect conductor for E and P.
//
// IC:  Prescribed functions for E, H and P.
//
// The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 4;                              
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;                        
// Time step.
const double time_step = 0.00001;                     
// Final time.
const double T_FINAL = 35.0;                       
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
ButcherTableType butcher_table = Implicit_RK_1;
//ButcherTableType butcher_table = Implicit_SDIRK_2_2;

// Problem parameters.
// Permeability of free space.
const double MU_0 = 1.0;
// Permittivity of free space.
const double EPS_0 = 1.0;
// Permittivity at infinite frequency.
const double EPS_INF = 1.0;
// Permittivity at zero frequency.
const double EPS_S = 2.0;
// EPS_Q.
const double EPS_Q = EPS_S / EPS_INF;
// Relaxation time of the medium.
const double TAU = 1.0;
// Angular frequency. Depends on wave number K. Must satisfy: 
// omega^3 - 2 omega^2 + K^2 M_PI^2 omega - K^2 M_PI^2 = 0.
// WARNING; Choosing wrong omega may lead to K**2 < 0.
const double OMEGA = 1.5;
// Wave vector direction (will be normalized to be compatible
// with omega).
double K_x = 1.0;
double K_y = 1.0;

int main(int argc, char* argv[])
{
  // Sanity check for omega. 
  double K_squared = Hermes::sqr(OMEGA/M_PI) * (OMEGA - 2) / (1 - OMEGA);
  if (K_squared <= 0) throw Hermes::Exceptions::Exception("Wrong choice of omega, K_squared < 0!");
  double K_norm_coeff = std::sqrt(K_squared) / std::sqrt(Hermes::sqr(K_x) + Hermes::sqr(K_y));
  Hermes::Mixins::Loggable::Static::info("Wave number K = %g", std::sqrt(K_squared));
  K_x *= K_norm_coeff;
  K_y *= K_norm_coeff;

  // Wave number.
  double K = std::sqrt(Hermes::sqr(K_x) + Hermes::sqr(K_y));

  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table);
  if (bt.is_explicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize solutions.
  double current_time = 0;
  CustomInitialConditionE E_time_prev(&mesh, current_time, OMEGA, K_x, K_y);
  CustomInitialConditionH H_time_prev(&mesh, current_time, OMEGA, K_x, K_y);
  CustomInitialConditionP P_time_prev(&mesh, current_time, OMEGA, K_x, K_y);
  Hermes::vector<Solution<double>*> slns_time_prev(&E_time_prev, &H_time_prev, &P_time_prev);
  Solution<double> E_time_new(&mesh), H_time_new(&mesh), P_time_new(&mesh);
  Hermes::vector<Solution<double>*> slns_time_new(&E_time_new, &H_time_new, &P_time_new);

  // Initialize the weak formulation.
  const CustomWeakFormMD wf(OMEGA, K_x, K_y, MU_0, EPS_0, EPS_INF, EPS_Q, TAU);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create x- and y- displacement space using the default H1 shapeset.
  HcurlSpace<double> E_space(&mesh, &bcs, P_INIT);
  H1Space<double> H_space(&mesh, NULL, P_INIT);
  //L2Space<double> H_space(&mesh, P_INIT);
  HcurlSpace<double> P_space(&mesh, &bcs, P_INIT);

  Hermes::vector<Space<double> *> spaces = Hermes::vector<Space<double> *>(&E_space, &H_space, &P_space);
  Hermes::vector<const Space<double> *> spaces_mutable = Hermes::vector<const Space<double> *>(&E_space, &H_space, &P_space);
  int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof = %d.", ndof);
  
  // Initialize views.
  ScalarView E1_view("Solution E1", new WinGeom(0, 0, 400, 350));
  E1_view.fix_scale_width(50);
  ScalarView E2_view("Solution E2", new WinGeom(410, 0, 400, 350));
  E2_view.fix_scale_width(50);
  ScalarView H_view("Solution H", new WinGeom(0, 410, 400, 350));
  H_view.fix_scale_width(50);
  ScalarView P1_view("Solution P1", new WinGeom(410, 410, 400, 350));
  P1_view.fix_scale_width(50);
  ScalarView P2_view("Solution P2", new WinGeom(820, 410, 400, 350));
  P2_view.fix_scale_width(50);

  // Visualize initial conditions.
  char title[100];
  sprintf(title, "E1 - Initial Condition");
  E1_view.set_title(title);
  E1_view.show(&E_time_prev, H2D_FN_VAL_0);
  sprintf(title, "E2 - Initial Condition");
  E2_view.set_title(title);
  E2_view.show(&E_time_prev, H2D_FN_VAL_1);

  sprintf(title, "H - Initial Condition");
  H_view.set_title(title);
  H_view.show(&H_time_prev);

  sprintf(title, "P1 - Initial Condition");
  P1_view.set_title(title);
  P1_view.show(&P_time_prev, H2D_FN_VAL_0);
  sprintf(title, "P2 - Initial Condition");
  P2_view.set_title(title);
  P2_view.show(&P_time_prev, H2D_FN_VAL_1);

  View::wait(HERMES_WAIT_KEYPRESS);

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(&wf, spaces_mutable, &bt);

  // Time stepping loop.
  int ts = 1;
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    Hermes::Mixins::Loggable::Static::info("Runge-Kutta time step (t = %g s, time_step = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool freeze_jacobian = false;
    bool block_diagonal_jacobian = false;
    bool verbose = true;
    double damping_coeff = 1.0;
    double max_allowed_residual_norm = 1e10;

    try
    {
      runge_kutta.setTime(current_time);
      runge_kutta.setTimeStep(time_step);
      runge_kutta.set_newton_max_iter(NEWTON_MAX_ITER);
      runge_kutta.set_newton_tol(NEWTON_TOL);
      runge_kutta.set_verbose_output(true);
      runge_kutta.rk_time_step_newton(slns_time_prev, slns_time_new);
    }
    catch(Exceptions::Exception& e)
    {
      e.printMsg();
      throw Hermes::Exceptions::Exception("Runge-Kutta time step failed");
    }

    // Visualize the solutions.
    char title[100];
    sprintf(title, "E1, t = %g", current_time + time_step);
    E1_view.set_title(title);
    E1_view.show(&E_time_new, H2D_FN_VAL_0);
    sprintf(title, "E2, t = %g", current_time + time_step);
    E2_view.set_title(title);
    E2_view.show(&E_time_new, H2D_FN_VAL_1);

    sprintf(title, "H, t = %g", current_time + time_step);
    H_view.set_title(title);
    H_view.show(&H_time_new);

    sprintf(title, "P1, t = %g", current_time + time_step);
    P1_view.set_title(title);
    P1_view.show(&P_time_new, H2D_FN_VAL_0);
    sprintf(title, "P2, t = %g", current_time + time_step);
    P2_view.set_title(title);
    P2_view.show(&P_time_new, H2D_FN_VAL_1);

    //View::wait();

    // Update solutions.
    E_time_prev.copy(&E_time_new);
    H_time_prev.copy(&H_time_new);
    P_time_prev.copy(&P_time_new);

    // Update time.
    current_time += time_step;
  
  } while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();

  return 0;
}
