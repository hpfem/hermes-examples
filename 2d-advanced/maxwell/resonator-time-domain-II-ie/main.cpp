#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// This example solves a time-domain resonator problem for the Maxwell's equation. 
// It is very similar to resonator-time-domain-I but B is eliminated from the 
// equations, thus converting the first-order system into one second -order
// equation in time. The second-order equation in time is decomposed back into 
// a first-order system in time in the standard way (see example wave-1). Time 
// discretization is performed using the implicit Euler method.
//
// PDE: \frac{1}{SPEED_OF_LIGHT**2}\frac{\partial^2 E}{\partial t^2} + curl curl E = 0,
// converted into
//
//      \frac{\partial E}{\partial t} - F = 0,
//      \frac{\partial F}{\partial t} + SPEED_OF_LIGHT**2 * curl curl E = 0.
//
// Approximated by
// 
//      \frac{E^{n+1} - E^{n}}{tau} - F^{n+1} = 0,
//      \frac{F^{n+1} - F^{n}}{tau} + SPEED_OF_LIGHT**2 * curl curl E^{n+1} = 0.
//
// Domain: Square (-pi/2, pi/2) x (-pi/2, pi/2)... See mesh file domain.mesh.
//
// BC:  E \times \nu = 0 on the boundary (perfect conductor),
//      F \times \nu = 0 on the boundary (E \times \nu = 0 => \partial E / \partial t \times \nu = 0).
//
// IC:  Prescribed wave for E, zero for F.
//
// The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 6;                              
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;                        
// Time step.
const double time_step = 0.05;                     
// Final time.
const double T_FINAL = 35.0;                       
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-8;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                  
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   

// Problem parameters.
// Square of wave speed.  
const double C_SQUARED = 1;                                         

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize solutions.
  CustomInitialConditionWave E_sln(&mesh);
  ZeroSolutionVector F_sln(&mesh);
  Hermes::vector<Solution<double>*> slns(&E_sln, &F_sln);

  // Initialize the weak formulation.
  CustomWeakFormWaveIE wf(time_step, C_SQUARED, &E_sln, &F_sln);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential("Perfect conductor", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create x- and y- displacement space using the default H1 shapeset.
  HcurlSpace<double> E_space(&mesh, &bcs, P_INIT);
  HcurlSpace<double> F_space(&mesh, &bcs, P_INIT);
  Hermes::vector<const Space<double> *> spaces = Hermes::vector<const Space<double> *>(&E_space, &F_space);
  int ndof = HcurlSpace<double>::get_num_dofs(spaces);
  info("ndof = %d.", ndof);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, spaces);

  // Project the initial condition on the FE space to obtain initial 
  // coefficient vector for the Newton's method.
  // NOTE: If you want to start from the zero vector, just define 
  // coeff_vec to be a vector of ndof zeros (no projection is needed).
  info("Projecting to obtain initial vector for the Newton's method.");
  double* coeff_vec = new double[ndof];
  OGProjection<double>::project_global(spaces, slns, coeff_vec, matrix_solver); 

  // Initialize Newton solver.
  NewtonSolver<double> newton(&dp, matrix_solver);

  // Initialize views.
  ScalarView E1_view("Solution E1", new WinGeom(0, 0, 400, 350));
  E1_view.fix_scale_width(50);
  ScalarView E2_view("Solution E2", new WinGeom(410, 0, 400, 350));
  E2_view.fix_scale_width(50);
  ScalarView F1_view("Solution F1", new WinGeom(0, 410, 400, 350));
  F1_view.fix_scale_width(50);
  ScalarView F2_view("Solution F2", new WinGeom(410, 410, 400, 350));
  F2_view.fix_scale_width(50);

  // Time stepping loop.
  double current_time = 0; int ts = 1;
  do
  {
    // Perform one implicit Euler time step.
    info("Implicit Euler time step (t = %g s, time_step = %g s).", current_time, time_step);

    // Perform Newton's iteration.
    try
    {
      newton.solve_keep_jacobian(coeff_vec, NEWTON_TOL, NEWTON_MAX_ITER);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    }

    // Translate the resulting coefficient vector into Solutions.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces, slns);

    // Visualize the solutions.
    char title[100];
    sprintf(title, "E1, t = %g", current_time + time_step);
    E1_view.set_title(title);
    E1_view.show(&E_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_0);
    sprintf(title, "E2, t = %g", current_time + time_step);
    E2_view.set_title(title);
    E2_view.show(&E_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_1);

    sprintf(title, "F1, t = %g", current_time + time_step);
    F1_view.set_title(title);
    F1_view.show(&F_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_0);
    sprintf(title, "F2, t = %g", current_time + time_step);
    F2_view.set_title(title);
    F2_view.show(&F_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_1);

    //View::wait();

    // Update time.
    current_time += time_step;

  } while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();

  return 0;
}
