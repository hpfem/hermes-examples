#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;

//  This example is similar to basic-ie-newton except it uses the 
//  Picard's method in each time step.
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

// If this is defined, use van Genuchten's constitutive relations, otherwise use Gardner's.
//#define CONSTITUTIVE_GENUCHTEN

const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 5;                   // Number of initial refinements towards boundary.
const int P_INIT = 2;                             // Initial polynomial degree.
double time_step = 5e-4;                          // Time step.
const double T_FINAL = 0.4;                       // Time interval length.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Picard's method.
const int PICARD_NUM_LAST_ITER_USED = 3;          // Number of last iterations used.
const double PICARD_ANDERSON_BETA = 1.0;          // Parameter for the Anderson acceleration. 
const double PICARD_TOL = 1e-6;                   // Stopping criterion for the Picard's method.
const int PICARD_MAX_ITER = 100;                  // Maximum allowed number of Picard iterations.

int main(int argc, char* argv[])
{
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

  // Initial condition vector is the zero vector. This is why we
  // use the H_OFFSET. 
  double* coeff_vec = new double[ndof];
  memset(coeff_vec, 0, ndof*sizeof(double));

  // Convert initial condition into a Solution.
  Solution<double> h_time_prev, h_iter_prev;
  Solution<double>::vector_to_solution(coeff_vec, &space, &h_time_prev);
  Solution<double>::vector_to_solution(coeff_vec, &space, &h_iter_prev);

  // Initialize views.
  ScalarView view("Initial condition", new WinGeom(0, 0, 600, 500));
  view.fix_scale_width(80);

  // Visualize the initial condition.
  view.show(&h_time_prev);

  // Initialize the weak formulation.
  double current_time = 0;
  CustomWeakFormRichardsIEPicard wf(time_step, &h_time_prev, &h_iter_prev);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Initialize the Picard solver.
  PicardSolver<double> picard(&dp, &h_iter_prev, matrix_solver);
  picard.set_verbose_output(true);

  // Time stepping:
  int ts = 1;
  do 
  {
    info("---- Time step %d, time %3.5f s", ts, current_time);

    // Perform the Picard's iteration (Anderson acceleration on by default).
    if (!picard.solve(PICARD_TOL, PICARD_MAX_ITER, PICARD_NUM_LAST_ITER_USED, 
                      PICARD_ANDERSON_BETA)) error("Picard's iteration failed.");

    // Translate the coefficient vector into a Solution. 
    Solution<double>::vector_to_solution(picard.get_sln_vector(), &space, &h_iter_prev);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;

    // Visualize the solution.
    char title[100];
    sprintf(title, "Time %g s", current_time);
    view.set_title(title);
    view.show(&h_iter_prev);

    // Save the next time level solution.
    h_time_prev.copy(&h_iter_prev);
  }
  while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

