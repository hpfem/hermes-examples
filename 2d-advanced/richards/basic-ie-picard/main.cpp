#include "definitions.h"

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

// Picard's method.
// Number of last iterations used.
const int PICARD_NUM_LAST_ITER_USED = 3;
// Parameter for the Anderson acceleration.
const double PICARD_ANDERSON_BETA = 1.0;
// Stopping criterion for the Picard's method.
const double PICARD_TOL = 1e-6;
// Maximum allowed number of Picard iterations.
const int PICARD_MAX_ITER = 100;

int main(int argc, char* argv[])
{
  // Load the mesh->
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Initial mesh refinements.
  for (int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh->refine_all_elements();
  mesh->refine_towards_boundary("Top", INIT_REF_NUM_BDY);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential(std::vector<std::string>({"Bottom", "Right", "Top", "Left"}));
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();
  Hermes::Mixins::Loggable::Static::info("ndof = %d.", ndof);

  // Zero initial solutions. This is why we use H_OFFSET.
  MeshFunctionSharedPtr<double> h_time_prev(new ConstantSolution<double>(mesh, 0.));

  // Initialize views.
  ScalarView view("Initial condition", new WinGeom(0, 0, 600, 500));
  view.fix_scale_width(80);

  // Initialize the constitutive relations.
  ConstitutiveRelations* constitutive_relations;
  if (constitutive_relations_type == CONSTITUTIVE_GENUCHTEN)
    constitutive_relations = new ConstitutiveRelationsGenuchten(ALPHA, M, N, THETA_S, THETA_R, K_S, STORATIVITY);
  else
    constitutive_relations = new ConstitutiveRelationsGardner(ALPHA, THETA_S, THETA_R, K_S);

  // Initialize the weak formulation.
  double current_time = 0;
  WeakFormSharedPtr<double> wf(new CustomWeakFormRichardsIEPicard(time_step, h_time_prev, constitutive_relations));

  // Initialize the Picard solver.
  PicardSolver<double> picard(wf, space);
  picard.set_verbose_output(true);

  // Time stepping:
  int ts = 1;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f s", ts, current_time);

    // Perform the Picard's iteration (Anderson acceleration on by default).
    picard.set_max_allowed_iterations(PICARD_MAX_ITER);
    picard.set_num_last_vector_used(PICARD_NUM_LAST_ITER_USED);
    picard.set_anderson_beta(PICARD_ANDERSON_BETA);

    try
    {
      picard.solve();
    }
    catch (std::exception& e)
    {
      std::cout << e.what();
    }

    // Translate the coefficient vector into a Solution.
    Solution<double>::vector_to_solution(picard.get_sln_vector(), space, h_time_prev);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;

    // Visualize the solution.
    char title[100];
    sprintf(title, "Time %g s", current_time);
    view.set_title(title);
    view.show(h_time_prev);
  } while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}