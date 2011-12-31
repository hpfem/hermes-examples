#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"

#include "definitions.h"

//  This example solves a simple version of the time-dependent
//  Richard's equation using the backward Euler method in time 
//  combined with the Newton's method in each time step. It describes
//  infiltration into an initially dry soil. The example has a exact 
//  solution that is given in terms of a Fourier series (see a paper 
//  by Tracy). The exact solution is not used here.
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
const double NEWTON_TOL = 1e-6;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                 
const double DAMPING_COEFF = 1.0;

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

  // Zero initial solutions. This is why we use H_OFFSET.
  ZeroSolution<double> h_time_prev(&mesh);

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
  double current_time = 0;
  CustomWeakFormRichardsIE wf(time_step, &h_time_prev, constitutive_relations);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Initialize Newton solver.
  NewtonSolver<double> newton(&dp, matrix_solver);
  newton.set_verbose_output(true);

  // Time stepping:
  int ts = 1;
  do 
  {
    info("---- Time step %d, time %3.5f s", ts, current_time);

    // Perform Newton's iteration.
    try
    {
      // NULL = zero initial coefficient vector.
      newton.solve(NULL, NEWTON_TOL, NEWTON_MAX_ITER);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    };

    // Translate the resulting coefficient vector into the Solution<double> sln.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &h_time_prev);

    // Visualize the solution.
    char title[100];
    sprintf(title, "Time %g s", current_time);
    view.set_title(title);
    view.show(&h_time_prev);

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

