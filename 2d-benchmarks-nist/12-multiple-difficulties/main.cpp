#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;

//  This is the twelfth in the series of NIST benchmarks with known exact solutions.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: -Laplace u + f = 0 
//
//  Known exact solution: 
//  See the class CustomExactSolution::value in file "definitions.h"
//
//  Domain: L-shaped domain (-1,1)x(-1,1)\(0,1)x(-1,0), see the file "lshape.mesh".
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

const double omega_c = 3.0 * M_PI / 2.0;
const double x_w = 0.0;             
const double y_w = -3.0 / 4.0;
const double r_0 = 3.0 / 4.0;
const double alpha_w = 200.0;
const double x_p = -Hermes::sqrt(5.0) / 4.0;
const double y_p = -1.0 / 4.0;
const double alpha_p = 1000.0;
const double epsilon = 1.0 / 100.0;

// Initial polynomial degree of mesh elements.
const int P_INIT = 2;  
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0; 
// This is a quantitative parameter of Adaptivity.
const double THRESHOLD = 0.8;
// This is a stopping criterion for Adaptivity.
const AdaptivityStoppingCriterion stoppingCriterion = AdaptStoppingCriterionLevels;   

// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Maximum allowed level of hanging nodes.
const int MESH_REGULARITY = -1;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;
const CalculatedErrorType errorType = RelativeErrorToGlobalNorm;

// Newton tolerance
const double NEWTON_TOLERANCE = 1e-6;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("lshape.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();
  basemesh->copy(mesh);

  // Set exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh, alpha_w, alpha_p, x_w, y_w, r_0, omega_c, epsilon, x_p, y_p));

  // Define right-hand side.
  CustomRightHandSide f(alpha_w, alpha_p, x_w, y_w, r_0, omega_c, epsilon, x_p, y_p);

  // Initialize the weak formulation.
  Hermes1DFunction<double> lambda(1.0);
  WeakFormsH1::DefaultWeakFormPoisson<double> wf(HERMES_ANY, &lambda, &f);

  // Initialize boundary conditions
  DefaultEssentialBCNonConst<double> bc_essential("Bdy", exact_sln);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Initialize approximate solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>), ref_sln(new Solution<double>);

  // Initialize refinement selector.
  MySelector selector(CAND_LIST);

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 440, 350));
  sview.show_mesh(false);
  sview.fix_scale_width(50);
  Views::OrderView  oview("Polynomial orders", new Views::WinGeom(450, 0, 420, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;

  // Adaptivity loop:
  NewtonSolver<double> newton;
  newton.set_weak_formulation(&wf);

  // Adaptivity loop.
  DefaultErrorCalculator<double, HERMES_H1_NORM> error_calculator(errorType, 1);
  Adapt<double> adaptivity(space, &error_calculator);
  adaptivity.set_strategy(stoppingCriterion);

  int number_of_steps = 20;
  for(int iteration = 0; iteration < number_of_steps; iteration++)
  {
    mesh->copy(basemesh);
    space->set_uniform_order(P_INIT);
    space->assign_dofs();

    /*
    Hermes::Mixins::Loggable::Static::info("Iteration: %i, Alpha: %g.", iteration, alpha);
    f.alpha = alpha;
    ((CustomExactSolution*)exact_sln.get())->alpha = alpha;
    */

    int as = 1;
    while (!adaptive_step_single_space(mesh, space, sln, selector, ref_sln, cpu_time,newton,sview,oview,graph_dof_est,graph_cpu_est, error_calculator, adaptivity,as, ERR_STOP))
    {
      Linearizer lin;
      char* filename = new char[1000];
      sprintf(filename, "Solution-%i.vtk", iteration);
      lin.save_solution_vtk(ref_sln, filename, "sln", false);
      Orderizer ord;
      sprintf(filename, "Orders-%i.vtk", iteration);
      ord.save_orders_vtk(newton.get_space(0), filename);
      sprintf(filename, "Mesh-%i.vtk", iteration);
      ord.save_mesh_vtk(newton.get_space(0), filename);
    }

    std::cout << std::endl << "------------------------------------------------" << std::endl;
  }
  return 0;
}
