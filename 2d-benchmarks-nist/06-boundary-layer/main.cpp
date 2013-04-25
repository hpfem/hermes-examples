#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;

//  This is the sixth in the series of NIST benchmarks with known exact solutions. It solves
//  a problem with boundary layer.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  The problem is made harder for adaptive algorithms by decreasing the (positive) parameter EPSILON.
//
//  PDE: -EPSILON Laplace u + 2du/dx + du/dy - f = 0
//
//  Known exact solution, see the class CustomExactSolution.
//
//  Domain: square (-1, 1) x (-1, 1), see the file square.mesh->
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

double epsilon = 1e1;

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// This is a quantitative parameter of Adaptivity.
const double THRESHOLD = 0.8;
// This is a stopping criterion for Adaptivity.
const AdaptivityStoppingCriterion stoppingCriterion = AdaptStoppingCriterionCumulative;

// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Maximum allowed level of hanging nodes.
const int MESH_REGULARITY = -1;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-3;
const CalculatedErrorType errorType = RelativeErrorToGlobalNorm;

// Newton tolerance
const double NEWTON_TOLERANCE = 1e-6;

bool HERMES_VISUALIZATION = false;
bool VTK_VISUALIZATION = false;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square_quad.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();
  
  basemesh->copy(mesh);

  // Set exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh, epsilon));

  // Define right-hand side.
  CustomRightHandSide f(epsilon);

  // Initialize weak formulation.
  CustomWeakForm wf(&f);

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

  // Assemble the discrete problem.    
  NewtonSolver<double> newton;
  newton.set_weak_formulation(&wf);

  // Adaptivity loop.
  DefaultErrorCalculator<double, HERMES_H1_NORM> error_calculator(errorType, 1);
  Adapt<double> adaptivity(space, &error_calculator);
  adaptivity.set_strategy(stoppingCriterion, THRESHOLD);

  int number_of_steps = 14;
  for(int iteration = 0; iteration < number_of_steps; iteration++)
  {
    mesh->copy(basemesh);
    space->set_uniform_order(P_INIT);
    space->assign_dofs();

    if(iteration < number_of_steps / 2)
      epsilon = epsilon / 2.;
    else if(iteration >= number_of_steps / 2 && iteration < number_of_steps)
      epsilon = epsilon * 2.;
    else if(iteration >= number_of_steps && iteration < number_of_steps * 3/2)
      epsilon = epsilon / 2.;
    else
      epsilon = epsilon * 2.;

    Hermes::Mixins::Loggable::Static::info("Iteration: %i, Epsilon: %g.", iteration, epsilon);
    f.epsilon = epsilon;
    ((CustomExactSolution*)exact_sln.get())->epsilon = epsilon;

    int as = 1;
    while (!adaptive_step_single_space(mesh, space, sln, selector, ref_sln, cpu_time,newton,sview,oview,graph_dof_est,graph_cpu_est, error_calculator, adaptivity,as, ERR_STOP))
    {
      Linearizer lin;
      char* filename = new char[1000];
      sprintf(filename, "Solution-%i.vtk", iteration);
      lin.save_solution_vtk(ref_sln, filename, "sln", true);
      Orderizer ord;
      sprintf(filename, "Orders-%i.vtk", iteration);
      ord.save_orders_vtk(newton.get_space(0), filename);
      sprintf(filename, "Mesh-%i.vtk", iteration);
      ord.save_mesh_vtk(newton.get_space(0), filename);
    }
    
    std::cout << std::endl << "------------------------------------------------" << std::endl;
  }
  
  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
