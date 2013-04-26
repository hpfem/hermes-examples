#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;

//  This is the fourth in the series of NIST benchmarks with known exact solutions.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: -Laplace u + f = 0.
//
//  Known exact solution: exp(-alpha * (pow(x - x_loc, 2) + pow(y - y_loc, 2))).
//  See functions CustomExactSolution::value and CustomExactSolution::derivatives in "exact_solution.cpp".
//
//  Domain: unit square (0, 1)x(0, 1), see the file "square_tri" or "square_quad.mesh".
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

// This problem has and exponential peak in the interior of the domain.
// (x_loc, y_loc) is the location of the peak, and alpha determines the strenghth of the peak. 
double alpha = 10;        
double x_loc = 0.5;         
double y_loc = 0.5;

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;
// This is a quantitative parameter of Adaptivity.
const double THRESHOLD = 0.75;
// This is a stopping criterion for Adaptivity.
const AdaptivityStoppingCriterion stoppingCriterion = AdaptStoppingCriterionSingleElement;

// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Maximum allowed level of hanging nodes.
const int MESH_REGULARITY = -1;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-2;
const CalculatedErrorType errorType = RelativeErrorToGlobalNorm;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  // Quadrilaterals.
  mloader.load("square_quad.mesh", mesh);
  // Triangles.
  // mloader.load("square_tri.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i<INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  basemesh->copy(mesh);

  // Set exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh, alpha, x_loc, y_loc));

  // Define right-hand side.
  CustomRightHandSide f(alpha, x_loc, y_loc);

  // Initialize weak formulation.
  WeakFormsH1::DefaultWeakFormPoissonLinear<double> wf(HERMES_ANY, &f);

  // Initialize boundary conditions
  DefaultEssentialBCNonConst<double> bc_essential("Bdy", exact_sln);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Initialize approximate solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>());
  MeshFunctionSharedPtr<double> ref_sln(new Solution<double>());
  
  // Initialize refinement selector.
  MySelector selector(CAND_LIST);

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 440, 350));
  sview.set_3d_mode(true);
  sview.fix_scale_width();
  Views::OrderView oview("Polynomial orders", new Views::WinGeom(450, 0, 420, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;

  // Assemble the discrete problem.    
  LinearSolver<double> newton;
  newton.set_weak_formulation(&wf);

  // Adaptivity loop.
  DefaultErrorCalculator<double, HERMES_H1_NORM> error_calculator(errorType, 1);
  Adapt<double> adaptivity(space, &error_calculator);
  adaptivity.set_strategy(stoppingCriterion, THRESHOLD);

  int period = 5;
  int number_of_steps = 100;

  for(int iteration = 0; iteration < number_of_steps; iteration++)
  {
    mesh->copy(basemesh);
    space->set_uniform_order(P_INIT);
    space->assign_dofs();

    double position = period - std::abs((iteration % (2*period)) - period);
    double factor = std::abs(std::sin( 0.5 * M_PI * std::pow(position / period, 4.0)));

    alpha = 10. + factor * 1000.;
    Hermes::Mixins::Loggable::Static::info("Iteration: %i, Position: %g, Factor: %g, Alpha: %g%.", iteration, position, factor, alpha);

    f.alpha = alpha;
    ((CustomExactSolution*)exact_sln.get())->alpha = alpha;

    int as = 1;
    while (!adaptive_step_single_space(mesh, space, sln, selector, ref_sln, cpu_time,newton,sview,oview,graph_dof_est,graph_cpu_est, error_calculator, adaptivity, as, ERR_STOP, exact_sln, graph_dof_exact, graph_cpu_exact))
    {
      Linearizer lin;
      char* filename = new char[1000];
      sprintf(filename, "Solution-%i.vtk", iteration);
      lin.save_solution_vtk(ref_sln, filename, "sln", true);
      Orderizer ord;
      sprintf(filename, "Orders-%i.vtk", iteration);
      ord.save_orders_vtk(space, filename);
      sprintf(filename, "Mesh-%i.vtk", iteration);
      ord.save_mesh_vtk(space, filename);
    }
    
    std::cout << std::endl << "------------------------------------------------" << std::endl;
  }
return 0;
}
