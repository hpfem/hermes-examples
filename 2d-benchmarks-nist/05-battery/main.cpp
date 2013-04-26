#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;

//  This is the fifth in the series of NIST benchmarks with unknown exact solution.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: -\frac{\partial }{\partial x}\left(p(x, y)\frac{\partial u}{\partial x}\right)
//       -\frac{\partial }{\partial y}\left(q(x, y)\frac{\partial u}{\partial y}\right) - f = 0.
//
//  Exact solution: unknown.
//
//  Domain: square (0, 8.4) x (0, 24), see the file "battery.mesh".
//
//  BC: Zero Neumann on left edge, Newton on the rest of the boundary:
//
//  The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// This is a quantitative parameter of Adaptivity.
const double THRESHOLD = 0.3;
// This is a stopping criterion for Adaptivity.
const AdaptivityStoppingCriterion stoppingCriterion = AdaptStoppingCriterionSingleElement;   

// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO_H;
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
  mloader.load("battery.mesh", mesh);
  MeshView m;
  m.show(mesh);
  View::wait();

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  basemesh->copy(mesh);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, P_INIT));

  // Initialize weak formulation.
  CustomWeakFormPoisson wf("e1", "e2", "e3", "e4", "e5", 
                           "Bdy_left", "Bdy_top", "Bdy_right", "Bdy_bottom", mesh);

  // Initialize coarse and fine mesh solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>), ref_sln(new Solution<double>);
  
  // Initialize refinement selector.
  MySelector selector(CAND_LIST);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 320, 600));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  OrderView  oview("Polynomial orders", new WinGeom(330, 0, 300, 600));
  
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

  int period = 5;
  int number_of_steps = 100;
  for(int iteration = 0; iteration < number_of_steps; iteration++)
  {
    mesh->copy(basemesh);
    space->set_uniform_order(P_INIT);
    space->assign_dofs();

    double position = period - std::abs((iteration % (2*period)) - period);
    double factor = std::abs(std::sin( 0.5 * M_PI * std::pow(position / period, 4.0)));

    /*
    alpha = 10. + factor * 1000.;
    Hermes::Mixins::Loggable::Static::info("Iteration: %i, Position: %g, Factor: %g, Alpha: %g%.", iteration, position, factor, alpha);

    f.alpha = alpha;
    ((CustomExactSolution*)exact_sln.get())->alpha = alpha;
    */

    int as = 1;
    while (!adaptive_step_single_space(mesh, space, sln, selector, ref_sln, cpu_time,newton,sview,oview,graph_dof_est,graph_cpu_est, error_calculator, adaptivity,as, ERR_STOP))
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

  
  Views::View::wait();

  return 0;
}

