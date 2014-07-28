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
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// This is a quantitative parameter of Adaptivity.
const double THRESHOLD = 0.3;
// This is a stopping criterion for Adaptivity.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);

// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO_H;
// Maximum allowed level of hanging nodes.
const int MESH_REGULARITY = -1;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1.0;
const CalculatedErrorType errorType = RelativeErrorToGlobalNorm;

// Newton tolerance
const double NEWTON_TOLERANCE = 1e-6;

bool HERMES_VISUALIZATION = false;
bool VTK_VISUALIZATION = false;

int main(int argc, char* argv[])
{
  // Load the mesh->
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("battery.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

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

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Time measurement.
    cpu_time.tick();

    // Construct globally refined mesh and setup fine mesh space->
    Mesh::ReferenceMeshCreator refMeshCreator(mesh);
    MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh);
    SpaceSharedPtr<double> ref_space = refSpaceCreator.create_ref_space();
    int ndof_ref = ref_space->get_num_dofs();

    // Initialize fine mesh problem.
    Hermes::Mixins::Loggable::Static::info("Solving on fine mesh->");
    DiscreteProblem<double> dp(wf, ref_space);

    NewtonSolver<double> newton(&dp);
    newton.set_verbose_output(true);

    // Perform Newton's iteration.
    try
    {
      newton.solve();
    }
    catch (Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    }

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh->
    Hermes::Mixins::Loggable::Static::info("Projecting fine mesh solution on coarse mesh->");
    OGProjection<double> ogProjection; ogProjection.project_global(space, ref_sln, sln);

    // Time measurement.
    cpu_time.tick();

    // VTK output.
    if (VTK_VISUALIZATION)
    {
      // Output solution in VTK format.
      Views::Linearizer lin(FileExport);
      char* title = new char[100];
      sprintf(title, "sln-%d.vtk", as);
      lin.save_solution_vtk(sln, title, "Potential", false);
      Hermes::Mixins::Loggable::Static::info("Solution in VTK format saved to file %s.", title);

      // Output mesh and element orders in VTK format.
      Views::Orderizer ord;
      sprintf(title, "ord-%d.vtk", as);
      ord.save_orders_vtk(space, title);
      Hermes::Mixins::Loggable::Static::info("Element orders in VTK format saved to file %s.", title);
    }

    // View the coarse mesh solution and polynomial orders.
    if (HERMES_VISUALIZATION)
    {
      sview.show(sln);
      oview.show(space);
    }

    // Skip visualization time.
    cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
    DefaultErrorCalculator<double, HERMES_H1_NORM> error_calculator(errorType, 1);
    error_calculator.calculate_errors(sln, ref_sln);
    double err_est_rel = error_calculator.get_total_error_squared() * 100.0;

    Adapt<double> adaptivity(space, &error_calculator);
    adaptivity.set_strategy(&stoppingCriterion);

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", space->get_num_dofs(), ref_space->get_num_dofs(), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    cpu_time.tick();
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    graph_dof.add_values(space->get_num_dofs(), err_est_rel);
    graph_dof.save("conv_dof_est.dat");

    // Skip the time spent to save the convergence graphs.
    cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

    // If err_est too large, adapt the mesh->
    if (err_est_rel < ERR_STOP)
      done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh->");
      done = adaptivity.adapt(&selector);

      // Increase the counter of performed adaptivity steps.
      if (done == false)
        as++;
    }
  } while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the fine mesh solution - final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(ref_sln);

  // Wait for all views to be closed.
  Views::View::wait();

  return 0;
}