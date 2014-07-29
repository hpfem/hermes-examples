#include "definitions.h"

//  This example solves a linear advection diffusion problem using optional
//  variational multiscale stabilization. To use the stabilization, you must
//  uncomment the definition of H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h
//  and rebuild this example. Note that in our experience, the stabilization
//  does only work for linear elements. Nevertheless, we were able to solve the
//  problem without stabilization using adaptive hp-FEM.
//
//  PDE: div(bu - \epsilon \nabla u) = 0 where b = (b1, b2) is a constant vector.
//
//  Domain: Square (0, 1)x(0, 1).
//
//  BC:  Dirichlet, see the function double essential_bc_values() below.
//
//  The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;
// Stabilization on/off (assumes that H2D_SECOND_DERIVATIVES_ENABLED is defined).
const bool STABILIZATION_ON = false;
// Shock capturing on/off.
const bool SHOCK_CAPTURING_ON = true;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;
// Number of initial uniform mesh refinements in the boundary layer region.
const int INIT_REF_NUM_BDY = 1;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;
// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;

// Problem parameters.
// Diffusivity.
const double EPSILON = 0.01;
// Advection direction, div(B) = 0.
const double B1 = 1., B2 = 1.;

int main(int argc, char* argv[])
{
  // Load the mesh->
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square_quad.mesh", mesh);
  // mloader.load("square_tri.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();
  mesh->refine_towards_boundary("Layer", INIT_REF_NUM_BDY);

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf(new WeakFormLinearAdvectionDiffusion(STABILIZATION_ON, SHOCK_CAPTURING_ON, B1, B2, EPSILON));

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_rest("Rest", 1.0);
  EssentialBCNonConst bc_layer("Layer");

  EssentialBCs<double> bcs({&bc_rest, &bc_layer});

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  WinGeom* sln_win_geom = new WinGeom(0, 0, 440, 350);
  WinGeom* mesh_win_geom = new WinGeom(450, 0, 400, 350);

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>), ref_sln(new Solution<double>);

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  OrderView  oview("Polynomial orders", new WinGeom(450, 0, 400, 350));

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator refMeshCreator(mesh);
    MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh);
    SpaceSharedPtr<double> ref_space = refSpaceCreator.create_ref_space();

    // Assemble the reference problem.
    Hermes::Mixins::Loggable::Static::info("Solving on reference mesh->");
    LinearSolver<double> solver(wf, ref_space);

    // Time measurement.
    cpu_time.tick();

    // Solve the linear system of the reference problem.
    // If successful, obtain the solution.
    solver.solve();
    Solution<double>::vector_to_solution(solver.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh->
    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh->");
    OGProjection<double> ogProjection; ogProjection.project_global(space, ref_sln, sln);

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution and polynomial orders.
    sview.show(sln);
    oview.show(space);

    // Skip visualization time.
    cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
    adaptivity.set_space(space);
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      Space<double>::get_num_dofs(space), Space<double>::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space<double>::get_num_dofs(space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh->
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh->");
      done = adaptivity.adapt(&selector);

      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
  } while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(ref_sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}