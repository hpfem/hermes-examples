#include "definitions.h"

using namespace RefinementSelectors;

//  This examples is an analogy to the 2D benchmark moving-front-space-adapt.
//  It has an exact solution that contains a moving front with arbitrary
//  steepness S. The mesh changes dynamically in time. Arbitrary RK method
//  can be used for time integration.
//
//  PDE: time-dependent heat transfer equation, du/dt = u'' + f.
//
//  Domain: Interval (0, 10).
//
//  BC: Zero Dirichlet.
//
//  IC: Zero.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// Initial polynomial degree of mesh elements.
const int P_INIT = 1;
// Time step.
const double time_step = 0.1;
// Time interval length.
const double T_FINAL = 10.0;

// Adaptivity
// Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_FREQ = 1;
// 1... mesh reset to basemesh and poly degrees to P_INIT.
// 2... one ref. layer shaved off, poly degrees reset to P_INIT.
// 3... one ref. layer shaved off, poly degrees decreased by one.
const int UNREF_METHOD = 3;
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

// Newton's method
// Stopping criterion for Newton on fine mesh->
const double NEWTON_TOL = 1e-5;
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 20;

// Choose one of the following time-integration methods, or define your own Butcher's table. The last number
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods:
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2,
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4,
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_5_4.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
//   Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded,
//   Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded,
//   Implicit_DIRK_ISMAIL_7_45_embedded.
ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

// Problem parameters.
double x_0 = 0.0;
double x_1 = 10.0;
double y_0 = -5.0;
double y_1 = 5.0;
// Steepness of the moving front.
double s = 20.0;
double c = 1000.0;

// Current time.
double current_time = 0.0;

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh->
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH1DXML mloader;
  try
  {
    mloader.load("domain.xml", basemesh);
  }
  catch (Hermes::Exceptions::MeshLoadFailureException& e)
  {
    e.print_msg();
    return -1;
  }

  // Perform initial mesh refinements.
  int refinement_type = 2;                        // Split elements vertically.
  for (int i = 0; i < INIT_REF_NUM; i++) basemesh->refine_all_elements(refinement_type, true);
  mesh->copy(basemesh);

  // Exact solution.
  CustomExactSolution exact_sln(mesh, x_0, x_1, y_0, y_1, &current_time, s, c);

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_essential({"Left", "Right"}, 0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof_coarse = space->get_num_dofs();

  // Initialize the weak formulation
  CustomFunction f(x_0, x_1, y_0, y_1, s, c);
  CustomWeakFormPoisson wf(new Hermes::Hermes1DFunction<double>(-1.0), &f);

  // Previous and next time level solution.
  MeshFunctionSharedPtr<double>  sln_time_prev(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double> sln_time_new(new Solution<double>(mesh));

  // Create a refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  // Visualize initial condition.
  char title[100];
  Views::ScalarView sview("Initial condition", new Views::WinGeom(0, 0, 1200, 200));
  Views::OrderView oview("Initial mesh", new Views::WinGeom(0, 260, 1200, 200));
  sview.show(sln_time_prev);
  oview.show(space);

  // Graph for dof history.
  SimpleGraph dof_history_graph;

  // Time stepping loop.
  int ts = 1;
  do
  {
    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0)
    {
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
      case 1: mesh->copy(basemesh);
        space->set_uniform_order(P_INIT);
        break;
      case 2: mesh->unrefine_all_elements();
        space->set_uniform_order(P_INIT);
        break;
      case 3: mesh->unrefine_all_elements();
        space->adjust_element_order(-1, -1, P_INIT, P_INIT);
        break;
      default: throw Hermes::Exceptions::Exception("Wrong global derefinement method.");
      }

      space->assign_dofs();
      ndof_coarse = space->get_num_dofs();
    }

    // Spatial adaptivity loop. Note: sln_time_prev must not be changed
    // during spatial adaptivity.
    bool done = false; int as = 1;
    double err_est;
    do {
      Hermes::Mixins::Loggable::Static::info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh and setup reference space.
      // FIXME: This should be increase in the x-direction only.
      int order_increase = 1;
      // FIXME: This should be '2' but that leads to a segfault.
      int refinement_type = 0;
      Mesh::ReferenceMeshCreator refMeshCreator(mesh);
      MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

      Space<double>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh);
      SpaceSharedPtr<double> ref_space = refSpaceCreator.create_ref_space();
      int ndof_ref = ref_space->get_num_dofs();

      // Initialize Runge-Kutta time stepping.
      RungeKutta<double> runge_kutta(wf, ref_space, &bt);

      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      Hermes::Mixins::Loggable::Static::info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).",
        current_time, time_step, bt.get_size());
      bool freeze_jacobian = true;
      bool block_diagonal_jacobian = false;
      bool verbose = true;

      try
      {
        runge_kutta.set_time(current_time);
        runge_kutta.set_time_step(time_step);
        runge_kutta.set_max_allowed_iterations(NEWTON_MAX_ITER);
        runge_kutta.set_tolerance(NEWTON_TOL);
        runge_kutta.rk_time_step_newton(sln_time_prev, sln_time_new);
      }
      catch (Exceptions::Exception& e)
      {
        e.print_msg();
        throw Hermes::Exceptions::Exception("Runge-Kutta time step failed");
      }

      // Project the fine mesh solution onto the coarse mesh->
      MeshFunctionSharedPtr<double> sln_coarse(new Solution<double>);
      Hermes::Mixins::Loggable::Static::info("Projecting fine mesh solution on coarse mesh for error estimation.");
      OGProjection<double> ogProjection; ogProjection.project_global(space, sln_time_new, sln_coarse);

      // Calculate element errors and total error estimate.
      Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
      adaptivity.set_space(space);
      errorCalculator.calculate_errors(sln_coarse, sln_time_new);
      double err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

      // Report results.
      Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_ref: %d, err_est_rel: %g%%",
        space->get_num_dofs(), ref_space->get_num_dofs(), err_est_rel_total);

      // If err_est too large, adapt the mesh->
      if (err_est_rel_total < ERR_STOP) done = true;
      else
      {
        Hermes::Mixins::Loggable::Static::info("Adapting the coarse mesh->");
        done = adaptivity.adapt(&selector);

        // Increase the counter of performed adaptivity steps.
        as++;
      }
    } while (done == false);

    // Visualize the solution and mesh->
    char title[100];
    sprintf(title, "Solution, time %g", current_time);
    sview.set_title(title);
    sview.show_mesh(false);
    sview.show(sln_time_new);
    sprintf(title, "Mesh, time %g", current_time);
    oview.set_title(title);
    oview.show(space);

    // Copy last reference solution into sln_time_prev->
    sln_time_prev->copy(sln_time_new);

    dof_history_graph.add_values(current_time, space->get_num_dofs());
    dof_history_graph.save("dof_history.dat");

    // Increase current time and counter of time steps.
    current_time += time_step;
    ts++;
  } while (current_time < T_FINAL);

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}