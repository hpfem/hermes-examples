#include "definitions.h"

//  This example shows how to combine automatic adaptivity with the Newton's
//  method for a nonlinear complex-valued time-dependent PDE (the Gross-Pitaevski
//  equation describing the behavior of Einstein-Bose quantum gases)
//  discretized implicitly in time (via implicit Euler or Crank-Nicolson).
//
//  PDE: non-stationary complex Gross-Pitaevski equation
//  describing resonances in Bose-Einstein condensates.
//
//  ih \partial \psi/\partial t = -h^2/(2m) \Delta \psi +
//  g \psi |\psi|^2 + 1/2 m \omega^2 (x^2 + y^2) \psi.
//
//  Domain: square (-1, 1)^2.
//
//  BC:  homogeneous Dirichlet everywhere on the boundary.
//
//  Time-stepping: either implicit Euler or Crank-Nicolson.
//
//  The following parameters can be changed:

// Number of initial uniform refinements.
const int INIT_REF_NUM = 2;
// Initial polynomial degree.
const int P_INIT = 2;
// Time interval length.
const double T_FINAL = 2.0;
// Time step.
double time_step = 0.005;

// Adaptivity.
// Every UNREF_FREQ time step the mesh is unrefined.
const int UNREF_FREQ = 1;
// 1... mesh reset to basemesh and poly degrees to P_INIT.
// 2... one ref. layer shaved off, poly degrees reset to P_INIT.
// 3... one ref. layer shaved off, poly degrees decreased by one.
const int UNREF_METHOD = 3;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;
// Error calculation & adaptivity.
DefaultErrorCalculator<::complex, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<::complex> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<::complex> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;

// Temporal adaptivity.
// This flag decides whether adaptive time stepping will be done.
// The methods for the adaptive and fixed-step versions are set
// below. An embedded method must be used with adaptive time stepping.
bool ADAPTIVE_TIME_STEP_ON = true;
// If rel. temporal error is greater than this threshold, decrease time
// step size and repeat time step.
const double TIME_ERR_TOL_UPPER = 1.0;
// If rel. temporal error is less than this threshold, increase time step
// but do not repeat time step (this might need further research).
const double TIME_ERR_TOL_LOWER = 0.1;
// Time step increase ratio (applied when rel. temporal error is too small).
const double TIME_STEP_INC_RATIO = 1.1;
// Time step decrease ratio (applied when rel. temporal error is too large).
const double TIME_STEP_DEC_RATIO = 0.8;

// Newton's method.
// Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_COARSE = 0.01;
// Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL_FINE = 0.05;
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 50;

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
ButcherTableType butcher_table_type = Implicit_SDIRK_CASH_3_23_embedded;

// Problem parameters.
// Planck constant 6.626068e-34.
const double h = 1;
// Mass of boson.
const double m = 1;
// Coupling constant.
const double g = 1;
// Frequency.
const double omega = 1;

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Turn off adaptive time stepping if R-K method is not embedded.
  if (bt.is_embedded() == false && ADAPTIVE_TIME_STEP_ON == true) {
    Hermes::Mixins::Loggable::Static::warn("R-K method not embedded, turning off adaptive time stepping.");
    ADAPTIVE_TIME_STEP_ON = false;
  }

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", basemesh);
  mesh->copy(basemesh);

  // Initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Convert initial condition into a Solution<::complex>.
  MeshFunctionSharedPtr<::complex> psi_time_prev(new CustomInitialCondition(mesh));

  // Initialize the weak formulation.
  double current_time = 0;

  // Initialize weak formulation.
  WeakFormSharedPtr<::complex> wf(new CustomWeakFormGPRK(h, m, g, omega));

  // Initialize boundary conditions.
  DefaultEssentialBCConst<::complex> bc_essential("Bdy", 0.0);
  EssentialBCs<::complex> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<::complex> space(new H1Space<::complex>(mesh, &bcs, P_INIT));
  adaptivity.set_space(space);
  int ndof = space->get_num_dofs();
  Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

  // Create a refinement selector.
  H1ProjBasedSelector<::complex> selector(CAND_LIST);

  // Visualize initial condition.
  char title[100];

  ScalarView sview_real("Initial condition - real part", new WinGeom(0, 0, 600, 500));
  ScalarView sview_imag("Initial condition - imaginary part", new WinGeom(610, 0, 600, 500));

  sview_real.show_mesh(false);
  sview_imag.show_mesh(false);
  sview_real.fix_scale_width(50);
  sview_imag.fix_scale_width(50);
  OrderView ord_view("Initial mesh", new WinGeom(445, 0, 440, 350));
  ord_view.fix_scale_width(50);
  ScalarView time_error_view("Temporal error", new WinGeom(0, 400, 440, 350));
  time_error_view.fix_scale_width(50);
  time_error_view.fix_scale_width(60);
  ScalarView space_error_view("Spatial error", new WinGeom(445, 400, 440, 350));
  space_error_view.fix_scale_width(50);
  MeshFunctionSharedPtr<double> real(new RealFilter(psi_time_prev));

  MeshFunctionSharedPtr<double> imag(new ImagFilter(psi_time_prev));

  sview_real.show(real);
  sview_imag.show(imag);
  ord_view.show(space);

  // Graph for time step history.
  SimpleGraph time_step_graph;
  if (ADAPTIVE_TIME_STEP_ON) Hermes::Mixins::Loggable::Static::info("Time step history will be saved to file time_step_history.dat.");

  // Time stepping:
  int num_time_steps = (int)(T_FINAL / time_step + 0.5);
  for (int ts = 1; ts <= num_time_steps; ts++)
    // Time stepping loop.
    double current_time = 0.0; int ts = 1;
  do
  {
    Hermes::Mixins::Loggable::Static::info("Begin time step %d.", ts);
    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0)
    {
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
      case 1: mesh->copy(basemesh);
        space->set_uniform_order(P_INIT);
        break;
      case 2: space->unrefine_all_mesh_elements();
        space->set_uniform_order(P_INIT);
        break;
      case 3: space->unrefine_all_mesh_elements();
        space->adjust_element_order(-1, -1, P_INIT, P_INIT);
        break;
      default: throw Hermes::Exceptions::Exception("Wrong global derefinement method.");
      }

      ndof = Space<::complex>::get_num_dofs(space);
    }
    Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

    // Spatial adaptivity loop. Note: psi_time_prev must not be
    // changed during spatial adaptivity.
    MeshFunctionSharedPtr<::complex> ref_sln(new Solution<::complex>());
    MeshFunctionSharedPtr<::complex> time_error_fn(new Solution<::complex>);
    bool done = false;
    int as = 1;
    double err_est;
    do {
      // Construct globally refined reference mesh and setup reference space.
      Mesh::ReferenceMeshCreator refMeshCreator(mesh);
      MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

      Space<::complex>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh);
      SpaceSharedPtr<::complex> ref_space = refSpaceCreator.create_ref_space();

      RungeKutta<::complex> runge_kutta(wf, ref_space, &bt);

      // Runge-Kutta step on the fine mesh.
      Hermes::Mixins::Loggable::Static::info("Runge-Kutta time step on fine mesh (t = %g s, time step = %g s, stages: %d).",
        current_time, time_step, bt.get_size());
      bool verbose = true;

      try
      {
        runge_kutta.set_time(current_time);
        runge_kutta.set_time_step(time_step);
        runge_kutta.set_max_allowed_iterations(NEWTON_MAX_ITER);
        runge_kutta.set_tolerance(NEWTON_TOL_FINE);
        runge_kutta.rk_time_step_newton(psi_time_prev, ref_sln, time_error_fn);
      }
      catch (Exceptions::Exception& e)
      {
        e.print_msg();
        throw Hermes::Exceptions::Exception("Runge-Kutta time step failed");
      }

      /* If ADAPTIVE_TIME_STEP_ON == true, estimate temporal error.
      If too large or too small, then adjust it and restart the time step. */

      double rel_err_time = 0;
      if (bt.is_embedded() == true) {
        Hermes::Mixins::Loggable::Static::info("Calculating temporal error estimate.");

        // Show temporal error.
        char title[100];
        sprintf(title, "Temporal error est, spatial adaptivity step %d", as);
        time_error_view.set_title(title);
        time_error_view.show_mesh(false);
        MeshFunctionSharedPtr<double> abs_time(new RealFilter(time_error_fn));

        MeshFunctionSharedPtr<double> abs_tef(new AbsFilter(abs_time));

        time_error_view.show(abs_tef);

        DefaultNormCalculator<::complex, HERMES_H1_NORM> normCalculator(1);
        rel_err_time = 100. * normCalculator.calculate_norm(time_error_fn) / normCalculator.calculate_norm(ref_sln);
        if (ADAPTIVE_TIME_STEP_ON == false)
          Hermes::Mixins::Loggable::Static::info("rel_err_time: %g%%", rel_err_time);
      }

      if (ADAPTIVE_TIME_STEP_ON) {
        if (rel_err_time > TIME_ERR_TOL_UPPER) {
          Hermes::Mixins::Loggable::Static::info("rel_err_time %g%% is above upper limit %g%%", rel_err_time, TIME_ERR_TOL_UPPER);
          Hermes::Mixins::Loggable::Static::info("Decreasing time step from %g to %g s and restarting time step.",
            time_step, time_step * TIME_STEP_DEC_RATIO);
          time_step *= TIME_STEP_DEC_RATIO;

          continue;
        }
        else if (rel_err_time < TIME_ERR_TOL_LOWER) {
          Hermes::Mixins::Loggable::Static::info("rel_err_time = %g%% is below lower limit %g%%", rel_err_time, TIME_ERR_TOL_LOWER);
          Hermes::Mixins::Loggable::Static::info("Increasing time step from %g to %g s.", time_step, time_step * TIME_STEP_INC_RATIO);
          time_step *= TIME_STEP_INC_RATIO;

          continue;
        }
        else {
          Hermes::Mixins::Loggable::Static::info("rel_err_time = %g%% is in acceptable interval (%g%%, %g%%)",
            rel_err_time, TIME_ERR_TOL_LOWER, TIME_ERR_TOL_UPPER);
        }

        // Add entry to time step history graph.
        time_step_graph.add_values(current_time, time_step);
        time_step_graph.save("time_step_history.dat");
      }

      /* Estimate spatial errors and perform mesh refinement */

      Hermes::Mixins::Loggable::Static::info("Spatial adaptivity step %d.", as);

      // Project the fine mesh solution onto the coarse mesh.
      MeshFunctionSharedPtr<::complex> sln(new Solution<::complex>);
      Hermes::Mixins::Loggable::Static::info("Projecting fine mesh solution on coarse mesh for error estimation.");
      OGProjection<::complex> ogProjection; ogProjection.project_global(space, ref_sln, sln);

      // Show spatial error.
      sprintf(title, "Spatial error est, spatial adaptivity step %d", as);
      MeshFunctionSharedPtr<::complex> space_error_fn(new DiffFilter<::complex>(std::vector<MeshFunctionSharedPtr<::complex> >({ ref_sln, sln })));

      space_error_view.set_title(title);
      space_error_view.show_mesh(false);

      MeshFunctionSharedPtr<double> abs_space(new RealFilter(space_error_fn));
      MeshFunctionSharedPtr<double> abs_sef(new AbsFilter(abs_space));

      space_error_view.show(abs_sef);

      // Calculate element errors and spatial error estimate.
      Hermes::Mixins::Loggable::Static::info("Calculating spatial error estimate.");
      errorCalculator.calculate_errors(sln, ref_sln);
      double err_rel_space = errorCalculator.get_total_error_squared() * 100.;

      // Report results.
      Hermes::Mixins::Loggable::Static::info("ndof: %d, ref_ndof: %d, err_rel_space: %g%%",
        Space<::complex>::get_num_dofs(space), Space<::complex>::get_num_dofs(ref_space), err_rel_space);

      // If err_est too large, adapt the mesh.
      if (err_rel_space < ERR_STOP)
        done = true;
      else
      {
        Hermes::Mixins::Loggable::Static::info("Adapting the coarse mesh.");
        done = adaptivity.adapt(&selector);

        // Increase the counter of performed adaptivity steps.
        as++;
      }
    } while (done == false);

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution - real part, Time %3.2f s", current_time);
    sview_real.set_title(title);
    sprintf(title, "Solution - imaginary part, Time %3.2f s", current_time);
    sview_imag.set_title(title);
    MeshFunctionSharedPtr<double> real(new RealFilter(ref_sln));
    MeshFunctionSharedPtr<double> imag(new ImagFilter(ref_sln));
    sview_real.show(real);
    sview_imag.show(imag);
    sprintf(title, "Mesh, time %g s", current_time);
    ord_view.set_title(title);
    ord_view.show(space);

    // Copy last reference solution into psi_time_prev.
    psi_time_prev->copy(ref_sln);

    // Increase current time and counter of time steps.
    current_time += time_step;
    ts++;
  } while (current_time < T_FINAL);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}