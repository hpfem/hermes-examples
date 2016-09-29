#include "definitions.h"

// This example solves a simple linear wave equation by converting it
// into a system of two first-order equations in time. Time discretization
// is performed using arbitrary (explicit or implicit, low-order or higher-order)
// Runge-Kutta methods entered via their Butcher's tables.
// For a list of available R-K methods see the file hermes_common/tables.h.
//
// The function rk_time_step_newton() needs more optimisation, see a todo list at
// the beginning of file src/runge-kutta.h.
//
// PDE: \frac{1}{C_SQUARED}\frac{\partial^2 u}{\partial t^2} - \Delta u = 0,
// converted into
//
//      \frac{\partial u}{\partial t} = v,
//      \frac{\partial v}{\partial t} = C_SQUARED * \Delta u.
//
// BC:  u = 0 on the boundary,
//      v = 0 on the boundary (u = 0 => \partial u / \partial t = 0).
//
// IC:  smooth peak for u, zero for v.
//
// The following parameters can be changed:

// Initial polynomial degree of all elements.
const int P_INIT = 2;
// Refinement.
const int INIT_REF_NUM = 3;
// Time step.
const double INITIAL_TIME_STEP = 1e-2;
// If rel. temporal error is greater than this threshold, decrease time
// step size and repeat time step.
const double time_tol_upper = 100.0;
// If rel. temporal error is less than this threshold, increase time step
// but do not repeat time step (this might need further research).
const double time_tol_lower = 0.1;
// Timestep decrease ratio after unsuccessful nonlinear solve.
double time_step_dec = 0.5;
// Timestep increase ratio after successful nonlinear solve.
double time_step_inc = 2.0;
// Computation will stop if time step drops below this value.
double time_step_min = 1e-8;
// Final time.
const double T_FINAL = 20.;

// Choose one of the following time-integration methods, or define your own Butcher's table. The last number
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods:
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2,
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4,
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_5_4.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded,
//   Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded,
//   Implicit_DIRK_ISMAIL_7_45_embedded.
ButcherTableType butcher_table_type = Explicit_HEUN_EULER_2_12_embedded;

// Problem parameters.
// Square of wave speed.
const double C_SQUARED = 100;

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Refine towards boundary.
  mesh->refine_towards_boundary("Bdy", 1, true);

  // Refine once towards vertex #4.
  mesh->refine_towards_vertex(4, 1);

  // Refine all.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize solutions.
  MeshFunctionSharedPtr<double>  u_sln(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double>  v_sln(new ZeroSolution<double>(mesh));
  std::vector<MeshFunctionSharedPtr<double> > slns({ u_sln, v_sln });
  // - previous ones.
  MeshFunctionSharedPtr<double>  u_prev_sln(new CustomInitialConditionWave(mesh));
  MeshFunctionSharedPtr<double>  v_prev_sln(new ZeroSolution<double>(mesh));
  std::vector<MeshFunctionSharedPtr<double> > prev_slns({ u_prev_sln, v_prev_sln });
  // - time error functions.
  MeshFunctionSharedPtr<double> u_time_error_fn(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double> v_time_error_fn(new ZeroSolution<double>(mesh));
  std::vector<MeshFunctionSharedPtr<double> > time_error_fns({ u_time_error_fn, v_time_error_fn });

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf(new CustomWeakFormWave(INITIAL_TIME_STEP, C_SQUARED, u_prev_sln, v_prev_sln));

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  SpaceSharedPtr<double> u_space(new H1Space<double>(mesh, &bcs, P_INIT));
  SpaceSharedPtr<double> v_space(new H1Space<double>(mesh, &bcs, P_INIT));
  Hermes::Mixins::Loggable::Static::info("ndof = %d.", Space<double>::get_num_dofs({ u_space, v_space }));

  // Initialize views.
  ScalarView u_view("Solution u", new WinGeom(0, 0, 500, 400));
  u_view.fix_scale_width(50);
  ScalarView v_view("Solution v", new WinGeom(510, 0, 500, 400));
  v_view.fix_scale_width(50);
  ScalarView e_u_view("Temporal error u", new WinGeom(0, 410, 500, 400));
  e_u_view.fix_scale_width(50);
  ScalarView e_v_view("Temporal error v", new WinGeom(510, 410, 500, 400));
  e_v_view.fix_scale_width(50);
  // Visualize the solutions.
  char title[100];
  sprintf(title, "Initial Solution u");
  u_view.set_title(title);
  u_view.show(u_prev_sln);
  sprintf(title, "Initial Solution v");
  v_view.set_title(title);
  v_view.show(v_prev_sln);

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(wf, { u_space, v_space }, &bt);
  runge_kutta.set_verbose_output(true);

  // Time stepping loop.
  double current_time = INITIAL_TIME_STEP; int ts = 1; double time_step = INITIAL_TIME_STEP;
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    Hermes::Mixins::Loggable::Static::info("Runge-Kutta time step (t = %g s, time_step = %g s, stages: %d).",
      current_time, time_step, bt.get_size());
    try
    {
      runge_kutta.set_time(current_time);
      runge_kutta.set_time_step(time_step);
      // If we can, we will be interested in time error.
      if (bt.is_embedded())
        runge_kutta.rk_time_step_newton(prev_slns, slns, time_error_fns);
      else
        runge_kutta.rk_time_step_newton(prev_slns, slns);
    }
    catch (Exceptions::Exception& e)
    {
      Hermes::Mixins::Loggable::Static::info("Runge-Kutta time step failed, decreasing time step.");
      current_time -= time_step;
      time_step *= time_step_dec;
      continue;
      if (time_step < time_step_min)
        throw Hermes::Exceptions::Exception("Time step became too small.");

      e.print_msg();
    }

    // Time error handling
    if (bt.is_embedded())
    {
      // Show error functions.
      e_u_view.show(u_time_error_fn);
      e_v_view.show(v_time_error_fn);

      // Calculate relative time stepping error and decide whether the
      // time step can be accepted. If not, then the time step size is
      // reduced and the entire time step repeated. If yes, then another
      // check is run, and if the relative error is very low, time step
      // is increased.
      DefaultNormCalculator<double, HERMES_H1_NORM> normCalculator(2);
      normCalculator.calculate_norms(time_error_fns);
      double rel_err_time = normCalculator.get_total_norm_squared() * 100.;

      Hermes::Mixins::Loggable::Static::info("rel_err_time = %g%%", rel_err_time);
      if (rel_err_time > time_tol_upper) {
        Hermes::Mixins::Loggable::Static::info("rel_err_time above upper limit %g%% -> decreasing time step from %g to %g days and repeating time step.",
          time_tol_upper, time_step, time_step * time_step_dec);
        time_step *= time_step_dec;
        continue;
      }
      if (rel_err_time < time_tol_lower) {
        Hermes::Mixins::Loggable::Static::info("rel_err_time = below lower limit %g%% -> increasing time step from %g to %g days",
          time_tol_lower, time_step, time_step * time_step_inc);
        time_step *= time_step_inc;
      }
    }

    // Advance the solutions.
    u_prev_sln->copy(u_sln);
    v_prev_sln->copy(v_sln);

    // Visualize the solutions.
    char title[100];
    sprintf(title, "Solution u, t = %g", current_time);
    u_view.set_title(title);
    u_view.show(u_sln);
    sprintf(title, "Solution v, t = %g", current_time);
    v_view.set_title(title);
    v_view.show(v_sln);

    // Update time.
    current_time += time_step;
  } while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();

  return 0;
}