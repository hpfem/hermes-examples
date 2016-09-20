#include "definitions.h"

// This example solves a time-domain resonator problem for the Maxwell's equation.
// It is very similar to resonator-time-domain-I but B is eliminated from the
// equations, thus converting the first-order system into one second -order
// equation in time. The second-order equation in time is decomposed back into
// a first-order system in time in the standard way (see example wave-1). Time
// discretization is performed using arbitrary Runge-Kutta methods entered via
// their Butcher's tables. For a list of available R-K methods see the file
// hermes_common/tables.h.
//
// The function rk_time_step_newton() needs more optimisation, see a todo list at
// the beginning of file src/runge-kutta.h.
//
// PDE: \frac{1}{SPEED_OF_LIGHT**2}\frac{\partial^2 E}{\partial t^2} + curl curl E = 0,
// converted into
//
//      \frac{\partial E}{\partial t} = F,
//      \frac{\partial F}{\partial t} = -SPEED_OF_LIGHT**2 * curl curl E.
//
// Domain: Square (-pi/2, pi/2) x (-pi/2, pi/2)... See mesh file domain.mesh->
//
// BC:  E \times \nu = 0 on the boundary (perfect conductor),
//      F \times \nu = 0 on the boundary (E \times \nu = 0 => \partial E / \partial t \times \nu = 0).
//
// IC:  Prescribed wave for E, zero for F.
//
// The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 6;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// Time step.
const double time_step = 0.05;
// Final time.
const double T_FINAL = 35.0;
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-5;
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;

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
ButcherTableType butcher_table = Implicit_RK_1;
//ButcherTableType butcher_table = Implicit_SDIRK_2_2;

// Problem parameters.
// Square of wave speed.
const double C_SQUARED = 1;

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table);
  if (bt.is_explicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Initialize solutions.
  MeshFunctionSharedPtr<double> E_time_prev(new CustomInitialConditionWave(mesh));
  MeshFunctionSharedPtr<double>  F_time_prev(new ZeroSolutionVector<double>(mesh));

  std::vector<MeshFunctionSharedPtr<double> > slns_time_prev({ E_time_prev, F_time_prev });

  MeshFunctionSharedPtr<double> E_time_new(new Solution<double>);

  MeshFunctionSharedPtr<double> F_time_new(new Solution<double>);

  std::vector<MeshFunctionSharedPtr<double> > slns_time_new({ E_time_new, F_time_new });

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf(new CustomWeakFormWaveRK(C_SQUARED));

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential("Perfect conductor", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  SpaceSharedPtr<double> E_space(new HcurlSpace<double>(mesh, &bcs, P_INIT));
  SpaceSharedPtr<double> F_space(new HcurlSpace<double>(mesh, &bcs, P_INIT));
  std::vector<SpaceSharedPtr<double> > spaces({ E_space, F_space });

  int ndof = HcurlSpace<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof = %d.", ndof);

  // Initialize views.
  ScalarView E1_view("Solution E1", new WinGeom(0, 0, 400, 350));
  E1_view.fix_scale_width(50);
  ScalarView E2_view("Solution E2", new WinGeom(410, 0, 400, 350));
  E2_view.fix_scale_width(50);
  ScalarView F1_view("Solution F1", new WinGeom(0, 410, 400, 350));
  F1_view.fix_scale_width(50);
  ScalarView F2_view("Solution F2", new WinGeom(410, 410, 400, 350));
  F2_view.fix_scale_width(50);

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(wf, spaces, &bt);

  // Time stepping loop.
  double current_time = 0; int ts = 1;
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    Hermes::Mixins::Loggable::Static::info("Runge-Kutta time step (t = %g s, time_step = %g s, stages: %d).",
      current_time, time_step, bt.get_size());
    try
    {
      runge_kutta.set_time(current_time);
      runge_kutta.set_time_step(time_step);
      runge_kutta.set_newton_max_allowed_iterations(NEWTON_MAX_ITER);
      runge_kutta.set_newton_tolerance(NEWTON_TOL);
      runge_kutta.rk_time_step_newton(slns_time_prev, slns_time_new);
    }
    catch (Exceptions::Exception& e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Runge-Kutta time step failed");
    }

    // Visualize the solutions.
    char title[100];
    sprintf(title, "E1, t = %g", current_time + time_step);
    E1_view.set_title(title);
    E1_view.show(E_time_new, H2D_FN_VAL_0);
    sprintf(title, "E2, t = %g", current_time + time_step);
    E2_view.set_title(title);
    E2_view.show(E_time_new, H2D_FN_VAL_1);

    sprintf(title, "F1, t = %g", current_time + time_step);
    F1_view.set_title(title);
    F1_view.show(F_time_new, H2D_FN_VAL_0);
    sprintf(title, "F2, t = %g", current_time + time_step);
    F2_view.set_title(title);
    F2_view.show(F_time_new, H2D_FN_VAL_1);

    //View::wait();

    // Update solutions.
    E_time_prev->copy(E_time_new);
    F_time_prev->copy(F_time_new);

    // Update time.
    current_time += time_step;
  } while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();

  return 0;
}