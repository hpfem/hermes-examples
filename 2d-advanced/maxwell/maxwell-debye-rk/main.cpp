#include "definitions.h"

// This example is a simple test case for the Debye-Maxwell model solved in terms of
// E, H and P. Here E is electric field (vector), H magnetic field (scalar), and P
// electric polarization (vector). The example comes with a known exact solution.
// Time discretization is performed using an arbitrary Runge-Kutta method.
//
// PDE system:
//
//      \frac{partial H}{\partial t} + \frac{1}{\mu_0} \curl E = 0,
//      \frac{\partial E}{\partial t} - \frac{1}{\epsilon_0 \epsilon_\infty} \curl H
//          + \frac{\epsilon_q - 1}{\tau}E - \frac{1}{\tau} P = 0,
//      \frac{\partial P}{\partial t} - \frac{(\epsilon_q - 1)\epsilon_0 \epsilon_\infty}{\tau} E
//          + \frac{1}{\tau} P = 0.
//
// Domain: Square (0, 1) x (0, 1).
//
// BC:  Perfect conductor for E and P.
//
// IC:  Prescribed functions for E, H and P.
//
// The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;
// Time step.
const double time_step = 0.00001;
// Final time.
const double T_FINAL = 35.0;
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-4;
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
// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 5;

// Number of mesh refinements between two unrefinements.
// The mesh is not unrefined unless there has been a refinement since
// last unrefinement.
int REFINEMENT_COUNT = 0;

// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies (see below).
const double THRESHOLD = 0.5;
// Error calculation & adaptivity.
class CustomErrorCalculator : public ErrorCalculator < double >
{
public:
  // Two first components are in the H1 space - we can use the classic class for that, for the last component, we will manually add the L2 norm for pressure.
  CustomErrorCalculator(CalculatedErrorType errorType) : ErrorCalculator<double>(errorType)
  {
    this->add_error_form(new DefaultNormFormVol<double>(0, 0, HERMES_HCURL_NORM));
    this->add_error_form(new DefaultNormFormVol<double>(1, 1, HERMES_H1_NORM));
    this->add_error_form(new DefaultNormFormVol<double>(2, 2, HERMES_HCURL_NORM));
  }
} errorCalculator(RelativeErrorToGlobalNorm);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;

// Problem parameters.
// Permeability of free space->
const double MU_0 = 1.0;
// Permittivity of free space->
const double EPS_0 = 1.0;
// Permittivity at infinite frequency.
const double EPS_INF = 1.0;
// Permittivity at zero frequency.
const double EPS_S = 2.0;
// EPS_Q.
const double EPS_Q = EPS_S / EPS_INF;
// Relaxation time of the medium.
const double TAU = 1.0;
// Angular frequency. Depends on wave number K. Must satisfy:
// omega^3 - 2 omega^2 + K^2 M_PI^2 omega - K^2 M_PI^2 = 0.
// WARNING; Choosing wrong omega may lead to K**2 < 0.
const double OMEGA = 1.5;
// Wave vector direction (will be normalized to be compatible
// with omega).
double K_x = 1.0;
double K_y = 1.0;

int main(int argc, char* argv[])
{
  try
  {
    // Sanity check for omega.
    double K_squared = Hermes::sqr(OMEGA / M_PI) * (OMEGA - 2) / (1 - OMEGA);
    if (K_squared <= 0) throw Hermes::Exceptions::Exception("Wrong choice of omega, K_squared < 0!");
    double K_norm_coeff = std::sqrt(K_squared) / std::sqrt(Hermes::sqr(K_x) + Hermes::sqr(K_y));
    Hermes::Mixins::Loggable::Static::info("Wave number K = %g", std::sqrt(K_squared));
    K_x *= K_norm_coeff;
    K_y *= K_norm_coeff;

    // Wave number.
    double K = std::sqrt(Hermes::sqr(K_x) + Hermes::sqr(K_y));

    // Choose a Butcher's table or define your own.
    ButcherTable bt(butcher_table);
    if (bt.is_explicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage explicit R-K method.", bt.get_size());
    if (bt.is_diagonally_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
    if (bt.is_fully_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage fully implicit R-K method.", bt.get_size());

    // Load the mesh.
    MeshSharedPtr E_mesh(new Mesh), H_mesh(new Mesh), P_mesh(new Mesh);
    MeshReaderH2D mloader;
    mloader.load("domain.mesh", E_mesh);
    mloader.load("domain.mesh", H_mesh);
    mloader.load("domain.mesh", P_mesh);

    // Perform initial mesh refinemets.
    for (int i = 0; i < INIT_REF_NUM; i++)
    {
      E_mesh->refine_all_elements();
      H_mesh->refine_all_elements();
      P_mesh->refine_all_elements();
    }

    // Initialize solutions.
    double current_time = 0;
    MeshFunctionSharedPtr<double> E_time_prev(new CustomInitialConditionE(E_mesh, current_time, OMEGA, K_x, K_y));
    MeshFunctionSharedPtr<double> H_time_prev(new CustomInitialConditionH(H_mesh, current_time, OMEGA, K_x, K_y));
    MeshFunctionSharedPtr<double> P_time_prev(new CustomInitialConditionP(P_mesh, current_time, OMEGA, K_x, K_y));
    std::vector<MeshFunctionSharedPtr<double> > slns_time_prev({ E_time_prev, H_time_prev, P_time_prev });
    MeshFunctionSharedPtr<double> E_time_new(new Solution<double>(E_mesh)), H_time_new(new Solution<double>(H_mesh)), P_time_new(new Solution<double>(P_mesh));
    MeshFunctionSharedPtr<double> E_time_new_coarse(new Solution<double>(E_mesh)), H_time_new_coarse(new Solution<double>(H_mesh)), P_time_new_coarse(new Solution<double>(P_mesh));
    std::vector<MeshFunctionSharedPtr<double> > slns_time_new({ E_time_new, H_time_new, P_time_new });

    // Initialize the weak formulation.
    WeakFormSharedPtr<double> wf(new CustomWeakFormMD(OMEGA, K_x, K_y, MU_0, EPS_0, EPS_INF, EPS_Q, TAU));

    // Initialize boundary conditions
    DefaultEssentialBCConst<double> bc_essential("Bdy", 0.0);
    EssentialBCs<double> bcs(&bc_essential);

    SpaceSharedPtr<double> E_space(new HcurlSpace<double>(E_mesh, &bcs, P_INIT));
    SpaceSharedPtr<double> H_space(new H1Space<double>(H_mesh, NULL, P_INIT));
    //L2Space<double> H_space(mesh, P_INIT));
    SpaceSharedPtr<double> P_space(new HcurlSpace<double>(P_mesh, &bcs, P_INIT));

    std::vector<SpaceSharedPtr<double> > spaces = std::vector<SpaceSharedPtr<double> >({ E_space, H_space, P_space });

    // Initialize views.
    ScalarView E1_view("Solution E1", new WinGeom(0, 0, 400, 350));
    E1_view.fix_scale_width(50);
    ScalarView E2_view("Solution E2", new WinGeom(410, 0, 400, 350));
    E2_view.fix_scale_width(50);
    ScalarView H_view("Solution H", new WinGeom(0, 410, 400, 350));
    H_view.fix_scale_width(50);
    ScalarView P1_view("Solution P1", new WinGeom(410, 410, 400, 350));
    P1_view.fix_scale_width(50);
    ScalarView P2_view("Solution P2", new WinGeom(820, 410, 400, 350));
    P2_view.fix_scale_width(50);

    // Visualize initial conditions.
    char title[100];
    sprintf(title, "E1 - Initial Condition");
    E1_view.set_title(title);
    E1_view.show(E_time_prev, H2D_FN_VAL_0);
    sprintf(title, "E2 - Initial Condition");
    E2_view.set_title(title);
    E2_view.show(E_time_prev, H2D_FN_VAL_1);

    sprintf(title, "H - Initial Condition");
    H_view.set_title(title);
    H_view.show(H_time_prev);

    sprintf(title, "P1 - Initial Condition");
    P1_view.set_title(title);
    P1_view.show(P_time_prev, H2D_FN_VAL_0);
    sprintf(title, "P2 - Initial Condition");
    P2_view.set_title(title);
    P2_view.show(P_time_prev, H2D_FN_VAL_1);

    // Initialize Runge-Kutta time stepping.
    RungeKutta<double> runge_kutta(wf, spaces, &bt);
    runge_kutta.set_max_allowed_iterations(NEWTON_MAX_ITER);
    runge_kutta.set_tolerance(NEWTON_TOL);
    runge_kutta.set_verbose_output(true);

    // Initialize refinement selector.
    H1ProjBasedSelector<double> H1selector(CAND_LIST);
    HcurlProjBasedSelector<double> HcurlSelector(CAND_LIST);

    // Time stepping loop.
    int ts = 1;
    do
    {
      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      Hermes::Mixins::Loggable::Static::info("\nRunge-Kutta time step (t = %g s, time_step = %g s, stages: %d).",
        current_time, time_step, bt.get_size());

      // Periodic global derefinements.
      if (ts > 1 && ts % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0)
      {
        Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
        REFINEMENT_COUNT = 0;

        E_space->unrefine_all_mesh_elements(true);
        H_space->unrefine_all_mesh_elements(true);
        P_space->unrefine_all_mesh_elements(true);

        E_space->adjust_element_order(-1, P_INIT);
        H_space->adjust_element_order(-1, P_INIT);
        P_space->adjust_element_order(-1, P_INIT);

        E_space->assign_dofs();
        H_space->assign_dofs();
        P_space->assign_dofs();
      }

      // Adaptivity loop:
      int as = 1;
      bool done = false;
      do
      {
        Hermes::Mixins::Loggable::Static::info("Adaptivity step %d:", as);

        // Construct globally refined reference mesh and setup reference space.
        int order_increase = 1;

        Mesh::ReferenceMeshCreator refMeshCreatorE(E_mesh);
        Mesh::ReferenceMeshCreator refMeshCreatorH(H_mesh);
        Mesh::ReferenceMeshCreator refMeshCreatorP(P_mesh);
        MeshSharedPtr ref_mesh_E = refMeshCreatorE.create_ref_mesh();
        MeshSharedPtr ref_mesh_H = refMeshCreatorH.create_ref_mesh();
        MeshSharedPtr ref_mesh_P = refMeshCreatorP.create_ref_mesh();

        Space<double>::ReferenceSpaceCreator refSpaceCreatorE(E_space, ref_mesh_E, order_increase);
        SpaceSharedPtr<double> ref_space_E = refSpaceCreatorE.create_ref_space();
        Space<double>::ReferenceSpaceCreator refSpaceCreatorH(H_space, ref_mesh_H, order_increase);
        SpaceSharedPtr<double> ref_space_H = refSpaceCreatorH.create_ref_space();
        Space<double>::ReferenceSpaceCreator refSpaceCreatorP(P_space, ref_mesh_P, order_increase);
        SpaceSharedPtr<double> ref_space_P = refSpaceCreatorP.create_ref_space();
        std::vector<SpaceSharedPtr<double> > ref_spaces({ ref_space_E, ref_space_H, ref_space_P });

        int ndof = Space<double>::get_num_dofs(ref_spaces);
        Hermes::Mixins::Loggable::Static::info("ndof = %d.", ndof);

        try
        {
          runge_kutta.set_spaces(ref_spaces);
          runge_kutta.set_time(current_time);
          runge_kutta.set_time_step(time_step);
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

        sprintf(title, "H, t = %g", current_time + time_step);
        H_view.set_title(title);
        H_view.show(H_time_new);

        sprintf(title, "P1, t = %g", current_time + time_step);
        P1_view.set_title(title);
        P1_view.show(P_time_new, H2D_FN_VAL_0);
        sprintf(title, "P2, t = %g", current_time + time_step);
        P2_view.set_title(title);
        P2_view.show(P_time_new, H2D_FN_VAL_1);

        // Project the fine mesh solution onto the coarse mesh.
        Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
        OGProjection<double>::project_global({ E_space, H_space, P_space }, { E_time_new, H_time_new, P_time_new }, { E_time_new_coarse, H_time_new_coarse, P_time_new_coarse });

        // Calculate element errors and total error estimate.
        Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
        adaptivity.set_spaces({ E_space, H_space, P_space });
        errorCalculator.calculate_errors({ E_time_new_coarse, H_time_new_coarse, P_time_new_coarse }, { E_time_new, H_time_new, P_time_new });

        double err_est_rel_total = errorCalculator.get_total_error_squared() * 100.;

        // Report results.
        Hermes::Mixins::Loggable::Static::info("Error estimate: %g%%", err_est_rel_total);

        // If err_est too large, adapt the mesh.
        if (err_est_rel_total < ERR_STOP)
        {
          Hermes::Mixins::Loggable::Static::info("Error estimate under the specified threshold -> moving to next time step.");
          done = true;
        }
        else
        {
          Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
          REFINEMENT_COUNT++;
          done = adaptivity.adapt({ &HcurlSelector, &H1selector, &HcurlSelector });

          if (!done)
            as++;
        }
      } while (!done);

      //View::wait();
      E_time_prev->copy(E_time_new);
      H_time_prev->copy(H_time_new);
      P_time_prev->copy(P_time_new);

      // Update time.
      current_time += time_step;
      ts++;
    } while (current_time < T_FINAL);

    // Wait for the view to be closed.
    View::wait();
  }
  catch (std::exception& e)
  {
    std::cout << e.what();
  }

  return 0;
}