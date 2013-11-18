#include "definitions.h"

using namespace RefinementSelectors;

// This example solves adaptively a time-dependent coupled problem of heat and moisture 
// transfer in massive concrete walls of a nuclear reactor vessel (simplified axi-symmetric 
// geometry). 
//
// PDE: Lengthy. See the paper P. Solin, L. Dubcova, J. Kruis: Adaptive hp-FEM with Dynamical 
// Meshes for Transient Heat and Moisture Transfer Problems, J. Comput. Appl. Math. 233 (2010) 3103-3112.
//
// The following parameters can be changed:

// Scaling factor for moisture. Since temperature is in hundreds of Kelvins and moisture between
// (0, 1), adaptivity works better when moisture is scaled. 
const double W_SCALING_FACTOR = 100.;
// Initial polynomial degrees.
const int P_INIT = 2;                             
// MULTI = true  ... use multi-mesh,
// MULTI = false ... use single-mesh.
// Note: In the single mesh option, the meshes are
// forced to be geometrically the same but the
// polynomial degrees can still vary.
const bool MULTI = true;                          
// Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_FREQ = 1;                         
// 1... mesh reset to basemesh and poly degrees to P_INIT.   
// 2... one ref. layer shaved off, poly degrees reset to P_INIT.
// 3... one ref. layer shaved off, poly degrees decreased by one. 
// and just one polynomial degree subtracted.
const int UNREF_METHOD = 3;                       
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.9;                     
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e0;

// Problem parameters.
const double c_TT = 2.18e+6;
const double d_TT = 2.1;
const double d_Tw = 2.37e-2 / W_SCALING_FACTOR;
const double k_TT = 25;
const double c_ww = 24.9;
const double d_wT = 1.78e-10 * W_SCALING_FACTOR;
const double d_ww = 3.02e-8;
const double k_ww = 1.84e-7;

CustomErrorForm cef_0_0(0, 0, d_TT, c_TT);
CustomErrorForm cef_0_1(0, 1, d_Tw, c_TT);
CustomErrorForm cef_1_0(1, 0, d_wT, c_ww);
CustomErrorForm cef_1_1(1, 1, d_ww, c_ww);

// Error calculation & adaptivity.
class MyErrorCalculator : public Hermes::Hermes2D::ErrorCalculator<double>
{
public:
  MyErrorCalculator() : Hermes::Hermes2D::ErrorCalculator<double>(RelativeErrorToGlobalNorm)
  {
    this->add_error_form(&cef_0_0);
    this->add_error_form(&cef_0_1);
    this->add_error_form(&cef_1_0);
    this->add_error_form(&cef_1_1);
  }
}
errorCalculator;
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);

// Newton's method
// Stopping criterion for Newton on fine mesh->
const double NEWTON_TOL = 1e-5;                   
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
ButcherTableType butcher_table_type = Implicit_RK_1;

// Time step and simulation time.
// Time step: 10 days
const double time_step = 10.*24*60*60;                  
// Physical time [seconds].
const double SIMULATION_TIME = 10000*time_step + 0.001;  

// Initial and boundary conditions.
// (Kelvins)
const double T_INITIAL = 293.0;           
// (dimensionless)
const double W_INITIAL = 0.9 * W_SCALING_FACTOR;            
// (Kelvins)
const double T_EXTERIOR = 293.0;          
// (dimensionless)
const double W_EXTERIOR = 0.55 * W_SCALING_FACTOR;          
// (Kelvins)
const double T_REACTOR_MAX = 550.0;       
// How long does the reactor
// need to warm up linearly from T_INITIAL
// to T_REACTOR_MAX [seconds].
const double REACTOR_START_TIME = 3600*24;   

// Physical time in seconds.
double current_time = 0.0;

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) Hermes::Mixins::Loggable::Static::info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  MeshSharedPtr basemesh(new Mesh), T_mesh(new Mesh), w_mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", basemesh);

  // Create temperature and moisture meshes.
  // This also initializes the multimesh hp-FEM.
  T_mesh->copy(basemesh);
  w_mesh->copy(basemesh);

  // Initialize boundary conditions.
  EssentialBCNonConst temp_reactor("bdy_react", REACTOR_START_TIME, T_INITIAL, T_REACTOR_MAX);
  EssentialBCs<double> bcs_T(&temp_reactor);

  SpaceSharedPtr<double> T_space(new H1Space<double>(T_mesh, &bcs_T, P_INIT));
  SpaceSharedPtr<double> w_space(new H1Space<double>(MULTI ? w_mesh : T_mesh, P_INIT));
  Hermes::vector<SpaceSharedPtr<double> > spaces(T_space, w_space);
  adaptivity.set_spaces(spaces);

  // Define constant initial conditions.
  Hermes::Mixins::Loggable::Static::info("Setting initial conditions.");
  MeshFunctionSharedPtr<double> T_time_prev(new ConstantSolution<double>(T_mesh, T_INITIAL));
  MeshFunctionSharedPtr<double> w_time_prev(new ConstantSolution<double>(w_mesh, W_INITIAL));
  MeshFunctionSharedPtr<double> T_time_new(new Solution<double>(T_mesh));
  MeshFunctionSharedPtr<double> w_time_new(new Solution<double>(w_mesh));

  // Solutions.
  MeshFunctionSharedPtr<double> T_coarse(new Solution<double>), w_coarse(new Solution<double>);

  // Initialize the weak formulation.
  CustomWeakFormHeatMoistureRK wf(c_TT, c_ww, d_TT, d_Tw, d_wT, d_ww, 
    k_TT, k_ww, T_EXTERIOR, W_EXTERIOR, "bdy_ext");

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  // Geometry and position of visualization windows.
  WinGeom* T_sln_win_geom = new WinGeom(0, 0, 300, 450);
  WinGeom* w_sln_win_geom = new WinGeom(310, 0, 300, 450);
  WinGeom* T_mesh_win_geom = new WinGeom(620, 0, 280, 450);
  WinGeom* w_mesh_win_geom = new WinGeom(910, 0, 280, 450);

  // Initialize views.
  ScalarView T_sln_view("Temperature", T_sln_win_geom);
  ScalarView w_sln_view("Moisture (scaled)", w_sln_win_geom);
  OrderView T_order_view("Temperature mesh", T_mesh_win_geom);
  OrderView w_order_view("Moisture mesh", w_mesh_win_geom);

  // Show initial conditions.
  T_sln_view.show(T_time_prev);
  w_sln_view.show(w_time_prev);
  T_order_view.show(T_space);
  w_order_view.show(w_space);

  // Time stepping loop:
  int ts = 1;
  while (current_time < SIMULATION_TIME)
  {
    Hermes::Mixins::Loggable::Static::info("Simulation time = %g s (%d h, %d d, %d y)",
      current_time, (int) current_time / 3600,
      (int) current_time / (3600*24), (int) current_time / (3600*24*364));

    // Update time-dependent essential BCs.
    if (current_time <= REACTOR_START_TIME) {
      Hermes::Mixins::Loggable::Static::info("Updating time-dependent essential BC.");
      Space<double>::update_essential_bc_values(Hermes::vector<SpaceSharedPtr<double> >(T_space, w_space), current_time);
    }

    // Uniform mesh derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
      case 1: T_mesh->copy(basemesh);
        w_mesh->copy(basemesh);
        T_space->set_uniform_order(P_INIT);
        w_space->set_uniform_order(P_INIT);
        break;
      case 2: T_mesh->unrefine_all_elements();
        if(MULTI)
          w_mesh->unrefine_all_elements();
        T_space->set_uniform_order(P_INIT);
        w_space->set_uniform_order(P_INIT);
        break;
      case 3: T_mesh->unrefine_all_elements();
        if(MULTI)
          w_mesh->unrefine_all_elements();
        T_space->adjust_element_order(-1, -1, P_INIT, P_INIT);
        w_space->adjust_element_order(-1, -1, P_INIT, P_INIT);
        break;
      default: throw Hermes::Exceptions::Exception("Wrong global derefinement method.");
      }
      T_space->assign_dofs();
      w_space->assign_dofs();
      Space<double>::assign_dofs(spaces);
    }

    // Spatial adaptivity loop. Note: T_time_prev and w_time_prev must not be changed during 
    // spatial adaptivity.
    bool done = false; int as = 1; 
    do
    {
      Hermes::Mixins::Loggable::Static::info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh and setup reference space.
      Mesh::ReferenceMeshCreator refMeshCreatorU(T_mesh);
      MeshSharedPtr ref_T_mesh = refMeshCreatorU.create_ref_mesh();

      Space<double>::ReferenceSpaceCreator refSpaceCreatorT(T_space, ref_T_mesh);
      SpaceSharedPtr<double> ref_T_space = refSpaceCreatorT.create_ref_space();

      Mesh::ReferenceMeshCreator refMeshCreatorW(w_mesh);
      MeshSharedPtr ref_w_mesh = refMeshCreatorW.create_ref_mesh();

      Space<double>::ReferenceSpaceCreator refSpaceCreatorW(w_space, ref_w_mesh);
      SpaceSharedPtr<double> ref_w_space = refSpaceCreatorW.create_ref_space();

      Hermes::vector<SpaceSharedPtr<double> > ref_spaces(ref_T_space, ref_w_space);

      // Initialize discrete problem on reference meshes.
      DiscreteProblem<double> dp(&wf, ref_spaces);

      // Initialize Runge-Kutta time stepping.
      RungeKutta<double> runge_kutta(&wf, ref_spaces, &bt);

      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      Hermes::Mixins::Loggable::Static::info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).",
        current_time, time_step, bt.get_size());
      try
      {
        runge_kutta.set_time(current_time);
        runge_kutta.set_time_step(time_step);
        runge_kutta.set_max_allowed_iterations(NEWTON_MAX_ITER);
        runge_kutta.set_tolerance(NEWTON_TOL);
        runge_kutta.rk_time_step_newton(Hermes::vector<MeshFunctionSharedPtr<double> >(T_time_prev, w_time_prev), 
          Hermes::vector<MeshFunctionSharedPtr<double> >(T_time_new, w_time_new));
      }
      catch(Exceptions::Exception& e)
      {
        e.print_msg();
        throw Hermes::Exceptions::Exception("Runge-Kutta time step failed");
      }

      // Project the fine mesh solution onto the coarse meshes.
      Hermes::Mixins::Loggable::Static::info("Projecting fine mesh solutions on coarse meshes for error estimation.");
      OGProjection<double> ogProjection; ogProjection.project_global(Hermes::vector<SpaceSharedPtr<double> >(T_space, w_space), 
        Hermes::vector<MeshFunctionSharedPtr<double> >(T_time_new, w_time_new), 
        Hermes::vector<MeshFunctionSharedPtr<double> >(T_coarse, w_coarse)); 


      // Calculate element errors and total error estimate.
      Hermes::Mixins::Loggable::Static::info("Calculating error estimate."); 
      errorCalculator.calculate_errors(Hermes::vector<MeshFunctionSharedPtr<double> >(T_coarse, w_coarse), 
        Hermes::vector<MeshFunctionSharedPtr<double> >(T_time_new, w_time_new));
      double err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

      // Report results.
      Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
        Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(T_space, w_space)), 
        Space<double>::get_num_dofs(ref_spaces), err_est_rel_total);

      // If err_est too large, adapt the meshes.
      if (err_est_rel_total < ERR_STOP)
        done = true;
      else 
      {
        Hermes::Mixins::Loggable::Static::info("Adapting the coarse mesh.");
        done = adaptivity.adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector));

        // Increase the counter of performed adaptivity steps.
        as++;
      }
    }
    while (done == false);

    // Update time.
    current_time += time_step;

    // Show new coarse meshes and solutions.
    char title[100];
    sprintf(title, "Temperature, t = %g days", current_time/3600./24);
    T_sln_view.set_title(title);
    T_sln_view.show(T_coarse);
    sprintf(title, "Moisture (scaled), t = %g days", current_time/3600./24);
    w_sln_view.set_title(title);
    w_sln_view.show(w_coarse);
    T_order_view.show(T_space);
    w_order_view.show(w_space);

    // Save fine mesh solutions for the next time step.
    T_time_prev->copy(T_time_new);
    w_time_prev->copy(w_time_new);

    ts++;
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
