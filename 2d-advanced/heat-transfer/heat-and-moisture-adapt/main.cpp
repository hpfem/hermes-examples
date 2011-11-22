#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
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
const double THRESHOLD = 0.3;                     
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
const int STRATEGY = 0;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;          
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0. 
const double CONV_EXP = 1.0;                      
// Stopping criterion for adaptivity.
const double ERR_STOP = 2.0;                      
// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 100000;                     
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Newton's method
// Stopping criterion for Newton on fine mesh.
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

// Problem parameters.
const double c_TT = 2.18e+6;
const double d_TT = 2.1;
const double d_Tw = 2.37e-2 / W_SCALING_FACTOR;
const double k_TT = 25;
const double c_ww = 24.9;
const double d_wT = 1.78e-10 * W_SCALING_FACTOR;
const double d_ww = 3.02e-8;
const double k_ww = 1.84e-7;

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
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh basemesh, T_mesh, w_mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &basemesh);

  // Create temperature and moisture meshes.
  // This also initializes the multimesh hp-FEM.
  T_mesh.copy(&basemesh);
  w_mesh.copy(&basemesh);

  // Initialize boundary conditions.
  EssentialBCNonConst temp_reactor("bdy_react", REACTOR_START_TIME, T_INITIAL, T_REACTOR_MAX);
  EssentialBCs<double> bcs_T(&temp_reactor);

  // Create H1 spaces with default shapesets.
  H1Space<double> T_space(&T_mesh, &bcs_T, P_INIT);
  H1Space<double> w_space(MULTI ? &w_mesh : &T_mesh, P_INIT);

  // Define constant initial conditions.
  info("Setting initial conditions.");
  ConstantSolution<double> T_time_prev(&T_mesh, T_INITIAL);
  ConstantSolution<double> w_time_prev(&w_mesh, W_INITIAL);
  Solution<double> T_time_new(&T_mesh);
  Solution<double> w_time_new(&w_mesh);
  
  // Solutions.
  Solution<double> T_coarse, w_coarse;

  // Initialize the weak formulation.
  CustomWeakFormHeatMoistureRK wf(c_TT, c_ww, d_TT, d_Tw, d_wT, d_ww, 
				  k_TT, k_ww, T_EXTERIOR, W_EXTERIOR, "bdy_ext");

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

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
  T_sln_view.show(&T_time_prev);
  w_sln_view.show(&w_time_prev);
  T_order_view.show(&T_space);
  w_order_view.show(&w_space);

  // Time stepping loop:
  int ts = 1;
  while (current_time < SIMULATION_TIME)
  {
    info("Simulation time = %g s (%d h, %d d, %d y)",
        current_time, (int) current_time / 3600,
        (int) current_time / (3600*24), (int) current_time / (3600*24*364));

    // Update time-dependent essential BCs.
    if (current_time <= REACTOR_START_TIME) {
      info("Updating time-dependent essential BC.");
      Space<double>::update_essential_bc_values(Hermes::vector<Space<double>*>(&T_space, &w_space), current_time);
    }

    // Uniform mesh derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
        case 1: T_mesh.copy(&basemesh);
                w_mesh.copy(&basemesh);
                T_space.set_uniform_order(P_INIT);
                w_space.set_uniform_order(P_INIT);
                break;
        case 2: T_mesh.unrefine_all_elements();
                if(MULTI)
                  w_mesh.unrefine_all_elements();
                T_space.set_uniform_order(P_INIT);
                w_space.set_uniform_order(P_INIT);
                break;
        case 3: T_mesh.unrefine_all_elements();
                if(MULTI)
                  w_mesh.unrefine_all_elements();
                T_space.adjust_element_order(-1, -1, P_INIT, P_INIT);
                w_space.adjust_element_order(-1, -1, P_INIT, P_INIT);
                break;
        default: error("Wrong global derefinement method.");
      }
    }

    // Spatial adaptivity loop. Note: T_time_prev and w_time_prev must not be changed during 
    // spatial adaptivity.
    bool done = false; int as = 1; 
    do
    {
      info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh and setup reference space.
      Hermes::vector<Space<double> *>* ref_spaces 
          = Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&T_space, &w_space));
      Hermes::vector<const Space<double> *> ref_spaces_const((*ref_spaces)[0], (*ref_spaces)[1]);

      // Initialize discrete problem on reference meshes.
      DiscreteProblem<double> dp(&wf, ref_spaces_const);

      // Initialize Runge-Kutta time stepping.
      RungeKutta<double> runge_kutta(&wf, ref_spaces_const, &bt, matrix_solver);

      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      info("Runge-Kutta time step (t = %g s, tau = %g s, stages: %d).",
           current_time, time_step, bt.get_size());
      bool freeze_jacobian = true;
      bool block_diagonal_jacobian = true;
      bool verbose = true;
      
      try
      {
        runge_kutta.rk_time_step_newton(current_time, time_step, 
            Hermes::vector<Solution<double> *>(&T_time_prev, &w_time_prev), 
            Hermes::vector<Solution<double> *>(&T_time_new, &w_time_new), 
            freeze_jacobian, block_diagonal_jacobian, verbose, NEWTON_TOL, NEWTON_MAX_ITER);
      }
      catch(Exceptions::Exception& e)
      {
        e.printMsg();
        error("Runge-Kutta time step failed");
      }

      // Project the fine mesh solution onto the coarse meshes.
      info("Projecting fine mesh solutions on coarse meshes for error estimation.");
      OGProjection<double>::project_global(Hermes::vector<Space<double> *>(&T_space, &w_space), 
          Hermes::vector<Solution<double> *>(&T_time_new, &w_time_new), 
	  Hermes::vector<Solution<double> *>(&T_coarse, &w_coarse),
          matrix_solver); 

      // Initialize an instance of the Adapt class and register custom error forms.
      Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<Space<double> *>(&T_space, &w_space));
      /* ADAPT IN ENERGY NORM 
      CustomErrorForm cef_0_0(d_TT, c_TT);
      CustomErrorForm cef_0_1(d_Tw, c_TT);
      CustomErrorForm cef_1_0(d_wT, c_ww);
      CustomErrorForm cef_1_1(d_ww, c_ww);
      adaptivity->set_error_form(0, 0, &cef_0_0);
      adaptivity->set_error_form(0, 1, &cef_0_1);
      adaptivity->set_error_form(1, 0, &cef_1_0);
      adaptivity->set_error_form(1, 1, &cef_1_1);
      */

      // Calculate element errors and total error estimate.
      info("Calculating error estimate."); 
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution<double> *>(&T_coarse, &w_coarse), 
                                 Hermes::vector<Solution<double> *>(&T_time_new, &w_time_new)) * 100;

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
	   Space<double>::get_num_dofs(Hermes::vector<Space<double> *>(&T_space, &w_space)), 
           Space<double>::get_num_dofs(ref_spaces_const), err_est_rel_total);

      // If err_est too large, adapt the meshes.
      if (err_est_rel_total < ERR_STOP)
        done = true;
      else 
      {
        info("Adapting the coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector), THRESHOLD, STRATEGY, MESH_REGULARITY);

        if (Space<double>::get_num_dofs(Hermes::vector<Space<double> *>(&T_space, &w_space)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }
 
      // Clean up.
      delete adaptivity;
      if(!done)
      {
        delete ref_spaces;
        delete T_time_new.get_mesh();
        delete w_time_new.get_mesh();
      }

      // Increase counter.
      as++;
    }
    while (done == false);

    // Update time.
    current_time += time_step;

    // Show new coarse meshes and solutions.
    char title[100];
    sprintf(title, "Temperature, t = %g days", current_time/3600./24);
    T_sln_view.set_title(title);
    T_sln_view.show(&T_coarse);
    sprintf(title, "Moisture (scaled), t = %g days", current_time/3600./24);
    w_sln_view.set_title(title);
    w_sln_view.show(&w_coarse);
    T_order_view.show(&T_space);
    w_order_view.show(&w_space);

    // Save fine mesh solutions for the next time step.
    T_time_prev.copy(&T_time_new);
    w_time_prev.copy(&w_time_new);

    ts++;
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
