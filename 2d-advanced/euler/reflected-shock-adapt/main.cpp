#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// This example solves the compressible Euler equations using Discontinuous Galerkin method of higher order with adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: A rectangular channel, see channel.mesh.
//
// BC: Solid walls, inlet / outlet. In this case there are two inlets (the left, and the upper wall).
//
// IC: Constant state.
//
// The following parameters can be changed:

// Frequently changed parameters.
// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;
// Maximum polynomial degree used. -1 for unlimited.
const int MAX_P_ORDER = 0;
// Time interval length.
const double T_END = 0.05;
// Shock capturing.
bool SHOCK_CAPTURING = true;
// Stopping criterion for adaptivity.
double ERR_STOP = 0.95;  
// Predefined list of element refinement candidates.
CandList CAND_LIST = H2D_HP_ANISO;                

// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = true;

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;

// For saving/loading of solution.
bool REUSE_SOLUTION = false;

// Initial polynomial degree.      
const int P_INIT = 0;                                             

// Number of initial uniform mesh refinements.  
const int INIT_REF_NUM = 3;                                            

// CFL value.
double CFL_NUMBER = 0.1;                         

// Initial time step.
double time_step = 1E-4;

// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
const MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Adaptivity.
// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 5;

// Number of mesh refinements between two unrefinements.
// The mesh is not unrefined unless there has been a refinement since
// last unrefinement.
int REFINEMENT_COUNT = 0;

// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.5;                     

// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
const int STRATEGY = 1;                           

// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   

// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0. 
const double CONV_EXP = 1;                        

// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 100000;

// Equation parameters.
const double KAPPA = 1.4;                         
const double RHO_LEFT = 1.0;
const double RHO_TOP = 1.7;

const double V1_LEFT = 2.9;
const double V1_TOP = 2.619334;

const double V2_LEFT = 0.0;
const double V2_TOP = -0.5063;

const double PRESSURE_LEFT = 0.714286;
const double PRESSURE_TOP = 1.52819;

// Initial values
const double RHO_INIT = 1.0;
const double V1_INIT = 2.9;
const double V2_INIT = 0.0;
const double PRESSURE_INIT = 0.714286;

// Boundary markers.
const std::string BDY_SOLID_WALL = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_INLET_TOP = "3";
const std::string BDY_INLET_LEFT = "4";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  Hermes::Hermes2D::Hermes2DApi.setParamValue(Hermes::Hermes2D::numThreads, 8);

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("channel.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) 
    mesh.refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  L2Space<double> space_rho(&mesh, P_INIT);
  L2Space<double> space_rho_v_x(&mesh, P_INIT);
  L2Space<double> space_rho_v_y(&mesh, P_INIT);
  L2Space<double> space_e(&mesh, P_INIT);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  ConstantSolution<double> sln_rho(&mesh, RHO_INIT);
  ConstantSolution<double> sln_rho_v_x(&mesh, RHO_INIT * V1_INIT);
  ConstantSolution<double> sln_rho_v_y(&mesh, RHO_INIT * V2_INIT);
  ConstantSolution<double> sln_e(&mesh, QuantityCalculator::calc_energy(RHO_INIT, RHO_INIT * V1_INIT, RHO_INIT * V2_INIT, PRESSURE_INIT, KAPPA));

  ConstantSolution<double> prev_rho(&mesh, RHO_INIT);
  ConstantSolution<double> prev_rho_v_x(&mesh, RHO_INIT * V1_INIT);
  ConstantSolution<double> prev_rho_v_y(&mesh, RHO_INIT * V2_INIT);
  ConstantSolution<double> prev_e(&mesh, QuantityCalculator::calc_energy(RHO_INIT, RHO_INIT * V1_INIT, RHO_INIT * V2_INIT, PRESSURE_INIT, KAPPA));

  ConstantSolution<double> rsln_rho(&mesh, RHO_INIT);
  ConstantSolution<double> rsln_rho_v_x(&mesh, RHO_INIT * V1_INIT);
  ConstantSolution<double> rsln_rho_v_y(&mesh, RHO_INIT * V2_INIT);
  ConstantSolution<double> rsln_e(&mesh, QuantityCalculator::calc_energy(RHO_INIT, RHO_INIT * V1_INIT, RHO_INIT * V2_INIT, PRESSURE_INIT, KAPPA));

  // Initialize weak formulation.
  EulerEquationsWeakFormSemiImplicitTwoInflows wf(KAPPA, RHO_LEFT, V1_LEFT, V2_LEFT, PRESSURE_LEFT, RHO_TOP, V1_TOP, V2_TOP, PRESSURE_TOP, BDY_SOLID_WALL, BDY_INLET_LEFT, BDY_INLET_TOP, BDY_OUTLET,
    &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, (P_INIT == 0 && CAND_LIST == H2D_H_ANISO));

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), KAPPA, RHO_INIT, P_INIT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  OrderView space_view("Space", new WinGeom(0, 0, 1500, 700));

  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, MAX_P_ORDER);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Look for a saved solution on the disk.
  CalculationContinuity<double> continuity(CalculationContinuity<double>::onlyTime);
  int iteration = 0; double t = 0;
  bool loaded_now = false;

  if(REUSE_SOLUTION && continuity.have_record_available())
  {
    continuity.get_last_record()->load_mesh(&mesh);
    continuity.get_last_record()->load_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::vector<SpaceType>(HERMES_L2_SPACE, HERMES_L2_SPACE, HERMES_L2_SPACE, HERMES_L2_SPACE), Hermes::vector<Mesh *>(&mesh, &mesh, 
      &mesh, &mesh));
    continuity.get_last_record()->load_time_step_length(time_step);
    t = continuity.get_last_record()->get_time() + time_step;
    iteration = (continuity.get_num()) * EVERY_NTH_STEP + 1;
    loaded_now = true;
  }

  DiscreteProblemLinear<double> dp(&wf, Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e));
  LinearSolver<double> solver(&dp);

  if(P_INIT == 0 && CAND_LIST == H2D_H_ANISO)
    dp.set_fvm();

  Hermes::vector<Space<double>*> spaces_to_delete;
      
  // Time stepping loop.
  for(; t < T_END; t += time_step)
  {
    CFL.set_number(CFL_NUMBER + (t/4.0) * 1.0);
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0)
    {
      if(REFINEMENT_COUNT > 0) 
      {
        Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
        REFINEMENT_COUNT = 0;

        space_rho.unrefine_all_mesh_elements(true);

        space_rho.adjust_element_order(-1, P_INIT);
        space_rho_v_x.copy_orders(&space_rho);
        space_rho_v_y.copy_orders(&space_rho);
        space_e.copy_orders(&space_rho);
      }
    }

    // Adaptivity loop:
    int as = 1; 
    int ndofs_prev = 0;
    bool done = false;
    do
    {
      Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      int order_increase = 0;

      Hermes::vector<Space<double> *>* ref_spaces = Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), order_increase);

      Hermes::vector<const Space<double> *> ref_spaces_const((*ref_spaces)[0], (*ref_spaces)[1], 
        (*ref_spaces)[2], (*ref_spaces)[3]);

      if(ndofs_prev != 0)
        if(Space<double>::get_num_dofs(ref_spaces_const) == ndofs_prev)
          selector.set_error_weights(2.0 * selector.get_error_weight_h(), 1.0, 1.0);
        else
          selector.set_error_weights(1.0, 1.0, 1.0);

      ndofs_prev = Space<double>::get_num_dofs(ref_spaces_const);

      // Project the previous time level solution onto the new fine mesh.
      Hermes::Mixins::Loggable::Static::info("Projecting the previous time level solution onto the new fine mesh.");
      if(loaded_now)
      {
        loaded_now = false;

        continuity.get_last_record()->load_solutions(Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
            Hermes::vector<Space<double> *>((*ref_spaces)[0], (*ref_spaces)[1], (*ref_spaces)[2], (*ref_spaces)[3]));
      }
      else
      {
        OGProjection<double> ogProjection;
        ogProjection.project_global(ref_spaces_const, Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
            Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), Hermes::vector<Hermes::Hermes2D::ProjNormType>());
        if(iteration > 1)
        {
          for(int i = 0; i < spaces_to_delete.size(); i++)
          {
            delete spaces_to_delete[i]->get_mesh();
            delete spaces_to_delete[i];
          }
        }
      }

      spaces_to_delete.clear();
      for(int i = 0; i < ref_spaces->size(); i++)
        spaces_to_delete.push_back((*ref_spaces)[i]);

      // Report NDOFs.
      Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d.", 
        Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e)), Space<double>::get_num_dofs(ref_spaces_const));

      // Assemble the reference problem.
      Hermes::Mixins::Loggable::Static::info("Solving on reference mesh.");
      
      solver.set_spaces(ref_spaces_const);
      wf.set_time_step(time_step);
    
      solver.solve();
      if(!SHOCK_CAPTURING)
        Solution<double>::vector_to_solutions(solver.get_sln_vector(), ref_spaces_const, 
        Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
      else
      {      
        FluxLimiter flux_limiter(FluxLimiter::Kuzmin, solver.get_sln_vector(), ref_spaces_const, true);

        flux_limiter.limit_second_orders_according_to_detector(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e));

        flux_limiter.limit_according_to_detector(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e));

        flux_limiter.get_limited_solutions(Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
      }

      // Project the fine mesh solution onto the coarse mesh.
      Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
      OGProjection<double> ogProjection;
      ogProjection.project_global(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), 
        Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), 
        Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
      Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e),
        Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e)) * 100;

      CFL.calculate_semi_implicit(Hermes::vector<Solution<double> *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), (ref_spaces_const)[0]->get_mesh(), time_step);

      // Report results.
      Hermes::Mixins::Loggable::Static::info("err_est_rel: %g%%", err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (Space<double>::get_num_dofs(*ref_spaces) > NDOF_STOP || err_est_rel_total < THRESHOLD)
      {
        done = true;
      }
      else
      {
        Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
        REFINEMENT_COUNT++;
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector, &selector, &selector), 
          THRESHOLD, STRATEGY, MESH_REGULARITY);

        if(!done)
          as++;
      }

      // Visualization and saving on disk.
      //if(done)
      {
        // Hermes visualization.
        if(HERMES_VISUALIZATION) 
        {
          Mach_number.reinit();
          //pressure.reinit();
          //entropy.reinit();
          //pressure_view.show(&pressure);
          Mach_number_view.show(&Mach_number);
          //pressure_view.save_numbered_screenshot("Pressure-%u.bmp", iteration - 1, true);
          Mach_number_view.save_numbered_screenshot("Mach-%u.bmp", iteration - 1, true);
          space_view.show((*ref_spaces)[0]);
          space_view.save_numbered_screenshot("Space-%u.bmp", iteration - 1, true);
        }
        // Output solution in VTK format.
        if(VTK_VISUALIZATION) 
        {
          pressure.reinit();
          Mach_number.reinit();
          Linearizer lin;
          char filename[40];
          //sprintf(filename, "Pressure-%i.vtk", iteration - 1);
          //lin.save_solution_vtk(&pressure, filename, "Pressure", false);
          sprintf(filename, "Mach number-%i.vtk", iteration - 1);
          lin.save_solution_vtk(&Mach_number, filename, "MachNumber", false);
        }
      }

      // Clean up.
      delete adaptivity;
    }
    while (done == false);


    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&rsln_rho);
    prev_rho_v_x.copy(&rsln_rho_v_x);
    prev_rho_v_y.copy(&rsln_rho_v_y);
    prev_e.copy(&rsln_e);
  }

  if(HERMES_VISUALIZATION) 
  {
    pressure_view.close();
    entropy_production_view.close();
    Mach_number_view.close();
  }

  return 0;
}
