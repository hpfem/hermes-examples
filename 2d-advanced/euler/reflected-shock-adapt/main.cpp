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
// Domain: A rectangular channel, see channel.mesh->
//
// BC: Solid walls, inlet / outlet. In this case there are two inlets (the left, and the upper wall).
//
// IC: Constant state.
//
// The following parameters can be changed:

// Frequently changed parameters.
// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = false;
// Maximum polynomial degree used. -1 for unlimited.
const int MAX_P_ORDER = 1;
// Time interval length.
const double T_END = 1.0;
// Shock capturing.
bool SHOCK_CAPTURING = true;
// Stopping criterion for adaptivity.
double ERR_STOP = 0.95;  
// Predefined list of element refinement candidates.
CandList CAND_LIST = H2D_HP_ANISO;                

// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;

// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;

// Initial polynomial degree.      
const int P_INIT = 0;                                             

// Number of initial uniform mesh refinements.  
const int INIT_REF_NUM = 4;                                            

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
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("channel.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) 
    mesh->refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT));
  SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT));
  SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT));
  SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT));
  int ndof = Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double>  >(space_rho, space_rho_v_x, space_rho_v_y, space_e));
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  MeshFunctionSharedPtr<double> prev_rho(new ConstantSolution<double>(mesh, RHO_INIT));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double> (mesh, RHO_INIT * V1_INIT));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double> (mesh, RHO_INIT * V2_INIT));
  MeshFunctionSharedPtr<double> prev_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_INIT, RHO_INIT * V1_INIT, RHO_INIT * V2_INIT, PRESSURE_INIT, KAPPA)));

  MeshFunctionSharedPtr<double> sln_rho(new ConstantSolution<double>(mesh, RHO_INIT));
  MeshFunctionSharedPtr<double> sln_rho_v_x(new ConstantSolution<double> (mesh, RHO_INIT * V1_INIT));
  MeshFunctionSharedPtr<double> sln_rho_v_y(new ConstantSolution<double> (mesh, RHO_INIT * V2_INIT));
  MeshFunctionSharedPtr<double> sln_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_INIT, RHO_INIT * V1_INIT, RHO_INIT * V2_INIT, PRESSURE_INIT, KAPPA)));

  MeshFunctionSharedPtr<double> rsln_rho(new ConstantSolution<double>(mesh, RHO_INIT));
  MeshFunctionSharedPtr<double> rsln_rho_v_x(new ConstantSolution<double> (mesh, RHO_INIT * V1_INIT));
  MeshFunctionSharedPtr<double> rsln_rho_v_y(new ConstantSolution<double> (mesh, RHO_INIT * V2_INIT));
  MeshFunctionSharedPtr<double> rsln_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_INIT, RHO_INIT * V1_INIT, RHO_INIT * V2_INIT, PRESSURE_INIT, KAPPA)));

  // Initialize weak formulation.
  Hermes::vector<std::string> solid_wall_markers;
  solid_wall_markers.push_back(BDY_SOLID_WALL);
  
  Hermes::vector<std::string> inlet_markers(BDY_INLET_LEFT, BDY_INLET_TOP);
  Hermes::vector<double> rho_ext(RHO_LEFT, RHO_TOP);
  Hermes::vector<double> v1_ext(V1_LEFT, V1_TOP);
  Hermes::vector<double> v2_ext(V2_LEFT, V2_TOP);
  Hermes::vector<double> pressure_ext(PRESSURE_LEFT, PRESSURE_TOP);

  Hermes::vector<std::string> outlet_markers;
  outlet_markers.push_back(BDY_OUTLET);

  EulerEquationsWeakFormSemiImplicit wf(KAPPA, rho_ext, v1_ext, v2_ext, pressure_ext, solid_wall_markers, inlet_markers, outlet_markers, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, (P_INIT == 0));

  // Filters for visualization of Mach number, pressure and entropy.
  MeshFunctionSharedPtr<double> Mach_number(new MachNumberFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), KAPPA));
  MeshFunctionSharedPtr<double> pressure(new PressureFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), KAPPA));
  MeshFunctionSharedPtr<double> entropy(new EntropyFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), KAPPA, RHO_INIT, P_INIT));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  Orderizer orderizer;
  OrderView space_view("Space", new WinGeom(0, 0, 1500, 700));

  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, MAX_P_ORDER);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  DiscreteProblemLinear<double> dp(&wf, Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
      space_rho_v_y, space_e));
  LinearSolver<double> solver(&dp);

  // Time stepping loop.
  double t = 0.0;
  int iteration = 0;
  for(; t < T_END && iteration < 25; t += time_step)
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

        space_rho->unrefine_all_mesh_elements(true);

        space_rho->adjust_element_order(-1, P_INIT);
        space_rho_v_x->adjust_element_order(-1, P_INIT);
        space_rho_v_y->adjust_element_order(-1, P_INIT);
        space_e->adjust_element_order(-1, P_INIT);

        Space<double>::assign_dofs(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, space_rho_v_y, space_e));
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
      int order_increase = (CAND_LIST == H2D_HP_ANISO ? 1 : 0);

      Mesh::ReferenceMeshCreator refMeshCreatorFlow(mesh);
      MeshSharedPtr ref_mesh = refMeshCreatorFlow.create_ref_mesh();

      Space<double>::ReferenceSpaceCreator refSpaceCreatorRho(space_rho, ref_mesh, order_increase);
      SpaceSharedPtr<double> ref_space_rho = refSpaceCreatorRho.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVx(space_rho_v_x, ref_mesh, order_increase);
      SpaceSharedPtr<double> ref_space_rho_v_x = refSpaceCreatorRhoVx.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVy(space_rho_v_y, ref_mesh, order_increase);
      SpaceSharedPtr<double> ref_space_rho_v_y = refSpaceCreatorRhoVy.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorE(space_e, ref_mesh, order_increase);
      SpaceSharedPtr<double> ref_space_e = refSpaceCreatorE.create_ref_space();

      Hermes::vector<SpaceSharedPtr<double> > ref_spaces(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e);
      Hermes::vector<SpaceSharedPtr<double>  > ref_spaces_const(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e);

      if(ndofs_prev != 0)
        if(Space<double>::get_num_dofs(ref_spaces_const) == ndofs_prev)
          selector.set_error_weights(2.0 * selector.get_error_weight_h(), 1.0, 1.0);
        else
          selector.set_error_weights(1.0, 1.0, 1.0);

      ndofs_prev = Space<double>::get_num_dofs(ref_spaces_const);

      // Project the previous time level solution onto the new fine mesh->
      Hermes::Mixins::Loggable::Static::info("Projecting the previous time level solution onto the new fine mesh->");
        OGProjection<double> ogProjection;
        ogProjection.project_global(ref_spaces_const, Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 
            Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), Hermes::vector<Hermes::Hermes2D::ProjNormType>());

      // Report NDOFs.
      Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d.", 
        Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
        space_rho_v_y, space_e)), Space<double>::get_num_dofs(ref_spaces_const));

      // Assemble the reference problem.
      Hermes::Mixins::Loggable::Static::info("Solving on reference mesh->");
      
      solver.set_spaces(ref_spaces_const);
      wf.set_current_time_step(time_step);
    
      solver.solve();

      if(!SHOCK_CAPTURING)
        Solution<double>::vector_to_solutions(solver.get_sln_vector(), ref_spaces_const, 
        Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e));
      else
      {      
        FluxLimiter flux_limiter(FluxLimiter::Kuzmin, solver.get_sln_vector(), ref_spaces_const, true);

        flux_limiter.limit_second_orders_according_to_detector(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
          space_rho_v_y, space_e));

        flux_limiter.limit_according_to_detector(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
          space_rho_v_y, space_e));

        flux_limiter.get_limited_solutions(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e));
      }

      // Project the fine mesh solution onto the coarse mesh->
      Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh->");
      ogProjection.project_global(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
        space_rho_v_y, space_e), Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), 
        Hermes::vector<MeshFunctionSharedPtr<double> >(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e), 
        Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
      Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
        space_rho_v_y, space_e), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<MeshFunctionSharedPtr<double> >(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e),
        Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e)) * 100;

      CFL.calculate_semi_implicit(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), (ref_spaces_const)[0]->get_mesh(), time_step);

      // Report results.
      Hermes::Mixins::Loggable::Static::info("err_est_rel: %g%%", err_est_rel_total);

      // If err_est too large, adapt the mesh->
      if (Space<double>::get_num_dofs(ref_spaces) > NDOF_STOP || err_est_rel_total < THRESHOLD)
      {
        done = true;
      }
      else
      {
        Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh->");
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
          Mach_number->reinit();
          //pressure->reinit();
          //entropy->reinit();
          //pressure_view.show(&pressure);
          Mach_number_view.show(Mach_number);
          space_view.show(ref_space_rho);

        }
        // Output solution in VTK format.
        if(VTK_VISUALIZATION) 
        {
          pressure->reinit();
          Mach_number->reinit();
          Linearizer lin;
          char filename[40];
          //sprintf(filename, "Pressure-%i.vtk", iteration - 1);
          //lin.save_solution_vtk(&pressure, filename, "Pressure", false);
          sprintf(filename, "Mach number-%i.vtk", iteration - 1);
          lin.save_solution_vtk(Mach_number, filename, "MachNumber", false);
          sprintf(filename, "Mesh-%i.vtk", iteration - 1);
          orderizer.save_mesh_vtk(ref_space_rho, filename);
          sprintf(filename, "Space-%i.vtk", iteration - 1);
          orderizer.save_orders_vtk(ref_space_rho, filename);
        }
      }

      // Clean up.
      delete adaptivity;
    }
    while (done == false);


    // Copy the solutions into the previous time level ones.
    prev_rho->copy(rsln_rho);
    prev_rho_v_x->copy(rsln_rho_v_x);
    prev_rho_v_y->copy(rsln_rho_v_y);
    prev_e->copy(rsln_e);
  }

  if(HERMES_VISUALIZATION) 
  {
    pressure_view.close();
    entropy_production_view.close();
    Mach_number_view.close();
  }

  return 0;
}
