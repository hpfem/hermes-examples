#define HERMES_REPORT_INFO
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// This example solves the compressible Euler equations using Discontinuous Galerkin method of higher order with adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: Joukowski profile, see file domain-nurbs.xml.
//
// BC: Solid walls, inlet, outlet.
//
// IC: Constant state identical to inlet.
//
// The following parameters can be changed:

// Visualization.
// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;
// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;

// Shock capturing.
bool SHOCK_CAPTURING = false;
// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;

// For saving/loading of solution.
bool REUSE_SOLUTION = true;

// Initial polynomial degree.        
const int P_INIT = 0;                                           
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM_VERTEX = 3;                
// Number of initial mesh refinements towards the profile.
const int INIT_REF_NUM_BOUNDARY_ANISO = 4;        
// CFL value.
double CFL_NUMBER = 0.05;                                
// Initial time step.
double time_step = 1E-6;                                

// Adaptivity.
// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 5;

// Number of mesh refinements between two unrefinements.
// The mesh is not unrefined unless there has been a refinement since
// last unrefinement.
int REFINEMENT_COUNT = 0;                         

// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies (see below).
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

// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
CandList CAND_LIST = H2D_HP_ANISO;                

// Maximum polynomial degree used. -1 for unlimited.
const int MAX_P_ORDER = -1;                       

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

// Stopping criterion for adaptivity (rel. error tolerance between the
// fine mesh and coarse mesh solution in percent).
const double ERR_STOP = 1E-4;                     

// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 6200;                   

// Matrix solver for orthogonal projections: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 7142.8571428571428571428571428571;         
// Inlet density (dimensionless).   
const double RHO_EXT = 1.0;                       
// Inlet x-velocity (dimensionless).
const double V1_EXT = 0.01;       
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;                        
// Kappa.
const double KAPPA = 1.4;                         

// Boundary markers.
const std::string BDY_INLET = "Inlet";
const std::string BDY_OUTLET = "Outlet";
const std::string BDY_SOLID_WALL_PROFILE = "Solid Profile";
const std::string BDY_SOLID_WALL = "Solid";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("domain-arcs.xml", mesh);

  mesh->refine_towards_boundary(BDY_SOLID_WALL_PROFILE, INIT_REF_NUM_BOUNDARY_ANISO, true, true);
  mesh->refine_towards_vertex(0, INIT_REF_NUM_VERTEX, true);

  MeshView m;
  m.show(mesh);
  m.wait_for_close();

SpaceSharedPtr<double>space_rho(mesh, P_INIT);
  L2Space<double>space_rho_v_x(mesh, P_INIT);
  L2Space<double>space_rho_v_y(mesh, P_INIT);
  L2Space<double>space_e(new // Initialize boundary condition types and spaces with default shapesets.
  L2Space<double>(mesh, P_INIT));
  int ndof = Space<double>::get_num_dofs({&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e});
  Hermes::Mixins::Loggable::Static::info("Initial coarse ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  ConstantSolution<double> sln_rho(mesh, RHO_EXT);
  ConstantSolution<double> sln_rho_v_x(mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double> sln_rho_v_y(mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double> sln_e(mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  ConstantSolution<double> prev_rho(mesh, RHO_EXT);
  ConstantSolution<double> prev_rho_v_x(mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double> prev_rho_v_y(mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double> prev_e(mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  Solution<double> rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e;

  // Initialize weak formulation.
  std::vector<std::string> solid_wall_markers;
  solid_wall_markers.push_back(BDY_SOLID_WALL);
  solid_wall_markers.push_back(BDY_SOLID_WALL_PROFILE);
  
  std::vector<std::string> inlet_markers;
  inlet_markers.push_back(BDY_INLET);
  std::vector<std::string> outlet_markers;
  outlet_markers.push_back(BDY_OUTLET);

  EulerEquationsWeakFormSemiImplicit wf(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT,solid_wall_markers, 
    inlet_markers, outlet_markers, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number({&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e}, KAPPA);
  PressureFilter pressure({&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e}, KAPPA);
  EntropyFilter entropy({&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e}, KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  OrderView space_view("Space", new WinGeom(700, 400, 600, 300));
  
  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, MAX_P_ORDER);
  selector.set_error_weights(1.0, 1.0, 1.0);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Look for a saved solution on the disk.
  CalculationContinuity<double> continuity(CalculationContinuity<double>::onlyTime);
  int iteration = 0; double t = 0;
  bool loaded_now = false;

  if(REUSE_SOLUTION && continuity.have_record_available())
  {
    continuity.get_last_record()->load_mesh(mesh);
SpaceSharedPtr<double> *> spaceVector = continuity.get_last_record()->load_spaces(new   std::vector<Space<double>({mesh, mesh, mesh, mesh}));
    space_rho.copy(spaceVector[0], mesh);
    space_rho_v_x.copy(spaceVector[1], mesh);
    space_rho_v_y.copy(spaceVector[2], mesh);
    space_e.copy(spaceVector[3], mesh);
    continuity.get_last_record()->load_time_step_length(time_step);
    t = continuity.get_last_record()->get_time() + time_step;
    iteration = continuity.get_num() * EVERY_NTH_STEP + 1;
    loaded_now = true;
  }

  // Time stepping loop.
  for(; t < 5.0; t += time_step)
  {
    CFL.set_number(CFL_NUMBER + (t/5.0) * 10.0);
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0) 
    {
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      REFINEMENT_COUNT = 0;
      
      space_rho.unrefine_all_mesh_elements(true);
      
      space_rho.adjust_element_order(-1, P_INIT);
      space_rho_v_x.adjust_element_order(-1, P_INIT);
      space_rho_v_y.adjust_element_order(-1, P_INIT);
      space_e.adjust_element_order(-1, P_INIT);
    }

    // Adaptivity loop:
    int as = 1; 
    int ndofs_prev = 0;
    bool done = false;
    do
    {
      Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      int order_increase = 1;

      Mesh::ReferenceMeshCreator refMeshCreatorFlow(mesh);
      MeshSharedPtr ref_mesh_flow = refMeshCreatorFlow.create_ref_mesh();

      Space<double>::ReferenceSpaceCreator refSpaceCreatorRho(space_rho, ref_mesh_flow, order_increase);
SpaceSharedPtr<double>* ref_space_rho = refSpaceCreatorRho.create_ref_space(new     Space<double>());
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVx(space_rho_v_x, ref_mesh_flow, order_increase);
SpaceSharedPtr<double>* ref_space_rho_v_x = refSpaceCreatorRhoVx.create_ref_space(new     Space<double>());
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVy(space_rho_v_y, ref_mesh_flow, order_increase);
SpaceSharedPtr<double>* ref_space_rho_v_y = refSpaceCreatorRhoVy.create_ref_space(new     Space<double>());
      Space<double>::ReferenceSpaceCreator refSpaceCreatorE(space_e, ref_mesh_flow, order_increase);
SpaceSharedPtr<double>* ref_space_e = refSpaceCreatorE.create_ref_space(new     Space<double>());

SpaceSharedPtr<double>*> ref_spaces(new     std::vector<const Space<double>(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e));

      if(ndofs_prev != 0)
        if(Space<double>::get_num_dofs(ref_spaces) == ndofs_prev)
          selector.set_error_weights(2.0 * selector.get_error_weight_h(), 1.0, 1.0);
        else
          selector.set_error_weights(1.0, 1.0, 1.0);

      ndofs_prev = Space<double>::get_num_dofs(ref_spaces);

      // Project the previous time level solution onto the new fine mesh.
      Hermes::Mixins::Loggable::Static::info("Projecting the previous time level solution onto the new fine mesh.");
      if(loaded_now)
      {
        loaded_now = false;

        continuity.get_last_record()->load_solutions({&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e}, 
SpaceSharedPtr<double> *>(new         std::vector<Space<double>(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e)));
      }
      else
      {
      OGProjection<double>::project_global(ref_spaces,{&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e}, 
        std::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), std::vector<Hermes::Hermes2D::NormType>());
      }

      // Report NDOFs.
      Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d.", 
        Space<double>::get_num_dofs({&space_rho, &space_rho_v_x, 
        space_rho_v_y, &space_e}), Space<double>::get_num_dofs(ref_spaces));

      // Assemble the reference problem.
      Hermes::Mixins::Loggable::Static::info("Solving on reference mesh.");
      DiscreteProblem<double> dp(wf, ref_spaces);

      SparseMatrix<double>* matrix = create_matrix<double>();
      Vector<double>* rhs = create_vector<double>();
      Hermes::Solvers::LinearMatrixSolver<double>* solver = create_linear_solver<double>( matrix, rhs);

      wf.set_current_time_step(time_step);

      // Assemble the stiffness matrix and rhs.
      Hermes::Mixins::Loggable::Static::info("Assembling the stiffness matrix and right-hand side vector.");
      dp.assemble(matrix, rhs);

      // Solve the matrix problem.
      Hermes::Mixins::Loggable::Static::info("Solving the matrix problem.");
      if(solver->solve())
        if(!SHOCK_CAPTURING)
          Solution<double>::vector_to_solutions(solver->get_sln_vector(), ref_spaces, 
          std::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
        else
        {      
          FluxLimiter flux_limiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), ref_spaces, true);
          
          flux_limiter.limit_second_orders_according_to_detector({&space_rho, &space_rho_v_x, 
            &space_rho_v_y, &space_e});
          
          flux_limiter.limit_according_to_detector({&space_rho, &space_rho_v_x, 
            &space_rho_v_y, &space_e});

          flux_limiter.get_limited_solutions({&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e});
        }
      else
        throw Hermes::Exceptions::Exception("Matrix solver failed.\n");
      
      // Project the fine mesh solution onto the coarse mesh.
      Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
      OGProjection<double>::project_global({&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e},{&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e}, 
        std::vector<Solution<double>*>(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e), 
        std::vector<NormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
      Adapt<double>* adaptivity = new Adapt<double>({&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e},{HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM});
      double err_est_rel_total = adaptivity->calc_err_est({sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e},
        std::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e)) * 100;

      CFL.calculate_semi_implicit({rsln_rho, rsln_rho_v_x, &rsln_rho_v_y, &rsln_e}, ref_space_rho->get_mesh(), time_step);

      // Report results.
      Hermes::Mixins::Loggable::Static::info("err_est_rel: %g%%", err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP)
        done = true;
      else
      {
        Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
        if (Space<double>::get_num_dofs({&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e}) >= NDOF_STOP) 
          done = true;
        else
        {
          REFINEMENT_COUNT++;
          done = adaptivity->adapt({&selector, &selector, &selector, &selector}, 
          THRESHOLD, STRATEGY, MESH_REGULARITY);
        }

        if(!done)
          as++;
      }

      // Visualization and saving on disk.
      if(done && (iteration - 1) % EVERY_NTH_STEP == 0 && iteration > 1)
      {
        // Hermes visualization.
        if(HERMES_VISUALIZATION) 
        {
          Mach_number.reinit();
          pressure.reinit();
          entropy.reinit();
          pressure_view.show(&pressure);
          entropy_production_view.show(&entropy);
          Mach_number_view.show(&Mach_number);
          pressure_view.save_numbered_screenshot("Pressure-%u.bmp", iteration - 1, true);
          Mach_number_view.save_numbered_screenshot("Mach-%u.bmp", iteration - 1, true);
        }
        // Output solution in VTK format.
        if(VTK_VISUALIZATION) 
        {
          pressure.reinit();
          Mach_number.reinit();
          Linearizer lin;
          char filename[40];
          sprintf(filename, "Pressure-%i.vtk", iteration - 1);
          lin.save_solution_vtk(pressure, filename, "Pressure", false);
          sprintf(filename, "Mach number-%i.vtk", iteration - 1);
          lin.save_solution_vtk(Mach_number, filename, "MachNumber", false);
        }
      }

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
    }
    while (done == false);

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&rsln_rho);
    prev_rho_v_x.copy(&rsln_rho_v_x);
    prev_rho_v_y.copy(&rsln_rho_v_y);
    prev_e.copy(&rsln_e);
  }

  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();

  return 0;
}
