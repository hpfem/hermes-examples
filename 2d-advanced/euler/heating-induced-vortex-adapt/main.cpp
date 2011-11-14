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
// Domain: A square, see file square.mesh.
//
// BC: Solid walls, inlet, no outlet.
//
// IC: Constant state identical to inlet, only with higher pressure.
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
bool SHOCK_CAPTURING = true;
// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;

// For saving/loading of solution.
bool REUSE_SOLUTION = false;

// Initial polynomial degree.      
const int P_INIT = 0;                                             
// Number of initial uniform mesh refinements.     
const int INIT_REF_NUM = 1;                                               
// CFL value.
double CFL_NUMBER = 0.5;                                
 // Initial time step.
double time_step = 1E-6;                         

// Adaptivity.
const int NDOFS_MIN = 7000;

// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 5;

// Number of mesh refinements between two unrefinements.
// The mesh is not unrefined unless there has been a refinement since
// last unrefinement.
int REFINEMENT_COUNT = 0;                         

// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies (see below).
const double THRESHOLD = 0.3;                     

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
// See User Documentation for details.
const int MAX_P_ORDER = 1;                       

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

// Stopping criterion for adaptivity.
double ERR_STOP = 5.0;                     

// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 100000;                   

// Matrix solver for orthogonal projections: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 2.0;           
// Initial pressure (dimensionless).
const double P_INITIAL_HIGH = 1.5;  
// Initial pressure (dimensionless).
const double P_INITIAL_LOW = 1.0;   
// Inlet density (dimensionless).   
const double RHO_EXT = 1.0;                       
// Initial density (dimensionless).   
const double RHO_INITIAL_HIGH = 0.5;
// Initial density (dimensionless).   
const double RHO_INITIAL_LOW = 0.3; 
// Inlet x-velocity (dimensionless).
const double V1_EXT = 0.0;          
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;                        
// Kappa.
const double KAPPA = 1.4;                         

// Boundary markers.
const std::string BDY_INLET = "Inlet";
const std::string BDY_SOLID_WALL = "Solid";

// Area (square) size.
// Must be in accordance with the mesh file.
const double MESH_SIZE = 3.0;

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) 
    mesh.refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  L2Space<double>space_rho(&mesh, P_INIT);
  L2Space<double>space_rho_v_x(&mesh, P_INIT);
  L2Space<double>space_rho_v_y(&mesh, P_INIT);
  L2Space<double>space_e(&mesh, P_INIT);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  InitialSolutionLinearProgress sln_rho(&mesh, RHO_INITIAL_HIGH, RHO_INITIAL_LOW, MESH_SIZE);
  ConstantSolution<double> sln_rho_v_x(&mesh, 0.0);
  ConstantSolution<double> sln_rho_v_y(&mesh, 0.0);
  InitialSolutionLinearProgress sln_e(&mesh, QuantityCalculator::calc_energy(RHO_INITIAL_HIGH, RHO_INITIAL_HIGH * V1_EXT, RHO_INITIAL_HIGH * V2_EXT, P_INITIAL_HIGH, KAPPA), QuantityCalculator::calc_energy(RHO_INITIAL_LOW, RHO_INITIAL_LOW * V1_EXT, RHO_INITIAL_LOW * V2_EXT, P_INITIAL_LOW, KAPPA), MESH_SIZE);

  InitialSolutionLinearProgress prev_rho(&mesh, RHO_INITIAL_HIGH, RHO_INITIAL_LOW, MESH_SIZE);
  ConstantSolution<double> prev_rho_v_x(&mesh, 0.0);
  ConstantSolution<double> prev_rho_v_y(&mesh, 0.0);
  InitialSolutionLinearProgress prev_e(&mesh, QuantityCalculator::calc_energy(RHO_INITIAL_HIGH, RHO_INITIAL_HIGH * V1_EXT, RHO_INITIAL_HIGH * V2_EXT, P_INITIAL_HIGH, KAPPA), QuantityCalculator::calc_energy(RHO_INITIAL_LOW, RHO_INITIAL_LOW * V1_EXT, RHO_INITIAL_LOW * V2_EXT, P_INITIAL_LOW, KAPPA), MESH_SIZE);

  InitialSolutionLinearProgress rsln_rho(&mesh, RHO_INITIAL_HIGH, RHO_INITIAL_LOW, MESH_SIZE);
  ConstantSolution<double> rsln_rho_v_x(&mesh, 0.0);
  ConstantSolution<double> rsln_rho_v_y(&mesh, 0.0);
  InitialSolutionLinearProgress rsln_e(&mesh, QuantityCalculator::calc_energy(RHO_INITIAL_HIGH, RHO_INITIAL_HIGH * V1_EXT, RHO_INITIAL_HIGH * V2_EXT, P_INITIAL_HIGH, KAPPA), QuantityCalculator::calc_energy(RHO_INITIAL_LOW, RHO_INITIAL_LOW * V1_EXT, RHO_INITIAL_LOW * V2_EXT, P_INITIAL_LOW, KAPPA), MESH_SIZE);

  // Numerical flux.
  VijayasundaramNumericalFlux num_flux(KAPPA);
  
  // For saving to the disk.
  Continuity<double> continuity(Continuity<double>::onlyNumber);

  // Initialize weak formulation.
  EulerEquationsWeakFormSemiImplicitMultiComponent wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL, BDY_SOLID_WALL, 
    BDY_INLET, "Outlet marker not used", &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA, RHO_INITIAL_HIGH, P_INITIAL_HIGH);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  VectorView velocity_view("Velocity", new WinGeom(700, 400, 600, 300));
  
  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, MAX_P_ORDER);
  selector.set_error_weights(1.0, 1.0, 1.0);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Time stepping loop.
  int iteration = 0; double t = 0;
  for(; t < 7.0; t += time_step)
  {
    if(t > 0.25)
      ERR_STOP = 2.3;

    info("---- Time step %d, time %3.5f.", iteration++, t);

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0) 
    {
      info("Global mesh derefinement.");
      REFINEMENT_COUNT = 0;
      
      space_rho.unrefine_all_mesh_elements(true);
      
      space_rho.adjust_element_order(-1, P_INIT);
      space_rho_v_x.copy_orders(&space_rho);
      space_rho_v_y.copy_orders(&space_rho);
      space_e.copy_orders(&space_rho);
    }

    // Adaptivity loop:
    int as = 1; 
    int ndofs_prev = 0;
    bool done = false;
    do
    {
      info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      int order_increase = 1;

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
      info("Projecting the previous time level solution onto the new fine mesh.");
      OGProjection<double>::project_global(ref_spaces_const, Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
        Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), matrix_solver, Hermes::vector<Hermes::Hermes2D::ProjNormType>(), iteration > 1);

      
      // Report NDOFs.
      info("ndof_coarse: %d, ndof_fine: %d.", 
        Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e)), Space<double>::get_num_dofs(ref_spaces_const));

      // Assemble the reference problem.
      info("Solving on reference mesh.");
      DiscreteProblem<double> dp(&wf, ref_spaces_const);

      SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
      Vector<double>* rhs = create_vector<double>(matrix_solver);
      LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver, matrix, rhs);

      wf.set_time_step(time_step);

      dp.assemble(matrix, rhs);

      // Solve the matrix problem.
      info("Solving the matrix problem.");
      if(solver->solve())
        if(!SHOCK_CAPTURING)
          Solution<double>::vector_to_solutions(solver->get_sln_vector(), ref_spaces_const, 
          Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
        else
        {      
          FluxLimiter flux_limiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), ref_spaces_const);
          
          flux_limiter.limit_second_orders_according_to_detector(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
            &space_rho_v_y, &space_e));
          
          flux_limiter.limit_according_to_detector(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
            &space_rho_v_y, &space_e));

          flux_limiter.get_limited_solutions(Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
        }
      else
        error ("Matrix solver failed.\n");

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection<double>::project_global(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), 
        Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), matrix_solver, 
        Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e),
        Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e)) * 100;

      CFL.calculate_semi_implicit(Hermes::vector<Solution<double> *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), (ref_spaces_const)[0]->get_mesh(), time_step);

      // Report results.
      info("err_est_rel: %g%%", err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP && Space<double>::get_num_dofs(ref_spaces_const) > NDOFS_MIN)
        done = true;
      else
      {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector, &selector, &selector), 
          THRESHOLD, STRATEGY, MESH_REGULARITY);

        REFINEMENT_COUNT++;
        if (Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e)) >= NDOF_STOP) 
          done = true;
        else
          as++;
      }

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
      if(!done)
        for(unsigned int i = 0; i < ref_spaces->size(); i++)
          delete (ref_spaces_const)[i];
    }
    while (done == false);

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&rsln_rho);
    prev_rho_v_x.copy(&rsln_rho_v_x);
    prev_rho_v_y.copy(&rsln_rho_v_y);
    prev_e.copy(&rsln_e);

    delete rsln_rho.get_mesh();
    rsln_rho.own_mesh = false;
    delete rsln_rho_v_x.get_mesh();
    rsln_rho_v_x.own_mesh = false;
    delete rsln_rho_v_y.get_mesh();
    rsln_rho_v_y.own_mesh = false;
    delete rsln_e.get_mesh();
    rsln_e.own_mesh = false;

      // Visualization and saving on disk.
    if((iteration - 1) % EVERY_NTH_STEP == 0)
    {
      continuity.add_record((unsigned int)(iteration - 1));
      continuity.get_last_record()->save_mesh(prev_rho.get_mesh());
      continuity.get_last_record()->save_space(const_cast<Space<double>*>(prev_rho.get_space()));
      continuity.get_last_record()->save_time_step_length(time_step);

      // Hermes visualization.
      if(HERMES_VISUALIZATION)
      {        
        pressure.reinit();
        pressure_view.show(&pressure);
        velocity_view.show(&rsln_rho_v_x, &rsln_rho_v_y);
        
        pressure_view.save_numbered_screenshot("pressure %i.bmp", iteration);
        velocity_view.save_numbered_screenshot("Velocity %i.bmp", iteration);
      }
      // Output solution in VTK format.
      if(VTK_VISUALIZATION)
      {
        pressure.reinit();
        Linearizer lin;
        char filename[40];
        sprintf(filename, "Pressure-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&pressure, filename, "Pressure", false);
        sprintf(filename, "VelocityX-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_rho_v_x, filename, "VelocityX", false);
        sprintf(filename, "VelocityY-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_rho_v_y, filename, "VelocityY", false);
        sprintf(filename, "Rho-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_rho, filename, "Rho", false);
      }
    }
  }

  pressure_view.close();
  velocity_view.close();

  return 0;
}
