#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// This example solves the compressible Euler equations using FVM and automatic h-adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: GAMM channel, see mesh file GAMM-channel.mesh
//
// BC: Normal velocity component is zero on solid walls.
//     Subsonic state prescribed on inlet and outlet.
//
// IC: Constant subsonic state identical to inlet. 
//
// The following parameters can be changed:
// Visualization.
const bool HERMES_VISUALIZATION = true;           // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;              // Set to "true" to enable VTK output.
const unsigned int EVERY_NTH_STEP = 1;            // Set visual output for every nth step.

// Shock capturing.
bool SHOCK_CAPTURING = true;
// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;

const int P_INIT = 0;                             // Initial polynomial degree.                      
const int INIT_REF_NUM = 2;                             // Number of initial uniform mesh refinements.                       
double CFL_NUMBER = 1.0;                                // CFL value.
double time_step = 1E-6;                          // Initial time step.

// Adaptivity.
// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 1;

// Number of mesh refinements between two unrefinements.
// The mesh is not unrefined unless there has been a refinement since
// last unrefinement.
int REFINEMENT_COUNT = 0;                         

// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies (see below).
const double THRESHOLD = 0.7;                     

// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
// More adaptive strategies can be created in adapt_ortho_h1.cpp.
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

// Default value is 1.0. This parameter influences the selection of
// cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double CONV_EXP = 1;                        

// Stopping criterion for adaptivity (rel. error tolerance between the
// fine mesh and coarse mesh solution in percent).
const double ERR_STOP = 1.0;                     

// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 100000;                   

// Matrix solver for orthogonal projections.
// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  

// Equation parameters.
const double P_EXT = 2.5;                         // Exterior pressure (dimensionless).
const double RHO_EXT = 1.0;                       // Inlet density (dimensionless).   
const double V1_EXT = 1.25;                       // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;                        // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;                         // Kappa.

// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("GAMM-channel.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) 
    mesh.refine_all_elements(0, true);
  //mesh.refine_towards_boundary(BDY_SOLID_WALL_BOTTOM, 2);

  // Initialize boundary condition types and spaces with default shapesets.
  L2Space<double>space_rho(&mesh, P_INIT);
  L2Space<double>space_rho_v_x(&mesh, P_INIT);
  L2Space<double>space_rho_v_y(&mesh, P_INIT);
  L2Space<double>space_e(&mesh, P_INIT);

  // Initialize solutions, set initial conditions.
  InitialSolutionEulerDensity sln_rho(&mesh, RHO_EXT);
  InitialSolutionEulerDensityVelX sln_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY sln_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy sln_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  InitialSolutionEulerDensity prev_rho(&mesh, RHO_EXT);
  InitialSolutionEulerDensityVelX prev_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY prev_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy prev_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  Solution<double> rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e;

  // Numerical flux.
  OsherSolomonNumericalFlux num_flux(KAPPA);

  // Initialize weak formulation.
  EulerEquationsWeakFormSemiImplicitMultiComponent wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP, 
    BDY_INLET, BDY_OUTLET, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, (P_INIT == 0 && CAND_LIST == H2D_H_ANISO));

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));

  //OrderView<double> order_view_coarse("Orders - coarse", new WinGeom(700, 350, 600, 300));
  //OrderView<double> order_view_fine("Orders - fine", new WinGeom(700, 700, 600, 300));

  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, MAX_P_ORDER);
  selector.set_error_weights(1.0, 1.0, 1.0);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 3.0; t += time_step)
  {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0) {
      REFINEMENT_COUNT = 0;
      info("Global mesh derefinement.");
      space_rho.unrefine_all_mesh_elements();
      if(CAND_LIST == H2D_HP_ANISO)
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
      int order_increase = 0;
      if(CAND_LIST == H2D_HP_ANISO)
        order_increase = 1;

      Hermes::vector<Space<double> *>* ref_spaces = Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), order_increase);

      if(ndofs_prev != 0)
        if(Space<double>::get_num_dofs(*ref_spaces) == ndofs_prev)
          selector.set_error_weights(2.0 * selector.get_error_weight_h(), 1.0, 1.0);
        else
          selector.set_error_weights(1.0, 1.0, 1.0);

      ndofs_prev = Space<double>::get_num_dofs(*ref_spaces);

      // Project the previous time level solution onto the new fine mesh.
      info("Projecting the previous time level solution onto the new fine mesh.");
      OGProjection<double>::project_global(*ref_spaces, Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
        Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), matrix_solver_type, Hermes::vector<Hermes::Hermes2D::ProjNormType>(), iteration > 1);

      // Report NDOFs.
      info("ndof_coarse: %d, ndof_fine: %d.", 
        Space<double>::get_num_dofs(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e)), Space<double>::get_num_dofs(*ref_spaces));

      // Assemble the reference problem.
      info("Solving on reference mesh.");
      bool is_linear = true;
      DiscreteProblem<double> dp(&wf, *ref_spaces);

      SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver_type);
      Vector<double>* rhs = create_vector<double>(matrix_solver_type);
      LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver_type, matrix, rhs);

      wf.set_time_step(time_step);

      dp.assemble(matrix, rhs);

      dp.get_last_profiling_output(std::cout);

      // Solve the matrix problem.
      info("Solving the matrix problem.");
      if(solver->solve())
        if(!SHOCK_CAPTURING)
          Solution<double>::vector_to_solutions(solver->get_sln_vector(), *ref_spaces, 
          Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
        else
        {      
          FluxLimiter flux_limiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), *ref_spaces);
          
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
      OGProjection<double>::project_global(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), 
        Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), matrix_solver_type, 
        Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e),
        Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e)) * 100;

      CFL.calculate_semi_implicit(Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), &mesh, time_step);

      // Report results.
      info("err_est_rel: %g%%", err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP)
      {
        done = true;

        // Visualization of the refined space before it gets deleted.
        if((iteration - 1) % EVERY_NTH_STEP == 0) {
          // Hermes visualization.
          if(HERMES_VISUALIZATION)
          {
            //order_view_coarse.show(&space_rho);
            // order_view_fine.show((*ref_spaces)[0]);
          }
          if(VTK_VISUALIZATION)
          {            
            Orderizer ord;
            char filename[40];
            sprintf(filename, "Orders-coarse-%i.vtk", iteration - 1);
            ord.save_orders_vtk(&space_rho, filename);
            sprintf(filename, "Orders-fine-%i.vtk", iteration - 1);
            ord.save_orders_vtk((*ref_spaces)[0], filename);
          }
        }
      }
      else
      {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector, &selector, &selector), 
          THRESHOLD, STRATEGY, MESH_REGULARITY);

        REFINEMENT_COUNT++;
        if (Space<double>::get_num_dofs(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
      for(unsigned int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i];
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

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0)
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
      }
      // Output solution in VTK format.
      if(VTK_VISUALIZATION)
      {
        pressure.reinit();
        Mach_number.reinit();
        Linearizer lin_pressure(&pressure);
        char filename[40];
        sprintf(filename, "pressure-3D-%i.vtk", iteration - 1);
        lin_pressure.save_solution_vtk(filename, "Pressure", true);
        Linearizer lin_mach(&Mach_number);
        sprintf(filename, "Mach number-3D-%i.vtk", iteration - 1);
        lin_mach.save_solution_vtk(filename, "MachNumber", true);
      }
    }
  }

  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();

  return 0;
}
