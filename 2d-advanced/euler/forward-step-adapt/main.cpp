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
// Domain: forward facing step, see mesh file ffs.mesh
//
// BC: Solid walls, inlet, no outlet.
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
bool SHOCK_CAPTURING = true;
// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;

// For saving/loading of solution.
bool REUSE_SOLUTION = true;

// Initial polynomial degree.     
const int P_INIT = 0;                                              
// Number of initial uniform mesh refinements. 
const int INIT_REF_NUM = 1;                                             
// Number of initial localized mesh refinements.         
const int INIT_REF_NUM_STEP = 1;                                 
// CFL value.
double CFL_NUMBER = 0.5;                         
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
const int NDOF_STOP = 7000;                   

// Matrix solver for orthogonal projections: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Equation parameters.
const double P_EXT = 1.0;         // Exterior pressure (dimensionless).
const double RHO_EXT = 1.4;       // Inlet density (dimensionless).   
const double V1_EXT = 3.0;        // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;        // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;         // Kappa.

// Boundary markers.
const std::string BDY_SOLID_WALL_BOTTOM = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_TOP = "3";
const std::string BDY_INLET = "4";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

// Criterion for mesh refinement.
int refinement_criterion(Element* e)
{
  if(e->vn[2]->y <= 0.4 && e->vn[1]->x <= 0.6)
    return 0;
  else
    return -1;
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, base_mesh;
  MeshReaderH2D mloader;
  mloader.load("ffs.mesh", &base_mesh);

  base_mesh.refine_by_criterion(refinement_criterion, INIT_REF_NUM_STEP);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    base_mesh.refine_all_elements(0, true);

  mesh.copy(&base_mesh);

  // Initialize boundary condition types and spaces with default shapesets.
  L2Space<double>space_rho(&mesh, P_INIT);
  L2Space<double>space_rho_v_x(&mesh, P_INIT);
  L2Space<double>space_rho_v_y(&mesh, P_INIT);
  L2Space<double>space_e(&mesh, P_INIT);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  ConstantSolution<double> sln_rho(&mesh, RHO_EXT);
  ConstantSolution<double> sln_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double> sln_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double> sln_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  ConstantSolution<double> prev_rho(&mesh, RHO_EXT);
  ConstantSolution<double> prev_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double> prev_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double> prev_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  Solution<double> rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e;

  // Numerical flux.
  VijayasundaramNumericalFlux num_flux(KAPPA);

  // Initialize weak formulation.
  EulerEquationsWeakFormSemiImplicitMultiComponent wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP, 
    BDY_INLET, BDY_OUTLET, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));

  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, MAX_P_ORDER);
  selector.set_error_weights(1.0, 1.0, 1.0);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Look for a saved solution on the disk.
  Continuity<double> continuity(Continuity<double>::onlyTime);
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
    iteration = (continuity.get_num() - 1) * EVERY_NTH_STEP + 1;
    loaded_now = true;
  }

  // Time stepping loop.
  for(; t < 4.5; t += time_step)
  {
    if(t > 0.3)
      ERR_STOP = 2.5;

    CFL.set_number(CFL_NUMBER + (t/4.5) * 1.0);
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
      if(loaded_now)
      {
        loaded_now = false;

        continuity.get_last_record()->load_solutions(Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
          Hermes::vector<Space<double> *>((*ref_spaces)[0], (*ref_spaces)[1], (*ref_spaces)[2], (*ref_spaces)[3]));
      }
      else
      {
        OGProjection<double>::project_global(ref_spaces_const, Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
            Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), matrix_solver, Hermes::vector<Hermes::Hermes2D::ProjNormType>());
        if(iteration > std::max((int)(continuity.get_num() * EVERY_NTH_STEP + 2), 1) && as > 1)
        {
          delete rsln_rho.get_mesh();
          delete rsln_rho.get_space();
          rsln_rho.own_mesh = false;
          delete rsln_rho_v_x.get_mesh();
          delete rsln_rho_v_x.get_space();
          rsln_rho_v_x.own_mesh = false;
          delete rsln_rho_v_y.get_mesh();
          delete rsln_rho_v_y.get_space();
          rsln_rho_v_y.own_mesh = false;
          delete rsln_e.get_mesh();
          delete rsln_e.get_space();
          rsln_e.own_mesh = false;
        }
      }

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

      // Assemble the stiffness matrix and rhs.
      info("Assembling the stiffness matrix and right-hand side vector.");
      dp.assemble(matrix, rhs);

      // Solve the matrix problem.
      info("Solving the matrix problem.");
      if(solver->solve())
        if(!SHOCK_CAPTURING)
          Solution<double>::vector_to_solutions(solver->get_sln_vector(), ref_spaces_const, 
          Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
        else
        {      
          FluxLimiter flux_limiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), ref_spaces_const, true);

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

      CFL.calculate_semi_implicit(Hermes::vector<Solution<double> *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), (*ref_spaces)[0]->get_mesh(), time_step);

      // Report results.
      info("err_est_rel: %g%%", err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP)
        done = true;
      else
      {
        info("Adapting coarse mesh.");
        if (Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e)) >= NDOF_STOP) 
          done = true;
        else
        {
          REFINEMENT_COUNT++;
          done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector, &selector, &selector), 
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
          entropy.reinit();
          Linearizer lin;
          char filename[40];
          sprintf(filename, "Pressure-%i.vtk", iteration - 1);
          lin.save_solution_vtk(&pressure, filename, "Pressure", false);
          sprintf(filename, "Mach number-%i.vtk", iteration - 1);
          lin.save_solution_vtk(&Mach_number, filename, "MachNumber", false);
          sprintf(filename, "Entropy-%i.vtk", iteration - 1);
          lin.save_solution_vtk(&entropy, filename, "Entropy", false);
        }
        // Save a current state on the disk.
        if(iteration > 1)
        {
          continuity.add_record(t);
          continuity.get_last_record()->save_mesh(&mesh);
          continuity.get_last_record()->save_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
            &space_rho_v_y, &space_e));
          continuity.get_last_record()->save_solutions(Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
          continuity.get_last_record()->save_time_step_length(time_step);
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

    delete rsln_rho.get_mesh();
    delete rsln_rho.get_space();
    rsln_rho.own_mesh = false;
    delete rsln_rho_v_x.get_mesh();
    delete rsln_rho_v_x.get_space();
    rsln_rho_v_x.own_mesh = false;
    delete rsln_rho_v_y.get_mesh();
    delete rsln_rho_v_y.get_space();
    rsln_rho_v_y.own_mesh = false;
    delete rsln_e.get_mesh();
    delete rsln_e.get_space();
    rsln_e.own_mesh = false;
  }

  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();

  return 0;
}
