#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// This example solves the compressible Euler equations using 
// Discontinuous Galerkin method of higher order with adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: GAMM channel, see mesh file GAMM-channel.mesh
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
enum shockCapturingType
{
  FEISTAUER,
  KUZMIN,
  KRIVODONOVA
};
bool SHOCK_CAPTURING = true;
shockCapturingType SHOCK_CAPTURING_TYPE = FEISTAUER;
// Quantitative parameter of the discontinuity detector in case of Krivodonova.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;
// Quantitative parameter of the shock capturing in case of Feistauer.
const double NU_1 = 0.1;
const double NU_2 = 0.1;

// For saving/loading of solution.
bool REUSE_SOLUTION = true;

// Initial polynomial degree. 
const int P_INIT = 0;                                                  
// Number of initial uniform mesh refinements.  
const int INIT_REF_NUM = 2;
// CFL value.
double CFL_NUMBER = 0.5;                          
// Initial time step.
double time_step_n = 1E-6;                        
// Initial time step.
double time_step_n_minus_one = 1E-6;              

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
const double ERR_STOP = 0.55;                     

// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 100000;                   

// Matrix solver for orthogonal projections: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 2.5;                         
// Inlet density (dimensionless).   
const double RHO_EXT = 1.0;                       
// Inlet x-velocity (dimensionless).
const double V1_EXT = 1.25;                       
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;                        
// Kappa.
const double KAPPA = 1.4;                         

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

  ConstantSolution<double>* prev_rho = new ConstantSolution<double>(&mesh, RHO_EXT);
  ConstantSolution<double>* prev_rho_v_x = new ConstantSolution<double>(&mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double>* prev_rho_v_y = new ConstantSolution<double>(&mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double>* prev_e = new ConstantSolution<double>(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  ConstantSolution<double>* prev_rho2 = new ConstantSolution<double>(&mesh, RHO_EXT);
  ConstantSolution<double>* prev_rho_v_x2 = new ConstantSolution<double>(&mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double>* prev_rho_v_y2 = new ConstantSolution<double>(&mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double>* prev_e2 = new ConstantSolution<double>(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  ConstantSolution<double>* rsln_rho = new ConstantSolution<double>(&mesh, RHO_EXT);
  ConstantSolution<double>* rsln_rho_v_x = new ConstantSolution<double>(&mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double>* rsln_rho_v_y = new ConstantSolution<double>(&mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double>* rsln_e = new ConstantSolution<double>(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction<double>*>(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction<double>*>(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), KAPPA, RHO_EXT, P_EXT);

  // Numerical flux.
  VijayasundaramNumericalFlux num_flux(KAPPA);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  ScalarView s1("prev_rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("prev_rho_v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("prev_rho_v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 400, 600, 300));

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
    continuity.get_last_record()->load_time_step_length(time_step_n);
    continuity.get_last_record()->load_time_step_length_n_minus_one(time_step_n_minus_one);
    t = continuity.get_last_record()->get_time();
    iteration = continuity.get_num();
    loaded_now = true;
  }

  // Time stepping loop.
  for(; t < 3.5; t += time_step_n)
  {
    CFL.set_number(CFL_NUMBER + (t/3.5) * 1.0);
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
      int order_increase = CAND_LIST == H2D_HP_ANISO ? 1 : 0;

      Hermes::vector<Space<double> *>* ref_spaces = Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), order_increase);

      Hermes::vector<const Space<double> *> ref_spaces_const((*ref_spaces)[0], (*ref_spaces)[1], 
        (*ref_spaces)[2], (*ref_spaces)[3]);

      L2Space<double> refspace_stabilization((*ref_spaces)[0]->get_mesh(), 0);

      if(ndofs_prev != 0)
        if(Space<double>::get_num_dofs(ref_spaces_const) == ndofs_prev)
          selector.set_error_weights(2.0 * selector.get_error_weight_h(), 1.0, 1.0);
        else
          selector.set_error_weights(1.0, 1.0, 1.0);

      ndofs_prev = Space<double>::get_num_dofs(ref_spaces_const);

      // Project the previous time level solution onto the new fine mesh.
      info("Projecting the previous time level solution onto the new fine mesh.");
      if(iteration == 1)
      {
        delete prev_rho;
        delete prev_rho_v_x;
        delete prev_rho_v_y;
        delete prev_e;

        delete prev_rho2;
        delete prev_rho_v_x2;
        delete prev_rho_v_y2;
        delete prev_e2;

        prev_rho = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), RHO_EXT);
        prev_rho_v_x = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), RHO_EXT * V1_EXT);
        prev_rho_v_y = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), RHO_EXT * V2_EXT);
        prev_e = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
        
        prev_rho2 = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), RHO_EXT);
        prev_rho_v_x2 = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), RHO_EXT * V1_EXT);
        prev_rho_v_y2 = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), RHO_EXT * V2_EXT);
        prev_e2 = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
      }
      else if(loaded_now)
      {
        loaded_now = false;

        continuity.get_last_record()->load_solutions(Hermes::vector<Solution<double>*>(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, prev_rho2, prev_rho_v_x2, prev_rho_v_y2, prev_e2), 
            Hermes::vector<Space<double> *>((*ref_spaces)[0], (*ref_spaces)[1], (*ref_spaces)[2], (*ref_spaces)[3], (*ref_spaces)[0], (*ref_spaces)[1], (*ref_spaces)[2], (*ref_spaces)[3]));
      }
      else
      {
        OGProjection<double>::project_global(ref_spaces_const, Hermes::vector<Solution<double>*>(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 
          Hermes::vector<Solution<double>*>(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), matrix_solver, Hermes::vector<Hermes::Hermes2D::ProjNormType>(), iteration > 1);
        
        if(iteration == 2)
        {
          delete prev_rho2;
          delete prev_rho_v_x2;
          delete prev_rho_v_y2;
          delete prev_e2;

          prev_rho2 = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), RHO_EXT);
          prev_rho_v_x2 = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), RHO_EXT * V1_EXT);
          prev_rho_v_y2 = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), RHO_EXT * V2_EXT);
          prev_e2 = new ConstantSolution<double>((ref_spaces_const)[0]->get_mesh(), QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
        }
        else
          OGProjection<double>::project_global(ref_spaces_const, Hermes::vector<Solution<double>*>(prev_rho2, prev_rho_v_x2, prev_rho_v_y2, prev_e2), 
            Hermes::vector<Solution<double>*>(prev_rho2, prev_rho_v_x2, prev_rho_v_y2, prev_e2), matrix_solver, Hermes::vector<Hermes::Hermes2D::ProjNormType>());
      }

      // Initialize weak formulation.
      EulerEquationsWeakFormSemiImplicitMultiComponent2ndOrder wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP, 
        BDY_INLET, BDY_OUTLET, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, prev_rho2, prev_rho_v_x2, prev_rho_v_y2, prev_e2, (P_INIT == 0 && CAND_LIST == H2D_H_ANISO));

      EulerEquationsWeakFormStabilization wf_stabilization(prev_rho);

      if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
        wf.set_stabilization(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, prev_rho2, prev_rho_v_x2, prev_rho_v_y2, prev_e2, NU_1, NU_2);
      
      // Report NDOFs.
      info("ndof_coarse: %d, ndof_fine: %d.", 
        Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e)), Space<double>::get_num_dofs(ref_spaces_const));

      // Assemble the reference problem.
      info("Solving on reference mesh.");
      DiscreteProblem<double> dp(&wf, ref_spaces_const);
      DiscreteProblem<double> dp_stabilization(&wf_stabilization, &refspace_stabilization);
      bool* discreteIndicator = NULL;

      SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
      Vector<double>* rhs = create_vector<double>(matrix_solver);
      Vector<double>* rhs_stabilization = create_vector<double>(matrix_solver);
      LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver, matrix, rhs);

      if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
      {
        dp_stabilization.assemble(rhs_stabilization);
        if(discreteIndicator != NULL)
          delete [] discreteIndicator;
        discreteIndicator = new bool[refspace_stabilization.get_mesh()->get_max_element_id() + 1];
        for(unsigned int i = 0; i < refspace_stabilization.get_mesh()->get_max_element_id() + 1; i++)
          discreteIndicator[i] = false;
        Element* e;
        for_all_active_elements(e, refspace_stabilization.get_mesh())
        {
          AsmList<double> al;
          refspace_stabilization.get_element_assembly_list(e, &al);
          if(rhs_stabilization->get(al.get_dof()[0]) >= 1)
            discreteIndicator[e->id] = true;
        }
        wf.set_discreteIndicator(discreteIndicator);
      }

      // Set the current time step.
      wf.set_time_step(time_step_n, time_step_n_minus_one);

      // If the FE problem is in fact a FV problem.
      if(P_INIT == 0 && CAND_LIST == H2D_H_ANISO) 
        dp.set_fvm();

      dp.assemble(matrix, rhs);

      // Solve the matrix problem.
      info("Solving the matrix problem.");
      if(solver->solve())
      {
        if(iteration > 1)
        {
          prev_rho2->copy(prev_rho);
          prev_rho_v_x2->copy(prev_rho_v_x);
          prev_rho_v_y2->copy(prev_rho_v_y);
          prev_e2->copy(prev_e);
        }

        if(!SHOCK_CAPTURING || SHOCK_CAPTURING_TYPE == FEISTAUER)
        {
          Solution<double>::vector_to_solutions(solver->get_sln_vector(), ref_spaces_const, 
            Hermes::vector<Solution<double>*>(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e));
        }
        else
        {
          FluxLimiter* flux_limiter;
          if(SHOCK_CAPTURING_TYPE == KUZMIN)
            FluxLimiter flux_limiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), ref_spaces_const);
          else
            FluxLimiter flux_limiter(FluxLimiter::Krivodonova, solver->get_sln_vector(), ref_spaces_const);
          if(SHOCK_CAPTURING_TYPE == KUZMIN)
            flux_limiter->limit_second_orders_according_to_detector(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
            &space_rho_v_y, &space_e));

          flux_limiter->limit_according_to_detector(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
            &space_rho_v_y, &space_e));

          flux_limiter->get_limited_solutions(Hermes::vector<Solution<double>*>(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e));
        }
      }
      else
        error ("Matrix solver failed.\n");

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection<double>::project_global(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<Solution<double>*>(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), 
        Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), matrix_solver, 
        Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e),
        Hermes::vector<Solution<double>*>(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e)) * 100;

      time_step_n_minus_one = time_step_n;
      CFL.calculate_semi_implicit(Hermes::vector<Solution<double> *>(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), (ref_spaces_const)[0]->get_mesh(), time_step_n);

      // Report results.
      info("err_est_rel: %g%%", err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP)
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

      // Visualization and saving on disk.
      if(done && (iteration - 1) % EVERY_NTH_STEP == 0 && iteration > 1)
      {
        // Save a current state on the disk.
        if(iteration > 1)
        {
          continuity.add_record(t);
          continuity.get_last_record()->save_mesh(&mesh);
          continuity.get_last_record()->save_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
              &space_rho_v_y, &space_e));
          continuity.get_last_record()->save_solutions(Hermes::vector<Solution<double>*>(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, prev_rho2, prev_rho_v_x2, prev_rho_v_y2, prev_e2));
          continuity.get_last_record()->save_time_step_length(time_step_n);
          continuity.get_last_record()->save_time_step_length_n_minus_one(time_step_n_minus_one);
        }
  
        // Hermes visualization.
        if(HERMES_VISUALIZATION)
        {        
          Mach_number.reinit();
          pressure.reinit();
          entropy.reinit();
          pressure_view.show(&pressure, 1);
          entropy_production_view.show(&entropy, 1);
          Mach_number_view.show(&Mach_number, 1);

          pressure_view.save_numbered_screenshot("pressure %i.bmp", iteration);
          Mach_number_view.save_numbered_screenshot("Mach no %i.bmp", iteration);
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
      }

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete rhs_stabilization;
      delete adaptivity;
    }
    while (done == false);

    // Copy the solutions into the previous time level ones.
    prev_rho->copy(rsln_rho);
    prev_rho_v_x->copy(rsln_rho_v_x);
    prev_rho_v_y->copy(rsln_rho_v_y);
    prev_e->copy(rsln_e);

    delete rsln_rho->get_mesh();
    rsln_rho->own_mesh = false;
    delete rsln_rho_v_x->get_mesh();
    rsln_rho_v_x->own_mesh = false;
    delete rsln_rho_v_y->get_mesh();
    rsln_rho_v_y->own_mesh = false;
    delete rsln_e->get_mesh();
    rsln_e->own_mesh = false;
  }

  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();

  return 0;
}
