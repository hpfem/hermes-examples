#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

using namespace RefinementSelectors;
// This example solves the compressible Euler equations coupled with an advection-diffution equation
// using a basic piecewise-constant finite volume method for the flow and continuous FEM for the concentration
// being advected by the flow.
//
// Equations: Compressible Euler equations, perfect gas state equation, advection-diffusion equation.
//
// Domains: Various
//
// BC: Various.
//
// IC: Various.
//
// The following parameters can be changed:
// Semi-implicit scheme.
const bool SEMI_IMPLICIT = true;
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

// Stability for the concentration part.
double ADVECTION_STABILITY_CONSTANT = 1.0;
const double DIFFUSION_STABILITY_CONSTANT = 1.0;

// Polynomial degree for the Euler equations (for the flow).
const int P_INIT_FLOW = 0;                        
// Polynomial degree for the concentration.
const int P_INIT_CONCENTRATION = 1;               
// CFL value.
double CFL_NUMBER = 1.0;                          
// Initial and utility time step.
double time_step = 1E-5, util_time_step;          

// Adaptivity.
// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 5;                         
// Number of mesh refinements between two unrefinements.
// The mesh is not unrefined unless there has been a refinement since
// last unrefinement.
int REFINEMENT_COUNT_FLOW = 0;                         
// Number of mesh refinements between two unrefinements.
// The mesh is not unrefined unless there has been a refinement since
// last unrefinement.
int REFINEMENT_COUNT_CONCENTRATION = 0;                         
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.1;                     
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//  error is processed. If more elements have similar errors, refine
//  all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//  than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//  than THRESHOLD.
const int STRATEGY = 1;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST_FLOW = H2D_HP_ANISO, CAND_LIST_CONCENTRATION = H2D_HP_ANISO;     
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
// Stopping criterion for adaptivity.
double ERR_STOP_FLOW = 1.0;                 
// Stopping criterion for adaptivity.
double ERR_STOP_CONCENTRATION = 10.0;        
// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 100000;                     
// Matrix solver for orthogonal projections: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Number of initial uniform mesh refinements of the mesh for the flow.
unsigned int INIT_REF_NUM_FLOW = 3;               
// Number of initial uniform mesh refinements of the mesh for the concentration.
unsigned int INIT_REF_NUM_CONCENTRATION = 3;      
// Number of initial mesh refinements of the mesh for the concentration towards the 
// part of the boundary where the concentration is prescribed.
unsigned int INIT_REF_NUM_CONCENTRATION_BDY = 1;  

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
// Concentration on the boundary.
const double CONCENTRATION_EXT = 0.1;                  
// Start time of the concentration on the boundary.
const double CONCENTRATION_EXT_STARTUP_TIME = 0.0;     
// Diffusivity.
const double EPSILON = 0.001;                           

// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";
const std::string BDY_DIRICHLET_CONCENTRATION = "3";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  Hermes::vector<std::string> BDY_NATURAL_CONCENTRATION;
  BDY_NATURAL_CONCENTRATION.push_back("2");

  std::ofstream time_step_out("time_step");

  // Load the mesh.
  Mesh basemesh;
  MeshReaderH2D mloader;
  mloader.load("GAMM-channel.mesh", &basemesh);

  // Initialize the meshes.
  Mesh mesh_flow, mesh_concentration;
  mesh_flow.copy(&basemesh);
  mesh_concentration.copy(&basemesh);

  for(unsigned int i = 0; i < INIT_REF_NUM_CONCENTRATION; i++)
    mesh_concentration.refine_all_elements(0, true);

  mesh_concentration.refine_towards_boundary(BDY_DIRICHLET_CONCENTRATION, INIT_REF_NUM_CONCENTRATION_BDY, true);
  //mesh_flow.refine_towards_boundary(BDY_DIRICHLET_CONCENTRATION, INIT_REF_NUM_CONCENTRATION_BDY);

  for(unsigned int i = 0; i < INIT_REF_NUM_FLOW; i++)
    mesh_flow.refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  // For the concentration.
  EssentialBCs<double> bcs_concentration;

  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_DIRICHLET_CONCENTRATION, CONCENTRATION_EXT, CONCENTRATION_EXT_STARTUP_TIME));
  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_SOLID_WALL_TOP, 0.0, CONCENTRATION_EXT_STARTUP_TIME));
  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_INLET, 0.0, CONCENTRATION_EXT_STARTUP_TIME));

  L2Space<double>space_rho(&mesh_flow, P_INIT_FLOW);
  L2Space<double>space_rho_v_x(&mesh_flow, P_INIT_FLOW);
  L2Space<double>space_rho_v_y(&mesh_flow, P_INIT_FLOW);
  L2Space<double>space_e(&mesh_flow, P_INIT_FLOW);
  // Space<double> for concentration.
  H1Space<double> space_c(&mesh_concentration, &bcs_concentration, P_INIT_CONCENTRATION);

  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  ConstantSolution<double> sln_rho(&mesh_flow, RHO_EXT);
  ConstantSolution<double> sln_rho_v_x(&mesh_flow, RHO_EXT * V1_EXT);
  ConstantSolution<double> sln_rho_v_y(&mesh_flow, RHO_EXT * V2_EXT);
  ConstantSolution<double> sln_e(&mesh_flow, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
  ConstantSolution<double> sln_c(&mesh_concentration, 0.0);

  ConstantSolution<double> prev_rho(&mesh_flow, RHO_EXT);
  ConstantSolution<double> prev_rho_v_x(&mesh_flow, RHO_EXT * V1_EXT);
  ConstantSolution<double> prev_rho_v_y(&mesh_flow, RHO_EXT * V2_EXT);
  ConstantSolution<double> prev_e(&mesh_flow, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
  ConstantSolution<double> prev_c(&mesh_concentration, 0.0);

  Solution<double> rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e, rsln_c;

  // Numerical flux.
  OsherSolomonNumericalFlux num_flux(KAPPA);

  // Initialize weak formulation.
  WeakForm<double>* wf = NULL;
  if(SEMI_IMPLICIT)
    wf = new EulerEquationsWeakFormSemiImplicitCoupled(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM,
    BDY_SOLID_WALL_TOP, BDY_INLET, BDY_OUTLET, BDY_NATURAL_CONCENTRATION, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c, EPSILON, (P_INIT_FLOW == 0 && CAND_LIST_FLOW == H2D_H_ANISO));
  else
    wf = new EulerEquationsWeakFormExplicitCoupled(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM,
    BDY_SOLID_WALL_TOP, BDY_INLET, BDY_OUTLET, BDY_NATURAL_CONCENTRATION, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c, EPSILON, (P_INIT_FLOW == 0 && CAND_LIST_FLOW == H2D_H_ANISO));

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 400));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 400));
  ScalarView s5("Concentration", new WinGeom(700, 400, 600, 400));

  OrderView order_view_flow("Orders - flow", new WinGeom(700, 350, 600, 400));
  OrderView order_view_conc("Orders - concentration", new WinGeom(700, 700, 600, 400));

  // Initialize refinement selector.
  L2ProjBasedSelector<double> l2selector_flow(CAND_LIST_FLOW, CONV_EXP, H2DRS_DEFAULT_ORDER);
  L2ProjBasedSelector<double> l2selector_concentration(CAND_LIST_CONCENTRATION, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Set up Advection-Diffusion-Equation stability calculation class.
  ADEStabilityCalculation ADES(ADVECTION_STABILITY_CONSTANT, DIFFUSION_STABILITY_CONSTANT, EPSILON);

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 100.0; t += time_step)
  {
    time_step_out << time_step << std::endl;
    info("---- Time step %d, time %3.5f.", iteration++, t);

    if(iteration == 2) {
      ERR_STOP_FLOW = 0.55;
      ERR_STOP_CONCENTRATION = 8.3;
    }

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && (REFINEMENT_COUNT_FLOW > 0 || REFINEMENT_COUNT_CONCENTRATION > 0)) {
      info("Global mesh derefinement.");
      if(REFINEMENT_COUNT_FLOW > 0) {
        REFINEMENT_COUNT_FLOW = 0;
        space_rho.unrefine_all_mesh_elements();
        if(CAND_LIST_FLOW == H2D_HP_ANISO)
          space_rho.adjust_element_order(-1, P_INIT_FLOW);
        space_rho_v_x.copy_orders(&space_rho);
        space_rho_v_y.copy_orders(&space_rho);
        space_e.copy_orders(&space_rho);
      }
      if(REFINEMENT_COUNT_CONCENTRATION > 0) {
        REFINEMENT_COUNT_CONCENTRATION = 0;
        space_c.unrefine_all_mesh_elements();
        space_c.adjust_element_order(-1, P_INIT_CONCENTRATION);
      }
    }

    // Adaptivity loop:
    int as = 1; 
    bool done = false;
    do
    {
      info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      int order_increase = 0;
      if(CAND_LIST_FLOW == H2D_HP_ANISO)
        order_increase = 1;
      Hermes::vector<Space<double> *>* ref_spaces = Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e, &space_c), order_increase);
      if(CAND_LIST_FLOW != H2D_HP_ANISO)
        (*ref_spaces)[4]->adjust_element_order(+1, P_INIT_CONCENTRATION);

      Hermes::vector<const Space<double> *> ref_spaces_const((*ref_spaces)[0], (*ref_spaces)[1], 
        (*ref_spaces)[2], (*ref_spaces)[3], (*ref_spaces)[4]);

      // Project the previous time level solution onto the new fine mesh.
      info("Projecting the previous time level solution onto the new fine mesh.");
      OGProjection<double>::project_global(ref_spaces_const, Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c), 
        Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c), matrix_solver);

      /*
      if(iteration == 1) {
        if(CAND_LIST_FLOW == H2D_HP_ANISO)
        {
          prev_rho.set_const((*ref_spaces)[4]->get_mesh(), RHO_EXT);
          prev_rho_v_x.set_const((*ref_spaces)[4]->get_mesh(), RHO_EXT * V1_EXT);
          prev_rho_v_y.set_const((*ref_spaces)[4]->get_mesh(), RHO_EXT * V2_EXT);
          prev_e.set_const((*ref_spaces)[4]->get_mesh(), QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
        }
        prev_c.set_const((*ref_spaces)[4]->get_mesh(), 0.0);
      }
      */

      if(as > 1) {
        delete rsln_rho.get_mesh();
        rsln_rho.own_mesh = false;
        delete rsln_rho_v_x.get_mesh();
        rsln_rho_v_x.own_mesh = false;
        delete rsln_rho_v_y.get_mesh();
        rsln_rho_v_y.own_mesh = false;
        delete rsln_e.get_mesh();
        rsln_e.own_mesh = false;
        delete rsln_c.get_mesh();
        rsln_c.own_mesh = false;
      }

      // Report NDOFs.
      info("ndof_coarse: %d, ndof_fine: %d.", 
        Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e, &space_c)), Space<double>::get_num_dofs(ref_spaces_const));

      // Very imporant, set the meshes for the flow as the same.
      (*ref_spaces)[1]->get_mesh()->set_seq((*ref_spaces)[0]->get_mesh()->get_seq());
      (*ref_spaces)[2]->get_mesh()->set_seq((*ref_spaces)[0]->get_mesh()->get_seq());
      (*ref_spaces)[3]->get_mesh()->set_seq((*ref_spaces)[0]->get_mesh()->get_seq());

      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
      Vector<double>* rhs = create_vector<double>(matrix_solver);
      LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver, matrix, rhs);

      // Initialize the FE problem.
      bool is_linear = true;
      DiscreteProblem<double> dp(wf, ref_spaces_const);
      if(SEMI_IMPLICIT)
        static_cast<EulerEquationsWeakFormSemiImplicitCoupled*>(wf)->set_time_step(time_step);
      else
        static_cast<EulerEquationsWeakFormExplicitCoupled*>(wf)->set_time_step(time_step);

      // Assemble stiffness matrix and rhs.
      info("Assembling the stiffness matrix and right-hand side vector.");
      dp.assemble(matrix, rhs);

      // Solve the matrix problem.
      info("Solving the matrix problem.");
      if (solver->solve())
        Solution<double>::vector_to_solutions(solver->get_sln_vector(), ref_spaces_const, 
        Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e, &rsln_c));
      else
        error ("Matrix solver failed.\n");

      Hermes::vector<const Space<double>*> flow_spaces((*ref_spaces)[0], (*ref_spaces)[1], (*ref_spaces)[2], (*ref_spaces)[3]);

      double* flow_solution_vector = new double[Space<double>::get_num_dofs(flow_spaces)];

      OGProjection<double>::project_global(flow_spaces, Hermes::vector<MeshFunction<double> *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), flow_solution_vector);

			FluxLimiter flux_limiter(FluxLimiter::Krivodonova, flow_solution_vector, flow_spaces);

			flux_limiter.limit_according_to_detector();

			flux_limiter.get_limited_solutions(Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
			
      if(SHOCK_CAPTURING)
        flux_limiter.get_limited_solutions(Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
      
      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection<double>::project_global(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x,
        &space_rho_v_y, &space_e, &space_c), Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e, &rsln_c),
        Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e, &sln_c), matrix_solver,
        Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));

      util_time_step = time_step;
      if(SEMI_IMPLICIT)
        CFL.calculate_semi_implicit(Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), &mesh_flow, util_time_step);
      else
        CFL.calculate(Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), &mesh_flow, util_time_step);

      time_step = util_time_step;

      ADES.calculate(Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y), &mesh_concentration, util_time_step);

      // Calculate element errors and total error estimate.
      info("Calculating error estimates.");
      Adapt<double> adaptivity_flow(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x,
        &space_rho_v_y, &space_e), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));

      double err_est_rel_total_flow = adaptivity_flow.calc_err_est(Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e),
        Hermes::vector<Solution<double>*>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e)) * 100;

      Adapt<double> adaptivity_concentration(&space_c, HERMES_L2_NORM);

      double err_est_rel_total_concentration = adaptivity_concentration.calc_err_est(&sln_c, &rsln_c) * 100;

      // Report results.
      info("Error estimate for the flow part: %g%%", err_est_rel_total_flow);

      info("Error estimate for the concentration part: %g%%", err_est_rel_total_concentration);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total_flow < ERR_STOP_FLOW && err_est_rel_total_concentration < ERR_STOP_CONCENTRATION)
      {
        done = true;
        if(SHOCK_CAPTURING)
          flux_limiter.limit_according_to_detector(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
            &space_rho_v_y, &space_e));
      }

      else
      {
        info("Adapting coarse meshes.");
        if(err_est_rel_total_flow > ERR_STOP_FLOW)
        {
          done = adaptivity_flow.adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&l2selector_flow, &l2selector_flow, &l2selector_flow, &l2selector_flow), 
            THRESHOLD, STRATEGY, MESH_REGULARITY);
          REFINEMENT_COUNT_FLOW++;
        }
        else
          done = true;
        if(err_est_rel_total_concentration > ERR_STOP_CONCENTRATION)
        {
          if(!adaptivity_concentration.adapt(&l2selector_concentration, THRESHOLD, STRATEGY, MESH_REGULARITY))
            done = false;
          REFINEMENT_COUNT_CONCENTRATION++;
        }

        if (Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e, &space_c)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

      // Save orders.
      if((iteration - 1) % EVERY_NTH_STEP == 0 && done)
      {
        if(HERMES_VISUALIZATION)
        {
          Hermes::vector<Space<double> *>* ref_spaces_local = Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_c), 0);
          order_view_flow.show((*ref_spaces_local)[0]);
          order_view_conc.show((*ref_spaces_local)[1]);
          order_view_flow.save_numbered_screenshot("FlowMesh%i.bmp", (int)(iteration / 5), true);
          order_view_conc.save_numbered_screenshot("ConcentrationMesh%i.bmp", (int)(iteration / 5), true);
          for(unsigned int i = 0; i < ref_spaces_local->size(); i++) {
            delete (*ref_spaces_local)[i]->get_mesh();
            delete (*ref_spaces_local)[i];
          }
        }
        if(VTK_VISUALIZATION)
        {
          Orderizer ord;
          char filename[40];
          sprintf(filename, "Flow-mesh-%i.vtk", iteration - 1);
          ord.save_orders_vtk((*ref_spaces)[0], filename);
          sprintf(filename, "Concentration-mesh-%i.vtk", iteration - 1);
          ord.save_orders_vtk((*ref_spaces)[4], filename);
        }
      }

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
    }
    while (done == false);

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&rsln_rho);
    prev_rho_v_x.copy(&rsln_rho_v_x);
    prev_rho_v_y.copy(&rsln_rho_v_y);
    prev_e.copy(&rsln_e);
    prev_c.copy(&rsln_c);

    delete rsln_rho.get_mesh();
    rsln_rho.own_mesh = false;
    delete rsln_rho_v_x.get_mesh();
    rsln_rho_v_x.own_mesh = false;
    delete rsln_rho_v_y.get_mesh();
    rsln_rho_v_y.own_mesh = false;
    delete rsln_e.get_mesh();
    rsln_e.own_mesh = false;
    delete rsln_c.get_mesh();
    rsln_c.own_mesh = false;

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) {
      // Hermes visualization.
      if(HERMES_VISUALIZATION)
      {
        Mach_number.reinit();
        pressure.reinit();

        pressure_view.show_mesh(false);
        pressure_view.show(&pressure);
        pressure_view.set_scale_format("%1.3f");

        Mach_number_view.show_mesh(false);
        Mach_number_view.set_scale_format("%1.3f");
        Mach_number_view.show(&Mach_number);

        s5.show_mesh(false);
        s5.set_scale_format("%0.3f");
        s5.show(&prev_c);

        pressure_view.save_numbered_screenshot("pressure%i.bmp", (int)(iteration / 5), true);
        Mach_number_view.save_numbered_screenshot("Mach_number%i.bmp", (int)(iteration / 5), true);
        s5.save_numbered_screenshot("concentration%i.bmp", (int)(iteration / 5), true);

      }
      // Output solution in VTK format.
      if(VTK_VISUALIZATION)
      {
        pressure.reinit();
        Mach_number.reinit();
        Linearizer lin_pressure;
        char filename[40];
        sprintf(filename, "pressure-3D-%i.vtk", iteration - 1);
        lin_pressure.save_solution_vtk(&pressure, filename, "Pressure", true);
        Linearizer lin_mach;
        sprintf(filename, "Mach number-3D-%i.vtk", iteration - 1);
        lin_mach.save_solution_vtk(&Mach_number, filename, "MachNumber", true);
        Linearizer lin_concentration;
        sprintf(filename, "Concentration-%i.vtk", iteration - 1);
        lin_concentration.save_solution_vtk(&prev_c, filename, "Concentration", true);

      }
    }
  }

  /*
  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();
  s5.close();
  */

  return 0;
}
