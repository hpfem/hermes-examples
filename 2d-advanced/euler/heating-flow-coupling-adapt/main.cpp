#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "../coupling.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::RefinementSelectors;
using namespace Hermes::Hermes2D::Views;

// Visualization.
// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;
// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;


// Initial polynomial degree.
const int P_INIT_FLOW = 0;
const int P_INIT_HEAT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;

// Shock capturing.
enum shockCapturingType
{
  KUZMIN,
  KRIVODONOVA
};
bool SHOCK_CAPTURING = true;
shockCapturingType SHOCK_CAPTURING_TYPE = KUZMIN;
// Quantitative parameter of the discontinuity detector in case of Krivodonova.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 1.5;
// Initial pressure (dimensionless).
const double P_INITIAL_HIGH = 1.5;  
// Initial pressure (dimensionless).
const double P_INITIAL_LOW = 1.0;   
// Inlet density (dimensionless).   
const double RHO_EXT = 1.0;         
// Initial density (dimensionless).   
const double RHO_INITIAL_HIGH = 1.0;
// Initial density (dimensionless).   
const double RHO_INITIAL_LOW = 0.66;
// Inlet x-velocity (dimensionless).
const double V1_EXT = 0.0;          
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;          
// Kappa.
const double KAPPA = 1.4;
// Lambda.
const double LAMBDA = 1e3;
// heat_capacity.
const double C_P = 1e2;
// heat flux through the inlet.
const double HEAT_FLUX = 1e-3;

// CFL value.
const double CFL_NUMBER = 0.1;                               
// Initial time step.
double time_step = 1E-5;
// Stability for the concentration part.
double ADVECTION_STABILITY_CONSTANT = 1e16;
const double DIFFUSION_STABILITY_CONSTANT = 1e16;

// Boundary markers.
const std::string BDY_INLET = "Inlet";
const std::string BDY_SOLID_WALL = "Solid";

// Area (square) size.
// Must be in accordance with the mesh file.
const double MESH_SIZE = 3.0;

// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 5;

// Adaptivity
int REFINEMENT_COUNT = 0;                         
const double THRESHOLD = 0.3;                     
const int STRATEGY = 0;                           
CandList CAND_LIST = H2D_HP_ANISO;                
const int MAX_P_ORDER = -1;                       
const int MESH_REGULARITY = -1;                   
const double CONV_EXP = 1;                        
double ERR_STOP = 5.0;

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh->
  MeshSharedPtr mesh(new Mesh), mesh_heat(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements(0, true);

  mesh_heat->copy(mesh);
  //mesh_heat->refine_towards_boundary("Inlet", 3, false);
  //mesh->refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  Hermes2D::DefaultEssentialBCConst<double> bc_temp_zero("Solid", 0.0);
  Hermes2D::DefaultEssentialBCConst<double> bc_temp_nonzero("Inlet", 1.0);
  Hermes::vector<Hermes2D::EssentialBoundaryCondition<double>*> bc_vector(&bc_temp_zero, &bc_temp_nonzero);
  EssentialBCs<double> bcs(bc_vector);

  SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT_FLOW));
  SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT_FLOW));
  SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT_FLOW));
  SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT_FLOW));
  SpaceSharedPtr<double> space_temp(new H1Space<double>(mesh_heat, &bcs, P_INIT_HEAT));

  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e, space_temp);
  Hermes::vector<SpaceSharedPtr<double>  > cspaces(space_rho, space_rho_v_x, space_rho_v_y, space_e, space_temp);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double>  >(space_rho, space_rho_v_x, space_rho_v_y, space_e, space_temp));
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  MeshFunctionSharedPtr<double> prev_rho(new InitialSolutionLinearProgress(mesh, RHO_INITIAL_HIGH, RHO_INITIAL_LOW, MESH_SIZE));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double> (mesh, 0.0));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double> (mesh, 0.0));
  MeshFunctionSharedPtr<double> prev_e(new InitialSolutionLinearProgress (mesh, QuantityCalculator::calc_energy(RHO_INITIAL_HIGH, RHO_INITIAL_HIGH * V1_EXT, RHO_INITIAL_HIGH * V2_EXT, P_INITIAL_HIGH, KAPPA), QuantityCalculator::calc_energy(RHO_INITIAL_LOW, RHO_INITIAL_LOW * V1_EXT, RHO_INITIAL_LOW * V2_EXT, P_INITIAL_LOW, KAPPA), MESH_SIZE));
  MeshFunctionSharedPtr<double> prev_temp(new ConstantSolution<double> (mesh_heat, 0.0));
  Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, prev_temp);
  Hermes::vector<const MeshFunctionSharedPtr<double> > cprev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, prev_temp);

  MeshFunctionSharedPtr<double> sln_rho(new InitialSolutionLinearProgress(mesh, RHO_INITIAL_HIGH, RHO_INITIAL_LOW, MESH_SIZE));
  MeshFunctionSharedPtr<double> sln_rho_v_x(new ConstantSolution<double> (mesh, 0.0));
  MeshFunctionSharedPtr<double> sln_rho_v_y(new ConstantSolution<double> (mesh, 0.0));
  MeshFunctionSharedPtr<double> sln_e(new InitialSolutionLinearProgress (mesh, QuantityCalculator::calc_energy(RHO_INITIAL_HIGH, RHO_INITIAL_HIGH * V1_EXT, RHO_INITIAL_HIGH * V2_EXT, P_INITIAL_HIGH, KAPPA), QuantityCalculator::calc_energy(RHO_INITIAL_LOW, RHO_INITIAL_LOW * V1_EXT, RHO_INITIAL_LOW * V2_EXT, P_INITIAL_LOW, KAPPA), MESH_SIZE));
  MeshFunctionSharedPtr<double> sln_temp(new ConstantSolution<double> (mesh_heat, 0.0));

  Hermes::vector<MeshFunctionSharedPtr<double> > slns(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e, sln_temp);
  Hermes::vector<const MeshFunctionSharedPtr<double> > cslns(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e, sln_temp);

   MeshFunctionSharedPtr<double> rsln_rho(new InitialSolutionLinearProgress(mesh, RHO_INITIAL_HIGH, RHO_INITIAL_LOW, MESH_SIZE));
  MeshFunctionSharedPtr<double> rsln_rho_v_x(new ConstantSolution<double> (mesh, 0.0));
  MeshFunctionSharedPtr<double> rsln_rho_v_y(new ConstantSolution<double> (mesh, 0.0));
  MeshFunctionSharedPtr<double> rsln_e(new InitialSolutionLinearProgress (mesh, QuantityCalculator::calc_energy(RHO_INITIAL_HIGH, RHO_INITIAL_HIGH * V1_EXT, RHO_INITIAL_HIGH * V2_EXT, P_INITIAL_HIGH, KAPPA), QuantityCalculator::calc_energy(RHO_INITIAL_LOW, RHO_INITIAL_LOW * V1_EXT, RHO_INITIAL_LOW * V2_EXT, P_INITIAL_LOW, KAPPA), MESH_SIZE));
  MeshFunctionSharedPtr<double> rsln_temp(new ConstantSolution<double> (mesh_heat, 0.0));
  Hermes::vector<MeshFunctionSharedPtr<double> > rslns(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e, rsln_temp);
  Hermes::vector<const MeshFunctionSharedPtr<double> > crslns(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e, rsln_temp);

  Hermes::vector<SpaceSharedPtr<double> > flow_spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);
  Hermes::vector<MeshFunctionSharedPtr<double> > flow_slns(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e);
  Hermes::vector<MeshFunctionSharedPtr<double> > rflow_slns(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e);

  // Filters for visualization of Mach number, pressure and entropy.
  MeshFunctionSharedPtr<double>  pressure(new PressureFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), KAPPA));
  MeshFunctionSharedPtr<double>  vel_x(new VelocityFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x)));
  MeshFunctionSharedPtr<double>  vel_y(new VelocityFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_y)));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 400, 300));
  VectorView velocity_view("Velocity", new WinGeom(0, 400, 400, 300));
  ScalarView density_view("Density", new WinGeom(500, 0, 400, 300));
  ScalarView temperature_view("Temperature", new WinGeom(500, 400, 400, 300));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>();
  Vector<double>* rhs = create_vector<double>();
  Vector<double>* rhs_stabilization = create_vector<double>();
  LinearMatrixSolver<double>* solver = create_linear_solver<double>( matrix, rhs);

  // Set up stability calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);
  ADEStabilityCalculation ADES(ADVECTION_STABILITY_CONSTANT, DIFFUSION_STABILITY_CONSTANT, LAMBDA);

  // Initialize refinement selector.
  H1ProjBasedSelector<double> l2_selector(CAND_LIST, CONV_EXP, MAX_P_ORDER);
  H1ProjBasedSelector<double> h1_selector(CAND_LIST, CONV_EXP, MAX_P_ORDER);

  // Look for a saved solution on the disk.
  int iteration = 0; double t = 0;

  // Initialize weak formulation.
  Hermes::vector<std::string> solid_wall_markers;
  solid_wall_markers.push_back(BDY_SOLID_WALL);
  Hermes::vector<std::string> inlet_markers;
  inlet_markers.push_back(BDY_INLET);
  Hermes::vector<std::string> outlet_markers;

  EulerEquationsWeakFormSemiImplicitCoupledWithHeat wf(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, solid_wall_markers, 
    inlet_markers, outlet_markers, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, prev_temp, LAMBDA, C_P);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, cspaces);

  // Time stepping loop.
  for(; t < 10.0; t += time_step)
  {
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0) 
    {
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      REFINEMENT_COUNT = 0;

      space_rho->unrefine_all_mesh_elements(true);
      space_temp->unrefine_all_mesh_elements(true);

      space_rho->adjust_element_order(-1, P_INIT_FLOW);
      space_rho_v_x->adjust_element_order(-1, P_INIT_FLOW);
      space_rho_v_y->adjust_element_order(-1, P_INIT_FLOW);
      space_e->adjust_element_order(-1, P_INIT_FLOW);
      space_temp->adjust_element_order(-1, P_INIT_HEAT);
      Space<double>::assign_dofs(spaces);
    }

    // Set the current time step.
    wf.set_current_time_step(time_step);

    // Adaptivity loop:
    int as = 1;
    bool done = false;
    do
    {
      Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space->
      int order_increase = 1;

      Mesh::ReferenceMeshCreator refMeshCreatorFlow(mesh);
      MeshSharedPtr ref_mesh_flow = refMeshCreatorFlow.create_ref_mesh();
      Mesh::ReferenceMeshCreator refMeshCreatorTemperature(mesh_heat);
      MeshSharedPtr ref_mesh_heat = refMeshCreatorTemperature.create_ref_mesh();

      Space<double>::ReferenceSpaceCreator refSpaceCreatorRho(space_rho, ref_mesh_flow, order_increase);
      SpaceSharedPtr<double> ref_space_rho = refSpaceCreatorRho.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVx(space_rho_v_x, ref_mesh_flow, order_increase);
      SpaceSharedPtr<double> ref_space_rho_v_x = refSpaceCreatorRhoVx.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVy(space_rho_v_y, ref_mesh_flow, order_increase);
      SpaceSharedPtr<double> ref_space_rho_v_y = refSpaceCreatorRhoVy.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorE(space_e, ref_mesh_flow, order_increase);
      SpaceSharedPtr<double> ref_space_e = refSpaceCreatorE.create_ref_space();

      Space<double>::ReferenceSpaceCreator refSpaceCreatorT(space_temp, ref_mesh_heat, order_increase);
      SpaceSharedPtr<double> ref_space_temp = refSpaceCreatorT.create_ref_space();

      Hermes::vector<SpaceSharedPtr<double> > ref_spaces(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e, ref_space_temp);
      Hermes::vector<SpaceSharedPtr<double>  > cref_spaces(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e, ref_space_temp);

      Hermes::vector<SpaceSharedPtr<double>  > cref_flow_spaces(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e);

      // Project the previous time level solution onto the new fine mesh->
      Hermes::Mixins::Loggable::Static::info("Projecting the previous time level solution onto the new fine mesh->");
      OGProjection<double> ogProjection; ogProjection.project_global(cref_spaces, prev_slns, prev_slns, Hermes::vector<Hermes::Hermes2D::NormType>(), iteration > 1);

      // Report NDOFs.
      Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d.", 
        Space<double>::get_num_dofs(spaces), Space<double>::get_num_dofs(cref_spaces));

      // Assemble the reference problem.
      Hermes::Mixins::Loggable::Static::info("Solving on reference mesh->");
      DiscreteProblem<double> dp(&wf, cref_spaces);

      // Assemble the stiffness matrix and rhs.
      Hermes::Mixins::Loggable::Static::info("Assembling the stiffness matrix and right-hand side vector.");
      dp.assemble(matrix, rhs);

      // Solve the matrix problem.
      Hermes::Mixins::Loggable::Static::info("Solving the matrix problem.");
      if(solver->solve())
      {
        if(!SHOCK_CAPTURING)
        {
          Solution<double>::vector_to_solutions(solver->get_sln_vector(), cref_spaces, rslns);
        }
        else
        {
          FluxLimiter* flux_limiter;
          if(SHOCK_CAPTURING_TYPE == KUZMIN)
            flux_limiter = new FluxLimiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), cref_flow_spaces, true);
          else
            flux_limiter = new FluxLimiter(FluxLimiter::Krivodonova, solver->get_sln_vector(), cref_flow_spaces);

          if(SHOCK_CAPTURING_TYPE == KUZMIN)
            flux_limiter->limit_second_orders_according_to_detector(flow_spaces);

          flux_limiter->limit_according_to_detector(flow_spaces);

          flux_limiter->get_limited_solutions(rflow_slns);

          Solution<double>::vector_to_solution(solver->get_sln_vector() + Space<double>::get_num_dofs(cref_flow_spaces), ref_space_temp, rsln_temp);
        }
      }
      else
        throw Hermes::Exceptions::Exception("Matrix solver failed.\n");

      // Project the fine mesh solution onto the coarse mesh->
      Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh->");
      ogProjection.project_global(cspaces, rslns, slns, Hermes::vector<NormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_H1_NORM)); 

      // Calculate element errors and total error estimate.
      Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");

      Adapt<double>* adaptivity = new Adapt<double>(spaces, Hermes::vector<NormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_H1_NORM));
      adaptivity->set_error_form(new CouplingErrorFormVelocity(velX, C_P));
      adaptivity->set_error_form(new CouplingErrorFormVelocity(velY, C_P));
      adaptivity->set_norm_form(new CouplingErrorFormVelocity(velX, C_P));
      adaptivity->set_norm_form(new CouplingErrorFormVelocity(velY, C_P));

      adaptivity->set_error_form(new CouplingErrorFormTemperature(velX, C_P));
      adaptivity->set_error_form(new CouplingErrorFormTemperature(velY, C_P));
      adaptivity->set_norm_form(new CouplingErrorFormTemperature(velX, C_P));
      adaptivity->set_norm_form(new CouplingErrorFormTemperature(velY, C_P));
      Hermes::vector<double> component_errors;
      double error_value = adaptivity->calc_err_est(slns, rslns, &component_errors) * 100;
      
      std::cout << std::endl;
      for(int k = 0; k < 5; k ++)
        std::cout << k << ':' << component_errors[k] << std::endl;
      std::cout << std::endl;

      CFL.calculate_semi_implicit(rflow_slns, ref_mesh_flow, time_step);

      double util_time_step = time_step;

      ADES.calculate(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y), ref_mesh_heat, util_time_step);

      if(util_time_step < time_step)
        time_step = util_time_step;

      // Report results.
      Hermes::Mixins::Loggable::Static::info("Error: %g%%.", error_value);

      // If err_est too large, adapt the mesh->
      if (error_value < ERR_STOP)
        done = true;
      else
      {
        Hermes::Mixins::Loggable::Static::info("Adapting coarse space->");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&l2_selector, &l2_selector, &l2_selector, &l2_selector, &h1_selector), 
          THRESHOLD, STRATEGY, MESH_REGULARITY);
      }

      as++;

      // Visualization.
      if((iteration - 1) % EVERY_NTH_STEP == 0) 
      {
        // Hermes visualization.
        if(HERMES_VISUALIZATION) 
        {
          pressure->reinit();
          vel_x->reinit();
          vel_y->reinit();
          pressure_view.show(pressure);
          velocity_view.show(vel_x, vel_y);
          density_view.show(prev_rho);
          temperature_view.show(prev_temp, HERMES_EPS_HIGH);
        }
        // Output solution in VTK format.
        if(VTK_VISUALIZATION) 
        {
          pressure->reinit();
          Linearizer lin_pressure;
          char filename[40];
          sprintf(filename, "pressure-3D-%i.vtk", iteration - 1);
          lin_pressure.save_solution_vtk(pressure, filename, "Pressure", true);
        }
      }

      delete adaptivity;
    }
    while (!done);

    prev_rho->copy(rsln_rho);
    prev_rho_v_x->copy(rsln_rho_v_x);
    prev_rho_v_y->copy(rsln_rho_v_y);
    prev_e->copy(rsln_e);
    prev_temp->copy(rsln_temp);
  }

  pressure_view.close();
  velocity_view.close();
  density_view.close();
  temperature_view.close();

  return 0;
}
