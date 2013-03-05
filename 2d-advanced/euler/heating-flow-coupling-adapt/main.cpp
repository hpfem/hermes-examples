#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::RefinementSelectors;
using namespace Hermes::Hermes2D::Views;

// Visualization.
// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = false;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;
// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;


// Initial polynomial degree.
const int P_INIT_FLOW = 1;
const int P_INIT_HEAT = 1;
// Number of initial uniform mesh refinements.    
const int INIT_REF_NUM = 2;

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

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  Hermes2DApi.set_integral_param_value(numThreads, 1);

  // Load the mesh.
  Mesh mesh, mesh_heat;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh.refine_all_elements(0, true);

  mesh_heat.copy(&mesh);
  mesh_heat.refine_towards_boundary("Inlet", 3, false);
  mesh.refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  Hermes2D::DefaultEssentialBCConst<double> bc_temp_zero("Solid", 0.0);
  Hermes2D::DefaultEssentialBCConst<double> bc_temp_nonzero("Inlet", 1.0);
  Hermes::vector<Hermes2D::EssentialBoundaryCondition<double>*> bc_vector(&bc_temp_zero, &bc_temp_nonzero);
  EssentialBCs<double> bcs(bc_vector);

  L2Space<double> space_rho(&mesh, P_INIT_FLOW);
  L2Space<double> space_rho_v_x(&mesh, P_INIT_FLOW);
  L2Space<double> space_rho_v_y(&mesh, P_INIT_FLOW);
  L2Space<double> space_e(&mesh, P_INIT_FLOW);
  H1Space<double> space_temp(&mesh_heat, &bcs, P_INIT_HEAT);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_temp));
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  InitialSolutionLinearProgress prev_rho(&mesh, RHO_INITIAL_HIGH, RHO_INITIAL_LOW, MESH_SIZE);
  ConstantSolution<double> prev_rho_v_x(&mesh, 0.0);
  ConstantSolution<double> prev_rho_v_y(&mesh, 0.0);
  InitialSolutionLinearProgress prev_e(&mesh, QuantityCalculator::calc_energy(RHO_INITIAL_HIGH, RHO_INITIAL_HIGH * V1_EXT, RHO_INITIAL_HIGH * V2_EXT, P_INITIAL_HIGH, KAPPA), QuantityCalculator::calc_energy(RHO_INITIAL_LOW, RHO_INITIAL_LOW * V1_EXT, RHO_INITIAL_LOW * V2_EXT, P_INITIAL_LOW, KAPPA), MESH_SIZE);
  ConstantSolution<double> prev_temp(&mesh_heat, 0.0);

  // Filters for visualization of Mach number, pressure and entropy.
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  VelocityFilter vel_x(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x));
  VelocityFilter vel_y(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_y));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 800, 600));
  VectorView velocity_view("Velocity", new WinGeom(0, 700, 800, 600));
  ScalarView density_view("Density", new WinGeom(900, 0, 800, 600));
  ScalarView temperature_view("Temperature", new WinGeom(900, 700, 800, 600));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>();
  Vector<double>* rhs = create_vector<double>();
  Vector<double>* rhs_stabilization = create_vector<double>();
  LinearMatrixSolver<double>* solver = create_linear_solver<double>( matrix, rhs);

  // Set up stability calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);
  ADEStabilityCalculation ADES(ADVECTION_STABILITY_CONSTANT, DIFFUSION_STABILITY_CONSTANT, LAMBDA);

  // Look for a saved solution on the disk.
  int iteration = 0; double t = 0;

  // Initialize weak formulation.
  Hermes::vector<std::string> solid_wall_markers;
  solid_wall_markers.push_back(BDY_SOLID_WALL);
  Hermes::vector<std::string> inlet_markers;
  inlet_markers.push_back(BDY_INLET);
  Hermes::vector<std::string> outlet_markers;

  EulerEquationsWeakFormSemiImplicitCoupledWithHeat wf(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, solid_wall_markers, 
    inlet_markers, outlet_markers, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_temp, LAMBDA, C_P, HEAT_FLUX);
  
  // Initialize the FE problem.
  DiscreteProblemLinear<double> dp(&wf, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_temp));

  // Time stepping loop.
  for(; t < 10.0; t += time_step)
  {
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

    // Set the current time step.
    wf.set_current_time_step(time_step);

    // Assemble the stiffness matrix and rhs.
    Hermes::Mixins::Loggable::Static::info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the matrix problem.
    Hermes::Mixins::Loggable::Static::info("Solving the matrix problem.");
    if(solver->solve())
    {
      if(!SHOCK_CAPTURING)
      {
        Solution<double>::vector_to_solutions(solver->get_sln_vector(), Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e, &space_temp), Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_temp));
      }
      else
      {
        FluxLimiter* flux_limiter;
        if(SHOCK_CAPTURING_TYPE == KUZMIN)
          flux_limiter = new FluxLimiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e), true);
        else
          flux_limiter = new FluxLimiter(FluxLimiter::Krivodonova, solver->get_sln_vector(), Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e));

        if(SHOCK_CAPTURING_TYPE == KUZMIN)
          flux_limiter->limit_second_orders_according_to_detector();

        flux_limiter->limit_according_to_detector();

        flux_limiter->get_limited_solutions(Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));

        Solution<double>::vector_to_solution(solver->get_sln_vector() + Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e)), &space_temp, &prev_temp);
      }
    }
    else
      throw Hermes::Exceptions::Exception("Matrix solver failed.\n");

    CFL.calculate_semi_implicit(Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), &mesh, time_step);

    double util_time_step = time_step;

    ADES.calculate(Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y), &mesh, util_time_step);

    if(util_time_step < time_step)
      time_step = util_time_step;

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) 
    {
      // Hermes visualization.
      if(HERMES_VISUALIZATION) 
      {
        pressure.reinit();
        vel_x.reinit();
        vel_y.reinit();
        pressure_view.show(&pressure);
        velocity_view.show(&vel_x, &vel_y);
        density_view.show(&prev_rho);
        temperature_view.show(&prev_temp, HERMES_EPS_HIGH);
      }
      // Output solution in VTK format.
      if(VTK_VISUALIZATION) 
      {
        pressure.reinit();
        Linearizer lin_pressure;
        char filename[40];
        sprintf(filename, "pressure-3D-%i.vtk", iteration - 1);
        lin_pressure.save_solution_vtk(&pressure, filename, "Pressure", true);
      }
    }
  }

  pressure_view.close();
  velocity_view.close();
  density_view.close();
  temperature_view.close();
  
  return 0;
}
