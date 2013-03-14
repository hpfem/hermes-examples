#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// This example solves the compressible Euler equations using a basic
// piecewise-constant finite volume method, or Discontinuous Galerkin 
// method of higher order with no adaptivity.
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
shockCapturingType SHOCK_CAPTURING_TYPE = KUZMIN;
// Quantitative parameter of the discontinuity detector in case of Krivodonova.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;
// Quantitative parameter of the shock capturing in case of Feistauer.
const double NU_1 = 0.1;
const double NU_2 = 0.1;

// Initial polynomial degree.      
const int P_INIT = 1;                                                   
// Number of initial uniform mesh refinements.      
const int INIT_REF_NUM = 3;                                              
// CFL value.
double CFL_NUMBER = 0.1;                                
// Initial time step.
double time_step = 1E-4;                                
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
const MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

double KAPPA = 1.4;

// Equation parameters.
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
  // Load the mesh->
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
  SpaceSharedPtr<double> space_stabilization(new L2Space<double>(mesh, 0));
  int ndof = Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double>  >(space_rho, space_rho_v_x, space_rho_v_y, space_e));
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  MeshFunctionSharedPtr<double> prev_rho(new ConstantSolution<double>(mesh, RHO_INIT));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double> (mesh, RHO_INIT * V1_INIT));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double> (mesh, RHO_INIT * V2_INIT));
  MeshFunctionSharedPtr<double> prev_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_INIT, RHO_INIT * V1_INIT, RHO_INIT * V2_INIT, PRESSURE_INIT, KAPPA)));

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
  MeshFunctionSharedPtr<double> Mach_number(new MachNumberFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), KAPPA));
  MeshFunctionSharedPtr<double> pressure(new PressureFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), KAPPA));
 
  ScalarView pressure_view("Pressure", new WinGeom(700, 400, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 700, 600, 300));
  ScalarView s1("prev_rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("prev_rho_v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("prev_rho_v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 400, 600, 300));

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
    wf.set_stabilization(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, NU_1, NU_2);

  // Initialize the FE problem.
  DiscreteProblemLinear<double> dp(&wf, Hermes::vector<SpaceSharedPtr<double>  >(space_rho, space_rho_v_x, space_rho_v_y, space_e));
  LinearSolver<double> solver(&dp);
  solver.output_matrix();
  solver.output_rhs();

  int iteration = 0.;
  double t = 0.0;

  // Time stepping loop.
  for(; t < 6.0; t += time_step)
  {
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

    // Set the current time step.
    wf.set_current_time_step(time_step);

    // Assemble the stiffness matrix and rhs.
    Hermes::Mixins::Loggable::Static::info("Assembling the stiffness matrix and right-hand side vector.");

    // Solve the matrix problem.
    Hermes::Mixins::Loggable::Static::info("Solving the matrix problem.");
    solver.solve();
    if(!SHOCK_CAPTURING || SHOCK_CAPTURING_TYPE == FEISTAUER)
    {
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
        space_rho_v_y, space_e), Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e));
    }
    else
    {      
      FluxLimiter* flux_limiter;
      if(SHOCK_CAPTURING_TYPE == KUZMIN)
        flux_limiter = new FluxLimiter(FluxLimiter::Kuzmin, solver.get_sln_vector(), Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
        space_rho_v_y, space_e), true);
      else
        flux_limiter = new FluxLimiter(FluxLimiter::Krivodonova, solver.get_sln_vector(), Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
        space_rho_v_y, space_e), true);

      if(SHOCK_CAPTURING_TYPE == KUZMIN)
        flux_limiter->limit_second_orders_according_to_detector();

      flux_limiter->limit_according_to_detector();

      flux_limiter->get_limited_solutions(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e));
    }

    CFL.calculate_semi_implicit(Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), mesh, time_step);

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) 
    {
      // Hermes visualization.
      if(HERMES_VISUALIZATION) 
      {
        Mach_number->reinit();
        pressure->reinit();
        pressure_view.show(pressure);
        Mach_number_view.show(Mach_number);
      }
      // Output solution in VTK format.
      if(VTK_VISUALIZATION) 
      {
        pressure->reinit();
        Mach_number->reinit();
        Linearizer lin_pressure;
        char filename[40];
        sprintf(filename, "pressure-3D-%i.vtk", iteration - 1);
        lin_pressure.save_solution_vtk(pressure, filename, "Pressure", true);
        Linearizer lin_mach;
        sprintf(filename, "Mach number-3D-%i.vtk", iteration - 1);
        lin_mach.save_solution_vtk(Mach_number, filename, "MachNumber", true);
      }
    }
  }

  pressure_view.close();
  Mach_number_view.close();

  return 0;
}
