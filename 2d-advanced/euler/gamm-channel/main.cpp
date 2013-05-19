#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// This example solves the compressible Euler equations using a basic
// Discontinuous Galerkin method of higher order with no adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: GAMM channel, see mesh file GAMM-channel.mesh
//
// BC: Solid walls, inlet, outlet.
//
// IC: Constant state identical to inlet.
//

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
bool SHOCK_CAPTURING = false;
shockCapturingType SHOCK_CAPTURING_TYPE = FEISTAUER;
// Quantitative parameter of the discontinuity detector in case of Krivodonova.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;
// Quantitative parameter of the shock capturing in case of Feistauer.
const double NU_1 = 0.1;
const double NU_2 = 0.1;

// Initial polynomial degree.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.    
const int INIT_REF_NUM = 4;
// CFL value.
double CFL_NUMBER = 0.9;
// Initial time step.
double time_step = 1E-6;

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
std::string BDY_INLET = "1";
std::string BDY_OUTLET = "2";
std::string BDY_SOLID_WALL_BOTTOM = "3";
std::string BDY_SOLID_WALL_TOP = "4";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  #pragma region 1. Load mesh and initialize spaces.
    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    MeshReaderH2D mloader;
    mloader.load("GAMM-channel.mesh", mesh);

    // Perform initial mesh refinements.
    for (int i = 0; i < INIT_REF_NUM; i++) 
      mesh->refine_all_elements(0, true);

    // Initialize boundary condition types and spaces with default shapesets.
    SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT));
    SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT));
    SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT));
    SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT));
    Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);
    int ndof = Space<double>::get_num_dofs(spaces);
    Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);
  #pragma endregion

  #pragma region 2. Initialize solutions, set initial conditions.
    MeshFunctionSharedPtr<double> prev_rho(new ConstantSolution<double>(mesh, RHO_EXT));
    MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
    MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
    MeshFunctionSharedPtr<double> prev_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
    Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  #pragma endregion

  #pragma region 3. Filters for visualization of Mach number, pressure + visualization setup.
    MeshFunctionSharedPtr<double> Mach_number(new MachNumberFilter (prev_slns, KAPPA));
    MeshFunctionSharedPtr<double> pressure(new PressureFilter(prev_slns, KAPPA));

    ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
    pressure_view.show_contours(.1);
    pressure_view.show_mesh(false);
    ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
    Mach_number_view.show_contours(.02);
    Mach_number_view.show_mesh(false);
    VectorView velocity_view("Velocity", new WinGeom(0, 330, 600, 300));
  #pragma endregion

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  #pragma region 4. Initialize weak formulation -> EF problem -> linear solver.
    Hermes::vector<std::string> solid_wall_markers(BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP);
    Hermes::vector<std::string> inlet_markers;
    inlet_markers.push_back(BDY_INLET);
    Hermes::vector<std::string> outlet_markers;
    outlet_markers.push_back(BDY_OUTLET);

    EulerEquationsWeakFormSemiImplicit wf(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT,solid_wall_markers,
      inlet_markers, outlet_markers, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, (P_INIT == 0));

    // Initialize the FE problem.
    Space<double>::assign_dofs(spaces);

    LinearSolver<double> solver(&wf, spaces);
  #pragma endregion
  
  #pragma region 5.Time stepping loop.
    int iteration = 0;
    for(double t = 0.0; t < 10.0; t += time_step)
    {
      // Info.
      Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

      // Set the current time step.
      wf.set_current_time_step(time_step);
    
      try
      {
        // Solve.
        solver.solve();

        #pragma region *. Get the solution with optional shock capturing.
        if(!SHOCK_CAPTURING)
          Solution<double>::vector_to_solutions(solver.get_sln_vector(), spaces, prev_slns);
        else
        {
          FluxLimiter* flux_limiter;
          if(SHOCK_CAPTURING_TYPE == KUZMIN)
          flux_limiter = new FluxLimiter(FluxLimiter::Kuzmin, solver.get_sln_vector(), spaces, true);
          else
          flux_limiter = new FluxLimiter(FluxLimiter::Krivodonova, solver.get_sln_vector(), spaces);

          if(SHOCK_CAPTURING_TYPE == KUZMIN)
            flux_limiter->limit_second_orders_according_to_detector();

          flux_limiter->limit_according_to_detector();

          flux_limiter->get_limited_solutions(prev_slns);
        }
        #pragma endregion
      }
      catch(std::exception& e) { std::cout << e.what(); }

      // Calculate time step according to CFL condition.
      CFL.calculate(prev_slns, mesh, time_step);

      #pragma region 5.1. Visualization
        if((iteration - 1) % EVERY_NTH_STEP == 0) 
        {
          // Hermes visualization.
          if(HERMES_VISUALIZATION) 
          {
            Mach_number->reinit();
            pressure->reinit();
            pressure_view.show(pressure);
            Mach_number_view.show(Mach_number);
            velocity_view.show(prev_rho_v_x, prev_rho_v_y);
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
        #pragma endregion
    }
  #pragma endregion

  // Done.
  return 0;
}
