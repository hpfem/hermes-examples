#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// This example solves the compressible Euler equations using a basic
// piecewise-constant finite volume method, or Discontinuous Galerkin method of higher order with no adaptivity.
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
enum shockCapturingType
{
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

// For saving/loading of solution.
bool REUSE_SOLUTION = false;

// Initial polynomial degree.   
const int P_INIT = 0;                                                      
// Number of initial uniform mesh refinements.          
const int INIT_REF_NUM = 0;                                          
// Number of initial localized mesh refinements.   
const int INIT_REF_NUM_STEP = 2;                                            
// CFL value.
double CFL_NUMBER = 0.25;                                
// Initial time step.
double time_step = 1E-6;                                
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
const MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 1.0;         
// Inlet density (dimensionless).   
const double RHO_EXT = 1.4;       
// Inlet x-velocity (dimensionless).
const double V1_EXT = 3.0;        
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;        
// Kappa.
const double KAPPA = 1.4;         

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
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("ffs.mesh", &mesh);
  mesh.refine_by_criterion(refinement_criterion, INIT_REF_NUM_STEP);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) 
    mesh.refine_all_elements(0, true);

  MeshView m;
  m.show(&mesh);

  // Initialize boundary condition types and spaces with default shapesets.
  L2Space<double> space_rho(&mesh, P_INIT);
  L2Space<double> space_rho_v_x(&mesh, P_INIT);
  L2Space<double> space_rho_v_y(&mesh, P_INIT);
  L2Space<double> space_e(&mesh, P_INIT);
  L2Space<double> space_stabilization(&mesh, 0);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  ConstantSolution<double> prev_rho(&mesh, RHO_EXT);
  ConstantSolution<double> prev_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double> prev_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double> prev_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  ScalarView s1("prev_rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("prev_rho_v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("prev_rho_v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 400, 600, 300));

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Look for a saved solution on the disk.
  CalculationContinuity<double> continuity(CalculationContinuity<double>::onlyTime);
  int iteration = 0; double t = 0;

  if(REUSE_SOLUTION && continuity.have_record_available())
  {
    continuity.get_last_record()->load_mesh(&mesh);
    Hermes::vector<Space<double> *> spaceVector = continuity.get_last_record()->load_spaces(Hermes::vector<Mesh *>(&mesh, &mesh, &mesh, &mesh));
    space_rho.copy(spaceVector[0], &mesh);
    space_rho_v_x.copy(spaceVector[1], &mesh);
    space_rho_v_y.copy(spaceVector[2], &mesh);
    space_e.copy(spaceVector[3], &mesh);
    continuity.get_last_record()->load_solutions(Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e));
    continuity.get_last_record()->load_time_step_length(time_step);
    t = continuity.get_last_record()->get_time();
    iteration = continuity.get_num();
  }

  // Initialize weak formulation.
  Hermes::vector<std::string> solid_wall_markers(BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP);
  Hermes::vector<std::string> inlet_markers;
  inlet_markers.push_back(BDY_INLET);
  Hermes::vector<std::string> outlet_markers;
  outlet_markers.push_back(BDY_OUTLET);

  EulerEquationsWeakFormSemiImplicit wf(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT,solid_wall_markers, 
    inlet_markers, outlet_markers, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, (P_INIT == 0));
  EulerEquationsWeakFormStabilization wf_stabilization(&prev_rho);

  // Initialize the FE problem.
  DiscreteProblemLinear<double> dp(&wf, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  DiscreteProblem<double> dp_stabilization(&wf_stabilization, &space_stabilization);
  LinearSolver<double> solver(&dp);

  // If the FE problem is in fact a FV problem.
  if(P_INIT == 0) 
    dp.set_fvm();

  // Time stepping loop.
  for(; t < 10.0; t += time_step)
  {
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);
    CFL.set_number(0.1 + (t/7.0) * 1.0);

    // Set the current time step.
    wf.set_current_time_step(time_step);

    // Assemble the stiffness matrix and rhs.
    Hermes::Mixins::Loggable::Static::info("Assembling the stiffness matrix and right-hand side vector.");

    // Solve the matrix problem.
    Hermes::Mixins::Loggable::Static::info("Solving the matrix problem.");
    try
    {
      solver.solve();
      if(!SHOCK_CAPTURING)
      {
        Solution<double>::vector_to_solutions(solver.get_sln_vector(), Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
      }
      else
      {
        FluxLimiter* flux_limiter;
        if(SHOCK_CAPTURING_TYPE == KUZMIN)
          flux_limiter = new FluxLimiter(FluxLimiter::Kuzmin, solver.get_sln_vector(), Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e), true);
        else
          flux_limiter = new FluxLimiter(FluxLimiter::Krivodonova, solver.get_sln_vector(), Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e));

        if(SHOCK_CAPTURING_TYPE == KUZMIN)
          flux_limiter->limit_second_orders_according_to_detector();

        flux_limiter->limit_according_to_detector();

        flux_limiter->get_limited_solutions(Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
      }
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }

    CFL.calculate(Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), &mesh, time_step);

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
        pressure_view.save_numbered_screenshot("Pressure-%u.bmp", iteration - 1, true);
        Mach_number_view.save_numbered_screenshot("Mach-%u.bmp", iteration - 1, true);
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
      }
    }
  }

  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();

  return 0;
}
