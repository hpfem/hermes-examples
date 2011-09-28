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
// Domain: Joukowski profile, see file domain-nurbs.xml.
//
// BC: Solid walls, inlet, outlet.
//
// IC: Constant state identical to inlet.
//
// The following parameters can be changed:
// Visualization.
const bool HERMES_VISUALIZATION = true;           // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;              // Set to "true" to enable VTK output.
const unsigned int EVERY_NTH_STEP = 1;            // Set visual output for every nth step.

// Shock capturing.
bool SHOCK_CAPTURING = false;
// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;

bool REUSE_SOLUTION = false;

const int P_INIT = 2;                                   // Initial polynomial degree.                      
const int INIT_REF_NUM_VERTEX = 1;                      // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BOUNDARY_ANISO = 7;              // Number of initial mesh refinements towards the profile.
double CFL_NUMBER = 0.5;                                // CFL value.
double time_step = 1E-4;                                // Initial time step.
const MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Equation parameters.
const double P_EXT = 7142.8571428571428571428571428571;         // Exterior pressure (dimensionless).
const double RHO_EXT = 1.0;       // Inlet density (dimensionless).   
const double V1_EXT = 0.01;       // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;        // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;         // Kappa.

// Boundary markers.
const std::string BDY_INLET = "Inlet";
const std::string BDY_OUTLET = "Outlet";
const std::string BDY_SOLID_WALL_PROFILE = "Solid Profile";
const std::string BDY_SOLID_WALL = "Solid";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2DXML mloader;
  mloader.load("domain-arcs.xml", &mesh);
  
  mesh.refine_towards_boundary(BDY_SOLID_WALL_PROFILE, INIT_REF_NUM_BOUNDARY_ANISO);
  mesh.refine_towards_vertex(0, INIT_REF_NUM_VERTEX);

  MeshView m;
  m.show(&mesh);
  m.wait_for_close();

  // Initialize boundary condition types and spaces with default shapesets.
  L2Space<double> space_rho(&mesh, P_INIT);
  L2Space<double> space_rho_v_x(&mesh, P_INIT);
  L2Space<double> space_rho_v_y(&mesh, P_INIT);
  L2Space<double> space_e(&mesh, P_INIT);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  ConstantSolution<double> prev_rho(&mesh, RHO_EXT);
  ConstantSolution<double> prev_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double> prev_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double> prev_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  // Numerical flux.
  VijayasundaramNumericalFlux num_flux(KAPPA);

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

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver_type);
  Vector<double>* rhs = create_vector<double>(matrix_solver_type);
  LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver_type, matrix, rhs);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Look for a saved solution on the disk.
  Continuity<double> continuity(Continuity<double>::onlyTime);
  int iteration = 0; double t = 0;

  if(REUSE_SOLUTION && continuity.have_record_available())
  {
    continuity.get_last_record()->load_mesh(&mesh);
    continuity.get_last_record()->load_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<SpaceType>(HERMES_L2_SPACE, HERMES_L2_SPACE, HERMES_L2_SPACE, HERMES_L2_SPACE), Hermes::vector<Mesh *>(&mesh, &mesh, 
        &mesh, &mesh));
    continuity.get_last_record()->load_solutions(Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), Hermes::vector<Mesh *>(&mesh, &mesh, 
        &mesh, &mesh));
    continuity.get_last_record()->load_time_step_length(time_step);
    t = continuity.get_last_record()->get_time();
    iteration = continuity.get_num();
  }

  // Initialize weak formulation.
  EulerEquationsWeakFormSemiImplicitMultiComponent wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL, BDY_SOLID_WALL_PROFILE, 
    BDY_INLET, BDY_OUTLET, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, (P_INIT == 0));

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));

  // If the FE problem is in fact a FV problem.
  if(P_INIT == 0) 
    dp.set_fvm();

  // Time stepping loop.
  for(; t < 3.0; t += time_step)
  {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    // Set the current time step.
    wf.set_time_step(time_step);

    // Assemble the stiffness matrix and rhs.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the matrix problem.
    info("Solving the matrix problem.");
    if(solver->solve())
      if(!SHOCK_CAPTURING)
        Solution<double>::vector_to_solutions(solver->get_sln_vector(), Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
      else
        {      
          FluxLimiter flux_limiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
            &space_rho_v_y, &space_e));

          flux_limiter.limit_second_orders_according_to_detector();

          flux_limiter.limit_according_to_detector();

          flux_limiter.get_limited_solutions(Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
        }
    else  
      error ("Matrix solver failed.\n");

    CFL.calculate_semi_implicit(Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), &mesh, time_step);

    // Visualization.

    if((iteration - 1) % EVERY_NTH_STEP == 0) 
    {
      // Hermes visualization.
      if(HERMES_VISUALIZATION) 
      {
        Mach_number.reinit();
        pressure.reinit();
        pressure_view.show(&pressure);
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
    // Save a current state on the disk.
    if(REUSE_SOLUTION)
    {
      continuity.add_record(t);
      continuity.get_last_record()->save_mesh(&mesh);
      continuity.get_last_record()->save_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e));
      continuity.get_last_record()->save_solutions(Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
      continuity.get_last_record()->save_time_step_length(time_step);
    }
  }

  View::wait();

  return 0;
}
