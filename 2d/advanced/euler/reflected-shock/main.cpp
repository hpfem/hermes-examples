#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// This example solves the compressible Euler equations using a basic
// piecewise-constant finite volume method.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: channel.
//
// BC: essential boundary conditions at the left and the top boundaries
//     slip condition at the bottom part, no condition on the right.
//
// IC: Exact solution.
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

const int P_INIT = 0;                                   // Initial polynomial degree.                      
const int INIT_REF_NUM = 4;                             // Number of initial uniform mesh refinements.                       
double CFL_NUMBER = 1.0;                                // CFL value.
double time_step = 1E-4;                                // Initial time step.
const MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

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
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("channel.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) 
    mesh.refine_all_elements(0);
  
  // Initialize boundary condition types and spaces with default shapesets.
  L2Space<double> space_rho(&mesh, P_INIT);
  L2Space<double> space_rho_v_x(&mesh, P_INIT);
  L2Space<double> space_rho_v_y(&mesh, P_INIT);
  L2Space<double> space_e(&mesh, P_INIT);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  InitialSolutionEulerDensity prev_rho(&mesh, RHO_INIT);
  InitialSolutionEulerDensityVelX prev_rho_v_x(&mesh, RHO_INIT * V1_INIT);
  InitialSolutionEulerDensityVelY prev_rho_v_y(&mesh, RHO_INIT * V2_INIT);
  InitialSolutionEulerDensityEnergy prev_e(&mesh, QuantityCalculator::calc_energy(RHO_INIT, RHO_INIT * V1_INIT, RHO_INIT * V2_INIT, PRESSURE_INIT, KAPPA));

  // Numerical flux.
  OsherSolomonNumericalFlux num_flux(KAPPA);

  // Initialize weak formulation.
  EulerEquationsWeakFormSemiImplicitMultiComponentTwoInflows wf(&num_flux, KAPPA, RHO_LEFT, V1_LEFT, V2_LEFT, PRESSURE_LEFT, RHO_TOP, V1_TOP, V2_TOP, PRESSURE_TOP, BDY_SOLID_WALL, BDY_INLET_LEFT, BDY_INLET_TOP, BDY_OUTLET,
    &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, (P_INIT == 0));

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));

  // If the FE problem is in fact a FV problem.
  if(P_INIT == 0) 
    dp.set_fvm();  

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver_type);
  Vector<double>* rhs = create_vector<double>(matrix_solver_type);
  LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver_type, matrix, rhs);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 3.0; t += time_step)
  {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    // Set the current time step.
    wf.set_time_step(time_step);

    // Assemble the stiffness matrix and rhs.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    std::ofstream out("out");
    for(int i = 0; i < matrix->get_size(); i++)
      for(int j = 0; j < matrix->get_size(); j++)
        out << matrix->get(i, j) << std::endl;
    out.close();

    dp.get_last_profiling_output(std::cout);

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
  Mach_number_view.close();

  return 0;
}
