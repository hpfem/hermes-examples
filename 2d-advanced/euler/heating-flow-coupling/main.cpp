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
// Domain: A square, see file square.mesh.
//
// BC: Solid walls, inlet, no outlet.
//
// IC: Constant state identical to inlet, only with higher pressure.
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

// For saving/loading of solution.
bool REUSE_SOLUTION = false;

// Initial polynomial degree.
const int P_INIT = 1;                                                         
// Number of initial uniform mesh refinements.    
const int INIT_REF_NUM = 4;                                                
// CFL value.
double CFL_NUMBER = 0.7;                                
// Initial time step.
double time_step = 1E-4;                                
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
const MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 1.5;           
// Initial pressure (dimensionless).
const double P_INITIAL_HIGH = 1.5;  
// Initial pressure (dimensionless).
const double P_INITIAL_LOW = 1.0;   
// Inlet density (dimensionless).   
const double RHO_EXT = 0.5;         
// Initial density (dimensionless).   
const double RHO_INITIAL_HIGH = 0.5;
// Initial density (dimensionless).   
const double RHO_INITIAL_LOW = 0.3; 
// Inlet x-velocity (dimensionless).
const double V1_EXT = 0.0;          
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;          
// Kappa.
const double KAPPA = 1.4;           

// Stability for the concentration part.
double ADVECTION_STABILITY_CONSTANT = 0.1;
const double DIFFUSION_STABILITY_CONSTANT = 0.1;

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
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) 
    mesh.refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  Hermes2D::DefaultEssentialBCConst<double> bc_temp_zero("Solid", 0.0);
  EssentialBCs<double> bcs(&bc_temp_zero);

  L2Space<double> space_rho(&mesh, P_INIT);
  L2Space<double> space_rho_v_x(&mesh, P_INIT);
  L2Space<double> space_rho_v_y(&mesh, P_INIT);
  L2Space<double> space_e(&mesh, P_INIT);
  H1Space<double> space_temp(&mesh, &bcs, 1);
  L2Space<double> space_stabilization(&mesh, 0);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_temp));
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  InitialSolutionLinearProgress prev_rho(&mesh, RHO_INITIAL_HIGH, RHO_INITIAL_LOW, MESH_SIZE);
  ConstantSolution<double> prev_rho_v_x(&mesh, 0.0);
  ConstantSolution<double> prev_rho_v_y(&mesh, 0.0);
  InitialSolutionLinearProgress prev_e(&mesh, QuantityCalculator::calc_energy(RHO_INITIAL_HIGH, RHO_INITIAL_HIGH * V1_EXT, RHO_INITIAL_HIGH * V2_EXT, P_INITIAL_HIGH, KAPPA), QuantityCalculator::calc_energy(RHO_INITIAL_LOW, RHO_INITIAL_LOW * V1_EXT, RHO_INITIAL_LOW * V2_EXT, P_INITIAL_LOW, KAPPA), MESH_SIZE);
  ConstantSolution<double> prev_temp(&mesh, 0.0);

  // Filters for visualization of Mach number, pressure and entropy.
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  VelocityFilter vel_x(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x));
  VelocityFilter vel_y(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_y));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  VectorView velocity_view("Velocity", new WinGeom(0, 400, 600, 300));
  ScalarView velocity_view_x("Velocity - x", new WinGeom(0, 400, 600, 300));
  ScalarView velocity_view_y("Velocity - y", new WinGeom(700, 400, 600, 300));
  ScalarView density_view("Density", new WinGeom(1400, 0, 600, 300));
  ScalarView temperature_view("Temperature", new WinGeom(1400, 400, 600, 300));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>();
  Vector<double>* rhs = create_vector<double>();
  Vector<double>* rhs_stabilization = create_vector<double>();
  LinearMatrixSolver<double>* solver = create_linear_solver<double>( matrix, rhs);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);
  ADEStabilityCalculation ADES(ADVECTION_STABILITY_CONSTANT, DIFFUSION_STABILITY_CONSTANT, 1e2);

  // Look for a saved solution on the disk.
  int iteration = 0; double t = 0;

  // Initialize weak formulation.
  Hermes::vector<std::string> solid_wall_markers;
  solid_wall_markers.push_back(BDY_SOLID_WALL);
  Hermes::vector<std::string> inlet_markers;
  inlet_markers.push_back(BDY_INLET);
  Hermes::vector<std::string> outlet_markers;

  EulerEquationsWeakFormSemiImplicitCoupledWithHeat wf(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, solid_wall_markers, 
    inlet_markers, outlet_markers, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_temp);
  
  EulerEquationsWeakFormStabilization wf_stabilization(&prev_rho);

  if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
    wf.set_stabilization(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, NU_1, NU_2);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_temp));
  DiscreteProblem<double> dp_stabilization(&wf_stabilization, &space_stabilization);

  // Time stepping loop.
  for(; t < 10.0; t += time_step)
  {
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

    if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
    {
      assert(space_stabilization.get_num_dofs() == space_stabilization.get_mesh()->get_num_active_elements());
      dp_stabilization.assemble(rhs_stabilization);
      bool* discreteIndicator = new bool[space_stabilization.get_num_dofs()];
      memset(discreteIndicator, 0, space_stabilization.get_num_dofs() * sizeof(bool));
      Element* e;
      for_all_active_elements(e, space_stabilization.get_mesh())
      {
        AsmList<double> al;
        space_stabilization.get_element_assembly_list(e, &al);
        if(rhs_stabilization->get(al.get_dof()[0]) >= 1)
          discreteIndicator[e->id] = true;
      }
      wf.set_discreteIndicator(discreteIndicator, space_stabilization.get_num_dofs());
    }

    // Set the current time step.
    wf.set_current_time_step(time_step);

    // Assemble the stiffness matrix and rhs.
    Hermes::Mixins::Loggable::Static::info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the matrix problem.
    Hermes::Mixins::Loggable::Static::info("Solving the matrix problem.");
    if(solver->solve())
    {
      if(!SHOCK_CAPTURING || SHOCK_CAPTURING_TYPE == FEISTAUER)
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
        temperature_view.show(&prev_temp);
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
