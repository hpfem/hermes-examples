#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

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
const bool HERMES_VISUALIZATION = false;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = true;
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

// Stability for the concentration part.
double ADVECTION_STABILITY_CONSTANT = 1.0;
const double DIFFUSION_STABILITY_CONSTANT = 1.0;

// Polynomial degree for the Euler equations (for the flow).
const int P_FLOW = 1;                        
// Polynomial degree for the concentration.
const int P_CONCENTRATION = 2;               
// CFL value.
double CFL_NUMBER = 0.1;                          
// Initial and utility time step.
double time_step_n = 1E-5, util_time_step;          

// Matrix solver for orthogonal projections: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK; 

// Number of initial uniform mesh refinements of the mesh for the flow.
unsigned int INIT_REF_NUM_FLOW = 4;               
// Number of initial uniform mesh refinements of the mesh for the concentration.
unsigned int INIT_REF_NUM_CONCENTRATION = 4;
// Number of initial mesh refinements of the mesh for the concentration towards the 
// part of the boundary where the concentration is prescribed.
unsigned int INIT_REF_NUM_CONCENTRATION_BDY = 1;

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 2.5;
// Inlet density (dimensionless).   
const double RHO_EXT = 1.0;
// Inlet x-velocity (dimensionless).
const double V1_EXT = 1.0;
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;
// Kappa.
const double KAPPA = 1.4;
// Concentration on the boundary.
const double CONCENTRATION_EXT = 0.1;
// Start time of the concentration on the boundary.
const double CONCENTRATION_EXT_STARTUP_TIME = 0.0;
// Diffusivity.
const double EPSILON = 0.01;                           

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

  // Load the mesh.
  Mesh basemesh;
  MeshReaderH2D mloader;
  mloader.load("GAMM-channel-serial.mesh", &basemesh);

  // Initialize the meshes.
  Mesh mesh_flow, mesh_concentration;
  mesh_flow.copy(&basemesh);
  mesh_concentration.copy(&basemesh);

  for(unsigned int i = 0; i < INIT_REF_NUM_CONCENTRATION; i++)
    mesh_concentration.refine_all_elements(0, true);

  mesh_concentration.refine_towards_boundary(BDY_DIRICHLET_CONCENTRATION, INIT_REF_NUM_CONCENTRATION_BDY, false);
  //mesh_flow.refine_towards_boundary(BDY_DIRICHLET_CONCENTRATION, INIT_REF_NUM_CONCENTRATION_BDY);

  for(unsigned int i = 0; i < INIT_REF_NUM_FLOW; i++)
    mesh_flow.refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  // For the concentration.
  EssentialBCs<double> bcs_concentration;

  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_DIRICHLET_CONCENTRATION, CONCENTRATION_EXT, CONCENTRATION_EXT_STARTUP_TIME));
  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_SOLID_WALL_TOP, 0.0, CONCENTRATION_EXT_STARTUP_TIME));
  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_INLET, 0.0, CONCENTRATION_EXT_STARTUP_TIME));

  L2Space<double>space_rho(&mesh_flow, P_FLOW);
  L2Space<double>space_rho_v_x(&mesh_flow, P_FLOW);
  L2Space<double>space_rho_v_y(&mesh_flow, P_FLOW);
  L2Space<double>space_e(&mesh_flow, P_FLOW);
  L2Space<double> space_stabilization(&mesh_flow, 0);
  
  // Space<double> for concentration.
  H1Space<double> space_c(&mesh_concentration, &bcs_concentration, P_CONCENTRATION);

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

  // Initialize weak formulation.
  WeakForm<double>* wf = NULL;
  if(SEMI_IMPLICIT)
  {
    wf = new EulerEquationsWeakFormSemiImplicitCoupled(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM,
    BDY_SOLID_WALL_TOP, BDY_INLET, BDY_OUTLET, BDY_NATURAL_CONCENTRATION, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c, EPSILON, P_FLOW == 0);
    
    if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
      static_cast<EulerEquationsWeakFormSemiImplicitCoupled*>(wf)->set_stabilization(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, NU_1, NU_2);
  }
  else
    wf = new EulerEquationsWeakFormExplicitCoupled(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM,
    BDY_SOLID_WALL_TOP, BDY_INLET, BDY_OUTLET, BDY_NATURAL_CONCENTRATION, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c, EPSILON, P_FLOW == 0);

  EulerEquationsWeakFormStabilization wf_stabilization(&prev_rho);

  
  // Initialize the FE problem.
  DiscreteProblem<double> dp(wf, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c));
  DiscreteProblem<double> dp_stabilization(&wf_stabilization, &space_stabilization);
  bool* discreteIndicator = NULL;
 
  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 400));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 400));
  ScalarView s5("Concentration", new WinGeom(700, 400, 600, 400));
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
  Vector<double>* rhs = create_vector<double>(matrix_solver);
  Vector<double>* rhs_stabilization = create_vector<double>(matrix_solver);
  LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver, matrix, rhs);
  
  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Set up Advection-Diffusion-Equation stability calculation class.
  ADEStabilityCalculation ADES(ADVECTION_STABILITY_CONSTANT, DIFFUSION_STABILITY_CONSTANT, EPSILON);

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 100.0; t += time_step_n)
  {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    if(SEMI_IMPLICIT)
    {
      static_cast<EulerEquationsWeakFormSemiImplicitCoupled*>(wf)->set_time_step(time_step_n);
      static_cast<EulerEquationsWeakFormSemiImplicitCoupled*>(wf)->realloc_cache(&mesh_flow);
      if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
      {
        dp_stabilization.assemble(rhs_stabilization);
        if(discreteIndicator != NULL)
          delete [] discreteIndicator;
        discreteIndicator = new bool[space_stabilization.get_mesh()->get_max_element_id() + 1];
        for(unsigned int i = 0; i < space_stabilization.get_mesh()->get_max_element_id() + 1; i++)
          discreteIndicator[i] = false;
        Element* e;
        for_all_active_elements(e, space_stabilization.get_mesh())
        {
          AsmList<double> al;
          space_stabilization.get_element_assembly_list(e, &al);
          if(rhs_stabilization->get(al.get_dof()[0]) >= 1)
            discreteIndicator[e->id] = true;
        }
        static_cast<EulerEquationsWeakFormSemiImplicitCoupled*>(wf)->set_discreteIndicator(discreteIndicator);
      }
    }
    else
      static_cast<EulerEquationsWeakFormExplicitCoupled*>(wf)->set_time_step(time_step_n);

    // Assemble stiffness matrix and rhs.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the matrix problem.
    info("Solving the matrix problem.");
    if(solver->solve())
    {
      if(!SHOCK_CAPTURING || SHOCK_CAPTURING_TYPE == FEISTAUER)
      {
        Solution<double>::vector_to_solutions(solver->get_sln_vector(), Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c), 
          Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e, &sln_c));
      }
      else
      {
        FluxLimiter* flux_limiter;
        if(SHOCK_CAPTURING_TYPE == KUZMIN)
          flux_limiter = new FluxLimiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e));
        else
          flux_limiter = new FluxLimiter(FluxLimiter::Krivodonova, solver->get_sln_vector(), Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e));

        if(SHOCK_CAPTURING_TYPE == KUZMIN)
          flux_limiter->limit_second_orders_according_to_detector();

        flux_limiter->limit_according_to_detector();

        flux_limiter->get_limited_solutions(Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
      }
    }
    else
      error ("Matrix solver failed.\n");

    util_time_step = time_step_n;
    if(SEMI_IMPLICIT)
      CFL.calculate_semi_implicit(Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), &mesh_flow, util_time_step);
    else
      CFL.calculate(Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), &mesh_flow, util_time_step);

    time_step_n = util_time_step;

    ADES.calculate(Hermes::vector<Solution<double>*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y), &mesh_concentration, util_time_step);

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&sln_rho);
    prev_rho_v_x.copy(&sln_rho_v_x);
    prev_rho_v_y.copy(&sln_rho_v_y);
    prev_e.copy(&sln_e);
    prev_c.copy(&sln_c);

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) 
    {
      // Hermes visualization.
      if(HERMES_VISUALIZATION)
      {
        Mach_number.reinit();
        Mach_number_view.show(&Mach_number);
        Mach_number_view.save_numbered_screenshot("Mach_number%i.bmp", (int)(iteration / 5), true);

        s5.show(&prev_c);
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

  return 0;
}
