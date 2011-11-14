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
shockCapturingType SHOCK_CAPTURING_TYPE = KUZMIN;
// Quantitative parameter of the discontinuity detector in case of Krivodonova.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;
// Quantitative parameter of the shock capturing in case of Feistauer.
const double NU_1 = 0.1;
const double NU_2 = 0.1;

// For saving/loading of solution.
bool REUSE_SOLUTION = true;

// Initial polynomial degree.       
const int P_INIT = 1;                                                  
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM_VERTEX = 2;                      
// Number of initial mesh refinements towards the profile.
const int INIT_REF_NUM_BOUNDARY_ANISO = 6;              
// CFL value.
double CFL_NUMBER = 0.1;                                
// Initial time step.
double time_step_n = 1E-6;                                
// Initial time step.
double time_step_n_minus_one = 1E-6;                                

// Matrix solver for orthogonal projections: 
// SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 7142.8571428571428571428571428571;         
// Inlet density (dimensionless).   
const double RHO_EXT = 1.0;       
// Inlet x-velocity (dimensionless).
const double V1_EXT = 0.01;       
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;        
// Kappa.
const double KAPPA = 1.4;        

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

  // Initialize boundary condition types and spaces with default shapesets.
  L2Space<double> space_rho(&mesh, P_INIT);
  L2Space<double> space_rho_v_x(&mesh, P_INIT);
  L2Space<double> space_rho_v_y(&mesh, P_INIT);
  L2Space<double> space_e(&mesh, P_INIT);
  L2Space<double> space_stabilization(&mesh, 0);
  int ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  info("ndof: %d", ndof);
        
  // Initialize solutions, set initial conditions.
  ConstantSolution<double> prev_rho(&mesh, RHO_EXT);
  ConstantSolution<double> prev_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double> prev_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double> prev_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  ConstantSolution<double> prev_rho2(&mesh, RHO_EXT);
  ConstantSolution<double> prev_rho_v_x2(&mesh, RHO_EXT * V1_EXT);
  ConstantSolution<double> prev_rho_v_y2(&mesh, RHO_EXT * V2_EXT);
  ConstantSolution<double> prev_e2(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

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
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
  Vector<double>* rhs = create_vector<double>(matrix_solver);
  Vector<double>* rhs_stabilization = create_vector<double>(matrix_solver);
  LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver, matrix, rhs);

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
    continuity.get_last_record()->load_solutions(Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_rho2, &prev_rho_v_x2, &prev_rho_v_y2, &prev_e2), Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e, &space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
    continuity.get_last_record()->load_time_step_length(time_step_n);
    continuity.get_last_record()->load_time_step_length_n_minus_one(time_step_n_minus_one);
    t = continuity.get_last_record()->get_time();
    iteration = continuity.get_num();
  }

  // Initialize weak formulation.
  EulerEquationsWeakFormSemiImplicitMultiComponent2ndOrder wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL, BDY_SOLID_WALL_PROFILE, 
    BDY_INLET, BDY_OUTLET, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_rho2, &prev_rho_v_x2, &prev_rho_v_y2, &prev_e2, (P_INIT == 0));

  EulerEquationsWeakFormStabilization wf_stabilization(&prev_rho);

  if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
    wf.set_stabilization(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_rho2, &prev_rho_v_x2, &prev_rho_v_y2, &prev_e2, NU_1, NU_2);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<const Space<double>*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  DiscreteProblem<double> dp_stabilization(&wf_stabilization, &space_stabilization);

  // If the FE problem is in fact a FV problem.
  if(P_INIT == 0) 
    dp.set_fvm();

  // Time stepping loop.
  for(; t < 5.0; t += time_step_n)
  {
    info("---- Time step %d, time %3.5f.", iteration++, t);
    CFL.set_number(0.1 + (t/2.5) * 1000.0);

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
      wf.set_discreteIndicator(discreteIndicator);
    }
    
    // Set the current time step.
    wf.set_time_step(time_step_n, time_step_n_minus_one);

    // Assemble the stiffness matrix and rhs.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the matrix problem.
    info("Solving the matrix problem.");
    if(solver->solve())
      {
        if(iteration > 1)
        {
          prev_rho2.copy(&prev_rho);
          prev_rho_v_x2.copy(&prev_rho_v_x);
          prev_rho_v_y2.copy(&prev_rho_v_y);
          prev_e2.copy(&prev_e);
        }

      if(!SHOCK_CAPTURING || SHOCK_CAPTURING_TYPE == FEISTAUER)
      {
        Solution<double>::vector_to_solutions(solver->get_sln_vector(), Hermes::vector<const Space<double> *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e), Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
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

    time_step_n_minus_one = time_step_n;
    CFL.calculate_semi_implicit(Hermes::vector<Solution<double> *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), &mesh, time_step_n);
    
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
    // Save a current state on the disk.
    if(iteration > 1)
    {
    continuity.add_record(t);
    continuity.get_last_record()->save_mesh(&mesh);
    continuity.get_last_record()->save_spaces(Hermes::vector<Space<double> *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e, &space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
      continuity.get_last_record()->save_solutions(Hermes::vector<Solution<double>*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_rho2, &prev_rho_v_x2, &prev_rho_v_y2, &prev_e2));
    continuity.get_last_record()->save_time_step_length(time_step_n);
      continuity.get_last_record()->save_time_step_length_n_minus_one(time_step_n_minus_one);
    }
  }

  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();

  return 0;
}
