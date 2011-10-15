#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// This example shows the use of subdomains. It models a round graphite object that is 
// heated through internal heat sources and cooled with water flowing past it. This 
// model is semi-realistic, double-check all parameter values and equations before 
// using it for your applications.

// If true, then just Stokes equation will be considered, not Navier-Stokes.
const bool STOKES = false;         
// If defined, pressure elements will be discontinuous (L2),
// otherwise continuous (H1).   
#define PRESSURE_IN_L2                           
// Initial polynomial degree for velocity components.
const int P_INIT_VEL = 2;          
// Initial polynomial degree for pressure.
// Note: P_INIT_VEL should always be greater than
// P_INIT_PRESSURE because of the inf-sup condition.
const int P_INIT_PRESSURE = 1;     
// Initial polynomial degree for temperature.
const int P_INIT_TEMP = 2;         
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;        
// Number of initial mesh refinements towards the hole.
const int INIT_REF_NUM_WALL = 3;   

// Domain sizes (need to be compatible with mesh file!).
// Domain height (necessary to define the parabolic
// velocity profile at inlet).
const double H = 1.0;                               
// For the calculation of Reynolds number.
const double OBSTACLE_DIAMETER = 0.3 * std::sqrt(2.0);    
// For the definition of custom initial condition.     
const double HOLE_MID_X = 0.5;
const double HOLE_MID_Y = 0.5;

// Problem parameters.
// Inlet velocity (reached after STARTUP_TIME).
const double VEL_INLET = 0.5;              
// Initial temperature.
const double TEMP_INIT_WATER = 20.0;                       
const double TEMP_INIT_GRAPHITE = 100.0;                       
// Correct is 1.004e-6 (at 20 deg Celsius) but then RE = 2.81713e+06 which 
// is too much for this simple model, so we use a larger viscosity. Note
// that kinematic viscosity decreases with rising temperature.
const double KINEMATIC_VISCOSITY_WATER = 1.004e-2;   
// We found a range of 25 - 470, but this is another 
// number that one needs to be careful about.
const double THERMAL_CONDUCTIVITY_GRAPHITE = 450;    
// At 25 deg Celsius.
const double THERMAL_CONDUCTIVITY_WATER = 0.6;       
// Density of graphite from Wikipedia, one should be 
// careful about this number.    
const double RHO_GRAPHITE = 2220;                      
const double RHO_WATER = 1000;      
// Also found on Wikipedia.
const double SPECIFIC_HEAT_GRAPHITE = 711;        
const double SPECIFIC_HEAT_WATER = 4187;            
// Heat source in graphite. This value is not realistic.
const double HEAT_SOURCE_GRAPHITE = 1e6;          
// During this time, inlet velocity increases gradually
// from 0 to VEL_INLET, then it stays constant.
const double STARTUP_TIME = 1.0;                  
// Time step.
const double time_step = 0.1;                     
// Time interval length.
const double T_FINAL = 30000.0;                   
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-4;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                   
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;  

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh_whole_domain, mesh_with_hole;
  Hermes::vector<Mesh*> meshes (&mesh_whole_domain, &mesh_with_hole);
  MeshReaderH2DXML mloader;
  mloader.load("subdomains.xml", meshes);

  /* View both meshes.
  MeshView m1("Mesh for temperature"), m2("Mesh for flow");
  m1.show(&mesh_whole_domain);
  m2.show(&mesh_with_hole);
  View::wait();
  */

  // Perform initial mesh refinements (optional).
  // Uniform.
  for(int i = 0; i < INIT_REF_NUM; i++)
    for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      meshes[meshes_i]->refine_all_elements();
  // Towards the hole.
  for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
    meshes[meshes_i]->refine_towards_boundary("Inner Wall", INIT_REF_NUM_WALL);
  // Towards the top and bottom edges.
  for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
    meshes[meshes_i]->refine_towards_boundary("Outer Wall", INIT_REF_NUM_WALL);

  // Initialize boundary conditions.
  // Flow.
  EssentialBCNonConst bc_inlet_vel_x("Inlet", VEL_INLET, H, STARTUP_TIME);
  DefaultEssentialBCConst<double> bc_other_vel_x(Hermes::vector<std::string>("Outer Wall", "Inner Wall"), 0.0);
  EssentialBCs<double> bcs_vel_x(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_inlet_vel_x, &bc_other_vel_x));
  DefaultEssentialBCConst<double> bc_vel_y(Hermes::vector<std::string>("Inlet", "Outer Wall", "Inner Wall"), 0.0);
  EssentialBCs<double> bcs_vel_y(&bc_vel_y);
  EssentialBCs<double> bcs_pressure;

  // Temperature.
  DefaultEssentialBCConst<double> bc_temperature("Inlet", 20.0);
  EssentialBCs<double> bcs_temperature(&bc_temperature);

  // Spaces for velocity components and pressure.
  H1Space<double> xvel_space(&mesh_with_hole, &bcs_vel_x, P_INIT_VEL);
  H1Space<double> yvel_space(&mesh_with_hole, &bcs_vel_y, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space<double> p_space(&mesh_with_hole, P_INIT_PRESSURE);
#else
  H1Space<double> p_space(&mesh_with_hole, &bcs_pressure, P_INIT_PRESSURE);
#endif
  // Space<double> for temperature.
  H1Space<double> temperature_space(&mesh_whole_domain, &bcs_temperature, P_INIT_TEMP);

  // Calculate and report the number of degrees of freedom.
  int ndof = Space<double>::get_num_dofs(Hermes::vector<Space<double> *>(&xvel_space, &yvel_space, &p_space, &temperature_space));
  info("ndof = %d.", ndof);

  // Define projection norms.
  ProjNormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  ProjNormType p_proj_norm = HERMES_L2_NORM;
#else
  ProjNormType p_proj_norm = HERMES_H1_NORM;
#endif
  ProjNormType temperature_proj_norm = HERMES_H1_NORM;

  // Solutions for the Newton's iteration and time stepping.
  info("Setting initial conditions.");
  ZeroSolution xvel_prev_time(&mesh_with_hole), yvel_prev_time(&mesh_with_hole), p_prev_time(&mesh_with_hole);

  // Initial solution for temperature.
  CustomInitialCondition temperature_ic(&mesh_whole_domain, HOLE_MID_X, HOLE_MID_Y, 
      0.5*OBSTACLE_DIAMETER, TEMP_INIT_WATER, TEMP_INIT_GRAPHITE); 

  // Project it to obtain continuous initial condition for temperature.
  double *coeff_vec_temp = new double[temperature_space.get_num_dofs()];
  OGProjection<double>::project_global(&temperature_space, &temperature_ic, coeff_vec_temp, matrix_solver_type, HERMES_L2_NORM);
  Solution<double> temperature_prev_time;
  Solution<double>::vector_to_solution(coeff_vec_temp, &temperature_space, &temperature_prev_time);
  delete [] coeff_vec_temp;

  // Calculate Reynolds number.
  double reynolds_number = VEL_INLET * OBSTACLE_DIAMETER / KINEMATIC_VISCOSITY_WATER;
  info("RE = %g", reynolds_number);

  // Initialize weak formulation.
  CustomWeakFormHeatAndFlow wf(STOKES, reynolds_number, time_step, &xvel_prev_time, &yvel_prev_time, &temperature_prev_time, 
      HEAT_SOURCE_GRAPHITE, SPECIFIC_HEAT_GRAPHITE, SPECIFIC_HEAT_WATER, RHO_GRAPHITE, RHO_WATER, 
      THERMAL_CONDUCTIVITY_GRAPHITE, THERMAL_CONDUCTIVITY_WATER);
  
  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<Space<double> *>(&xvel_space, &yvel_space, &p_space, &temperature_space));

  // Initialize the Newton solver.
  NewtonSolver<double> newton(&dp, matrix_solver_type);

  // Initialize views.
  Views::VectorView vview("velocity [m/s]", new Views::WinGeom(0, 0, 700, 360));
  Views::ScalarView pview("pressure [Pa]", new Views::WinGeom(0, 415, 700, 350));
  Views::ScalarView tempview("temperature [C]", new Views::WinGeom(0, 795, 700, 350));
  //vview.set_min_max_range(0, 0.5);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(false);
  tempview.fix_scale_width(80);
  tempview.show_mesh(false);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  double* coeff_vec = new double[ndof];
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  OGProjection<double>::project_global(Hermes::vector<Space<double> *>(&xvel_space, &yvel_space, &p_space, &temperature_space), 
    Hermes::vector<MeshFunction<double> *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time, &temperature_prev_time), 
    coeff_vec, matrix_solver_type, 
    Hermes::vector<ProjNormType>(vel_proj_norm, vel_proj_norm, p_proj_norm, temperature_proj_norm));

  // Time-stepping loop:
  char title[100];
  int num_time_steps = T_FINAL / time_step;
  double current_time = 0.0;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    current_time += time_step;
    info("---- Time step %d, time = %g:", ts, current_time);

    // Update time-dependent essential BCs.
    if (current_time <= STARTUP_TIME) 
    {
      info("Updating time-dependent essential BC.");
      Space<double>::update_essential_bc_values(Hermes::vector<Space<double> *>(&xvel_space, &yvel_space, &p_space, 
                                                &temperature_space), current_time);
    }

    // Perform Newton's iteration.
    info("Solving nonlinear problem:");
    bool verbose = true;
    // Perform Newton's iteration and translate the resulting coefficient vector into previous time level solutions.
    newton.set_verbose_output(verbose);
    try
    {
      newton.solve(coeff_vec, NEWTON_TOL, NEWTON_MAX_ITER);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    };
    {
      Hermes::vector<Solution<double> *> tmp(&xvel_prev_time, &yvel_prev_time, &p_prev_time, &temperature_prev_time);
      Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<Space<double> *>(&xvel_space, 
          &yvel_space, &p_space, &temperature_space), tmp);
    }
    
    // Show the solution at the end of time step.
    sprintf(title, "Velocity [m/s], time %g s", current_time);
    vview.set_title(title);
    vview.show(&xvel_prev_time, &yvel_prev_time);
    sprintf(title, "Pressure [Pa], time %g s", current_time);
    pview.set_title(title);
    pview.show(&p_prev_time);
    sprintf(title, "Temperature [C], time %g s", current_time);
    tempview.set_title(title);
    tempview.show(&temperature_prev_time,  Views::HERMES_EPS_HIGH);
  }

  delete [] coeff_vec;

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
