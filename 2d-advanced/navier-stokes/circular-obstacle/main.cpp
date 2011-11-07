#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// The time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method. If NEWTON == true,
// the Newton's method is used to solve the nonlinear problem at each time
// step. We also show how to use discontinuous ($L^2$) elements for pressure 
// and thus make the velocity discreetely divergence free. Comparison to 
// approximating the pressure with the standard (continuous) Taylor-Hood elements 
// is enabled. The Reynolds number Re = 200 which is embarrassingly low. You
// can increase it but then you will need to make the mesh finer, and the
// computation will take more time.
//
// PDE: incompressible Navier-Stokes equations in the form
// \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
// div v = 0.
//
// BC: u_1 is a time-dependent constant and u_2 = 0 on Gamma_4 (inlet),
//     u_1 = u_2 = 0 on Gamma_1 (bottom), Gamma_3 (top) and Gamma_5 (obstacle),
//     "do nothing" on Gamma_2 (outlet).
//
// Geometry: Rectangular channel containing an off-axis circular obstacle. The
//           radius and position of the circle, as well as other geometry
//           parameters can be changed in the mesh file "domain.mesh".
//
// The following parameters can be changed:

// For application of Stokes flow (creeping flow).
const bool STOKES = false;                        
// If this is defined, the pressure is approximated using
                                                  // discontinuous L2 elements (making the velocity discreetely
                                                  // divergence-free, more accurate than using a continuous
                                                  // pressure approximation). Otherwise the standard continuous
                                                  // elements are used. The results are striking - check the
                                                  // tutorial for comparisons.
#define PRESSURE_IN_L2                            
// Initial polynomial degree for velocity components.
const int P_INIT_VEL = 2;                         
// Initial polynomial degree for pressure.
// Note: P_INIT_VEL should always be greater than
// P_INIT_PRESSURE because of the inf-sup condition.
const int P_INIT_PRESSURE = 1;                    
// Reynolds number.
const double RE = 200.0;                          
// Inlet velocity (reached after STARTUP_TIME).
const double VEL_INLET = 1.0;                     
// During this time, inlet velocity increases gradually
// from 0 to VEL_INLET, then it stays constant.
const double STARTUP_TIME = 1.0;                  
// Time step.
const double TAU = 0.1;                           
// Time interval length.
const double T_FINAL = 30000.0;                   
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-3;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 10;                   
// Domain height (necessary to define the parabolic
// velocity profile at inlet).
const double H = 5;                               
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Boundary markers.
const std::string BDY_BOTTOM = "b1";
const std::string BDY_RIGHT = "b2";
const std::string BDY_TOP = "b3";
const std::string BDY_LEFT = "b4";
const std::string BDY_OBSTACLE = "b5";

// Current time (used in weak forms).
double current_time = 0;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_OBSTACLE, 4, false);
  // '4' is the number of levels,
  mesh.refine_towards_boundary(BDY_TOP, 4, true);     
  // 'true' stands for anisotropic refinements.
  mesh.refine_towards_boundary(BDY_BOTTOM, 4, true);  

  // Initialize boundary conditions.
  EssentialBCNonConst bc_left_vel_x(BDY_LEFT, VEL_INLET, H, STARTUP_TIME);
  DefaultEssentialBCConst<double> bc_other_vel_x(Hermes::vector<std::string>(BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  EssentialBCs<double> bcs_vel_x(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_left_vel_x, &bc_other_vel_x));
  DefaultEssentialBCConst<double> bc_vel_y(Hermes::vector<std::string>(BDY_LEFT, BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  EssentialBCs<double> bcs_vel_y(&bc_vel_y);

  // Spaces for velocity components and pressure.
  H1Space<double> xvel_space(&mesh, &bcs_vel_x, P_INIT_VEL);
  H1Space<double> yvel_space(&mesh, &bcs_vel_y, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space<double> p_space(&mesh, P_INIT_PRESSURE);
#else
  H1Space<double> p_space(&mesh, P_INIT_PRESSURE);
#endif
  Hermes::vector<Space<double>* > spaces(&xvel_space, &yvel_space, &p_space);
  Hermes::vector<const Space<double>* > spaces_const(&xvel_space, &yvel_space, &p_space);

  // Calculate and report the number of degrees of freedom.
  int ndof = Space<double>::get_num_dofs(spaces_const);
  info("ndof = %d.", ndof);

  // Define projection norms.
  ProjNormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  ProjNormType p_proj_norm = HERMES_L2_NORM;
#else
  ProjNormType p_proj_norm = HERMES_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  info("Setting zero initial conditions.");
  ZeroSolution xvel_prev_time(&mesh);
  ZeroSolution yvel_prev_time(&mesh);
  ZeroSolution p_prev_time(&mesh);
  Hermes::vector<Solution<double>* > slns_prev_time = Hermes::vector<Solution<double>* >(&xvel_prev_time, &yvel_prev_time, &p_prev_time);

  // Initialize weak formulation.
  WeakForm<double>* wf;
  wf = new WeakFormNSNewton(STOKES, RE, TAU, &xvel_prev_time, &yvel_prev_time);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(wf, spaces_const);

  // Initialize views.
  VectorView vview("velocity [m/s]", new WinGeom(0, 0, 750, 240));
  ScalarView pview("pressure [Pa]", new WinGeom(0, 290, 750, 240));
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(true);

  // Time-stepping loop:
  char title[100];
  int num_time_steps = T_FINAL / TAU;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    current_time += TAU;
    info("---- Time step %d, time = %g:", ts, current_time);

    // Update time-dependent essential BCs.
    if (current_time <= STARTUP_TIME) {
      info("Updating time-dependent essential BC.");
      Space<double>::update_essential_bc_values(spaces, current_time);
    }

    // Perform Newton's iteration.
    info("Solving nonlinear problem:");
    Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver);
    try
    {
      newton.solve(NULL, NEWTON_TOL, NEWTON_MAX_ITER);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    };

    // Update previous time level solutions.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces_const, slns_prev_time);

    // Show the solution at the end of time step.
    sprintf(title, "Velocity, time %g", current_time);
    vview.set_title(title);
    vview.show(&xvel_prev_time, &yvel_prev_time, HERMES_EPS_LOW);
    sprintf(title, "Pressure, time %g", current_time);
    pview.set_title(title);
    pview.show(&p_prev_time);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
