#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// Flow inside a rotating circle. Both the flow and the circle are not moving 
// at the beginning. As the circle starts to rotate at increasing speed. also 
// the flow starts to move. We use time-dependent laminar incompressible Navier-Stokes 
// equations discretized in time via the implicit Euler method. The Newton's method 
// is used to solve the nonlinear problem at each time step.
//
// PDE: incompressible Navier-Stokes equations in the form
//     \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
//     div v = 0.
//
// BC: tangential velocity V on the boundary.
//
// Geometry: A circular domain.
//
// The following parameters can be changed:

// Number of initial uniform mesh refinements. 
const int INIT_REF_NUM = 2;                       
// Number of initial mesh refinements towards boundary. 
const int INIT_BDY_REF_NUM_INNER = 2;             
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
const double RE = 5000.0;                         
// Surface velocity of inner circle.
const double VEL = 0.1;                           
// During this time, surface velocity of the inner circle increases 
// gradually from 0 to VEL, then it stays constant.
const double STARTUP_TIME = 1.0;                  
// Time step.
const double TAU = 10.0;                          
// Time interval length.
const double T_FINAL = 3600.0;                    
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-5;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 10;                   
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Current time (used in weak forms).
double current_time = 0;

// Custom function to calculate drag coefficient.
double integrate_over_wall(MeshFunction<double>* meshfn, int marker)
{
  Quad2D* quad = &g_quad_2d_std;
  meshfn->set_quad_2d(quad);

  double integral = 0.0;
  Element* e;
  Mesh* mesh = meshfn->get_mesh();

  for_all_active_elements(e, mesh)
  {
    for(int edge = 0; edge < e->get_num_surf(); edge++)
    {
      if ((e->en[edge]->bnd) && (e->en[edge]->marker == marker))
      {
        update_limit_table(e->get_mode());
        RefMap* ru = meshfn->get_refmap();

        meshfn->set_active_element(e);
        int eo = quad->get_edge_points(edge);
        meshfn->set_quad_order(eo, H2D_FN_VAL);
        double *uval = meshfn->get_fn_values();
        double3* pt = quad->get_points(eo);
        double3* tan = ru->get_tangent(edge);
        for (int i = 0; i < quad->get_num_points(eo); i++)
          integral += pt[i][2] * uval[i] * tan[i][2];
      }
    }
  }
  return integral * 0.5;
}

int main(int argc, char* argv[])
{
  

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  //mloader.load("domain-concentric.mesh", &mesh);

  // Initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(Hermes::vector<std::string>("Bdy-1", "Bdy-2", "Bdy-3", "Bdy-4"), 
                               INIT_BDY_REF_NUM_INNER, false);  // True for anisotropic refinement.

  // Initialize boundary conditions.
  EssentialBCNonConstX bc_vel_x(Hermes::vector<std::string>("Bdy-1", "Bdy-2", "Bdy-3","Bdy-4"), VEL, STARTUP_TIME);
  EssentialBCNonConstY bc_vel_y(Hermes::vector<std::string>("Bdy-1", "Bdy-2", "Bdy-3","Bdy-4"), VEL, STARTUP_TIME);
  EssentialBCs<double> bcs_vel_x(&bc_vel_x);
  EssentialBCs<double> bcs_vel_y(&bc_vel_y);

  // Spaces for velocity components and pressure.
  H1Space<double> xvel_space(&mesh, &bcs_vel_x, P_INIT_VEL);
  H1Space<double> yvel_space(&mesh, &bcs_vel_y, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space<double> p_space(&mesh, P_INIT_PRESSURE);
#else
  H1Space<double> p_space(&mesh, P_INIT_PRESSURE);
#endif
  Hermes::vector<const Space<double>*> spaces = Hermes::vector<const Space<double>*>(&xvel_space, &yvel_space, &p_space);

  // Calculate and report the number of degrees of freedom.
  int ndof = Space<double>::get_num_dofs(spaces);
  info("ndof = %d.", ndof);

  // Define projection norms.
  ProjNormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  ProjNormType p_proj_norm = HERMES_L2_NORM;
#else
  ProjNormType p_proj_norm = HERMES_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  info("Setting initial conditions.");
  ZeroSolution xvel_prev_time(&mesh);
  ZeroSolution yvel_prev_time(&mesh);
  ZeroSolution p_prev_time(&mesh);
  Hermes::vector<Solution<double>*> slns_prev_time = 
      Hermes::vector<Solution<double>*>(&xvel_prev_time, &yvel_prev_time, &p_prev_time);

  // Initialize weak formulation.
  WeakForm<double>* wf = new WeakFormNSNewton(STOKES, RE, TAU, &xvel_prev_time, &yvel_prev_time);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(wf, spaces);

  // Initialize views.
  VectorView vview("velocity [m/s]", new WinGeom(0, 0, 600, 500));
  ScalarView pview("pressure [Pa]", new WinGeom(610, 0, 600, 500));
  //vview.set_min_max_range(0, 1.6);
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
    info("Updating time-dependent essential BC.");
    Space<double>::update_essential_bc_values(Hermes::vector<Space<double>*>(&xvel_space, &yvel_space), current_time);

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
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces, slns_prev_time);

    // Show the solution at the end of time step.
    sprintf(title, "Velocity, time %g", current_time);
    vview.set_title(title);
    vview.show(&xvel_prev_time, &yvel_prev_time);
    sprintf(title, "Pressure, time %g", current_time);
    pview.set_title(title);
    pview.show(&p_prev_time);
 }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
