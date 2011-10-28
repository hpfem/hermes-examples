#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This example is a very simple flame propagation model (laminar flame,
//  zero flow velocity), and its purpose is to show how the Newton's method
//  is applied to a time-dependent two-equation system.
//
//  PDEs:
//
//  dT/dt - laplace T = omega(T,Y),
//  dY/dt - 1/Le * laplace Y = - omega(T,Y).
//
//  Domain: rectangle with cooled rods.
//
//  BC:  T = 1, Y = 0 on the inlet,
//       dT/dn = - kappa T on cooled rods,
//       dT/dn = 0, dY/dn = 0 elsewhere.
//
//  Time-stepping: a second order BDF formula.

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;                       
// Initial polynomial degree.
const int P_INIT = 2;                             

// Newton's method.
const double DAMPING_COEFF = 1.0;
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-4;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 50;                   
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  
// If NEWTON == true then the Newton's iteration is performed.
// in every time step. Otherwise the convective term is linearized
// using the velocities from the previous time step.
const bool NEWTON = true;                         
// Preconditioning by jacobian (1) (less GMRES iterations, more time to create precond)
// or by approximation of jacobian (2) (less time for precond creation, more GMRES iters).
// in case of jfnk, default Ifpack proconditioner in case of Newton.
const int PRECOND = 2;                            

// Problem constants.
// Time step.
const double TAU   = 0.005;                        
// Time interval length.
const double T_FINAL = 60.0;                      
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> left_t("Left", 1.0);
  EssentialBCs<double> bcs_t(&left_t);

  DefaultEssentialBCConst<double> left_c("Left", 0.0);
  EssentialBCs<double> bcs_c(&left_c);

  // Create H1 spaces with default shapesets.
  H1Space<double>* t_space = new H1Space<double>(&mesh, &bcs_t, P_INIT);
  H1Space<double>* c_space = new H1Space<double>(&mesh, &bcs_c, P_INIT);
  Hermes::vector<Space<double> *> spaces = Hermes::vector<Space<double> *>(t_space, c_space);
  int ndof = Space<double>::get_num_dofs(spaces);
  info("ndof = %d.", ndof);

  // Define initial conditions.
  InitialSolutionTemperature t_prev_time_1(&mesh, x1);
  InitialSolutionConcentration c_prev_time_1(&mesh, x1, Le);
  Hermes::vector<MeshFunction<double>*> meshfns_prev_time_1 = Hermes::vector<MeshFunction<double>*>(&t_prev_time_1, &c_prev_time_1);
  Hermes::vector<Solution<double>*> slns_prev_time_1 = Hermes::vector<Solution<double>*>(&t_prev_time_1, &c_prev_time_1);
  InitialSolutionTemperature t_prev_time_2(&mesh, x1);
  InitialSolutionConcentration c_prev_time_2(&mesh, x1, Le);
  Solution<double> t_prev_newton(&mesh);
  Solution<double> c_prev_newton(&mesh);
  Hermes::vector<Solution<double>*> slns_prev_newton = Hermes::vector<Solution<double>*>(&t_prev_newton, &c_prev_newton);

  // Filters for the reaction rate omega and its derivatives.
  CustomFilter omega(slns_prev_time_1, Le, alpha, beta, kappa, x1, TAU);
  CustomFilterDt omega_dt(slns_prev_time_1, Le, alpha, beta, kappa, x1, TAU);
  CustomFilterDc omega_dc(slns_prev_time_1, Le, alpha, beta, kappa, x1, TAU);

  // Initialize visualization.
  ScalarView rview("Reaction rate", new WinGeom(0, 0, 800, 230));

  // Initialize weak formulation.
  CustomWeakForm wf(Le, alpha, beta, kappa, x1, TAU, &omega, &omega_dt, 
                    &omega_dc, &t_prev_time_1, &c_prev_time_1, &t_prev_time_2, &c_prev_time_2);

  // Project the functions "t_prev_time_1" and "c_prev_time_1" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solutions on the FE meshes.");
  double* coeff_vec = new double[ndof];
  OGProjection<double>::project_global(spaces, meshfns_prev_time_1, coeff_vec);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, spaces);

  // Time stepping:
  double current_time = 0.0;
  int ts = 1;
  bool jacobian_changed = true;
  do 
  {
    info("Time step %d, time %3.5f s", ts, current_time);

    if (NEWTON)
    {
      // Perform Newton's iteration.
      info("Solving nonlinear problem:");

      NewtonSolver<double> newton(&dp, matrix_solver);
      newton.set_verbose_output(true);
      try
      {
        newton.solve(coeff_vec, NEWTON_TOL, NEWTON_MAX_ITER);
      }
      catch(Hermes::Exceptions::Exception e)
      {
        e.printMsg();
        error("Newton's iteration failed.");
      };
      // Translate the resulting coefficient vector into the instance of Solution.
      Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces, slns_prev_newton);

      // Saving solutions for the next time step.
      if(ts > 1)
      {
        t_prev_time_2.copy(&t_prev_time_1);
        c_prev_time_2.copy(&c_prev_time_1);
      }

      t_prev_time_1.copy(&t_prev_newton);
      c_prev_time_1.copy(&c_prev_newton);

    }

    // Visualization.
    rview.set_min_max_range(0.0,2.0);
    rview.show(&omega);

    // Increase current time and time step counter.
    current_time += TAU;
    ts++;
  }
  while (current_time < T_FINAL);

  // Clean up.
  delete [] coeff_vec;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
