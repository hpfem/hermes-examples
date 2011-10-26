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
const double TAU   = 0.05;                        
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
  int ndof = Space<double>::get_num_dofs(Hermes::vector<Space<double>*>(t_space, c_space));
  info("ndof = %d.", ndof);

  // Define initial conditions.
  InitialSolutionTemperature t_prev_time_1(&mesh, x1);
  InitialSolutionConcentration c_prev_time_1(&mesh, x1, Le);
  InitialSolutionTemperature t_prev_time_2(&mesh, x1);
  InitialSolutionConcentration c_prev_time_2(&mesh, x1, Le);
  Solution<double> t_prev_newton;
  Solution<double> c_prev_newton;

  // Filters for the reaction rate omega and its derivatives.
  CustomFilter omega(Hermes::vector<Solution<double>*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);
  CustomFilterDt omega_dt(Hermes::vector<Solution<double>*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);
  CustomFilterDc omega_dc(Hermes::vector<Solution<double>*>(&t_prev_time_1, &c_prev_time_1), Le, alpha, beta, kappa, x1, TAU);

  // Initialize visualization.
  ScalarView rview("Reaction rate", new WinGeom(0, 0, 800, 230));

  double* coeff_vec = new double[Space<double>::get_num_dofs(Hermes::vector<Space<double>*>(t_space, c_space))];
  memset(coeff_vec, 0, ndof * sizeof(double));
  Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<Space<double>*>(t_space, c_space), 
                                Hermes::vector<Solution<double> *>(&t_prev_time_1, &c_prev_time_1));

  // Initialize weak formulation.
  CustomWeakForm wf(Le, alpha, beta, kappa, x1, TAU, &omega, &omega_dt, 
                    &omega_dc, &t_prev_time_1, &c_prev_time_1, &t_prev_time_2, &c_prev_time_2);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<Space<double>*>(t_space, c_space));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
  Vector<double>* rhs = create_vector<double>(matrix_solver);
  LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver, matrix, rhs);

  // Time stepping:
  double current_time = 0.0;
  int ts = 1;
  bool jacobian_changed = true;
  do 
  {
    info("---- Time step %d, time %3.5f s", ts, current_time);

    if (NEWTON)
    {
      // Perform Newton's iteration.
      info("Solving nonlinear problem:");

      NewtonSolver<double> newton(&dp, matrix_solver);
      newton.set_verbose_output(false);
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
      Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<Space<double>*>(t_space, c_space), 
                                    Hermes::vector<Solution<double>*>(&t_prev_newton, &c_prev_newton));

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
  delete matrix;
  delete rhs;
  delete solver;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
