#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This is a simple version of the quantum billiard problem, described with 
//  a complex-valued linear time-dependent wave equation. The equation is converted 
//  into a system of two PDE which are first-order in time. For time-discretization 
//  one can use either the first-order implicit Euler method or the second-order 
//  Crank-Nicolson method.
//
//  PDE: 
//
//  \frac{\partial^2}{\partial t^2}\Psi - \frac{\partial^2}{\partial x^2}\Psi
//  - \frac{\partial^2}{\partial y^2}\Psi = 0.
//  
//  Domain: square (-1, 1)^2.
//
//  BC: homogeneous Dirichlet everywhere on the boundary.
//
//  IC: Gaussian distribution \Psi(0, x, y) = 

const int INIT_REF_NUM = 5;                       // Number of initial uniform refinements.
const int P_INIT = 1;                             // Initial polynomial degree.
const double TAU = 0.05;                          // Time step.
const double T_FINAL = 1000;                      // Time interval length.
const int TIME_DISCR = 2;                         // 1 for implicit Euler, 2 for Crank-Nicolson.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem constants
std::complex<double> C = std::complex<double>(1./(30 * std::sqrt((double)3)), 0.0);
std::complex<double> C2 = std::complex<double>(200., 0.);

// Imaginary unit.
std::complex<double> ii = std::complex<double>(0.0, 1.0);

// Boundary markers.
const std::string BDY_BOTTOM = "1";
const std::string BDY_RIGHT = "2";
const std::string BDY_TOP = "3";
const std::string BDY_LEFT = "4";

// Initial condition for Psi.
std::complex<double> init_cond_psi(double x, double y, std::complex<double> & dx, std::complex<double> & dy)
{
  std::complex<double>  val = exp(-(x*x + y*y)/(2.*C*C)) * exp(C2 * ii * x);
  dx = (-x/(C*C)+ii*C2)*val;
  dy = (-y/(C*C))*val;
  return val;
}

// Initial condition for Phi.
std::complex<double> init_cond_phi(double x, double y, std::complex<double> & dx, std::complex<double> & dy)
{
  std::complex<double>  val = ii * C2 * exp(-(x*x + y*y)/(2.*C*C)) * exp(C2 * ii * x);
  dx = (-x/(C*C)+ii*C2)*val;
  dy = (-y/(C*C))*val;
  return val;
}

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<std::complex<double> > bc(Hermes::vector<std::string>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT), std::complex<double>(0.0, 0.0));

  EssentialBCs<std::complex<double> > bcs(&bc);

  // Create an H1 space.
  H1Space<std::complex<double> >* phi_space = new H1Space<std::complex<double> >(&mesh, &bcs, P_INIT);
  H1Space<std::complex<double> >* psi_space = new H1Space<std::complex<double> >(&mesh, &bcs, P_INIT);

  int ndof = Space<std::complex<double> >::get_num_dofs(Hermes::vector<Space<std::complex<double> > *>(phi_space, psi_space));
  info("ndof = %d.", ndof);

  // Initialize previous time level solutions.
  ConstantSolution<std::complex<double> > phi_prev_time(&mesh, init_cond_phi);
  ConstantSolution<std::complex<double> > psi_prev_time(&mesh, init_cond_psi);

  // Initialize the weak formulation.
  WeakForm<std::complex<double> > wf(2);
  wf.add_matrix_form(0, 0, callback(biform_euler_0_0));
  wf.add_matrix_form(0, 1, callback(biform_euler_0_1));
  wf.add_matrix_form(1, 0, callback(biform_euler_1_0));
  wf.add_matrix_form(1, 1, callback(biform_euler_1_1));
  wf.add_vector_form(0, callback(liform_euler_0), HERMES_ANY, &phi_prev_time);
  wf.add_vector_form(1, callback(liform_euler_1), HERMES_ANY, &psi_prev_time);

  // Initialize views.
  ScalarView view("Psi", new WinGeom(0, 0, 600, 500));
  view.fix_scale_width(80);

  // Time stepping loop:
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {

    info("Time step %d:", ts);

    info("Solving linear system.");
    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem<std::complex<double> > dp(&wf, Hermes::vector<Space<double>* *>(phi_space, psi_space), is_linear);
   
    SparseMatrix<std::complex<double> >* matrix = create_matrix<std::complex<double> >(matrix_solver_type);
    Vector<std::complex<double> >* rhs = create_vector<std::complex<double> >(matrix_solver_type);
    LinearSolver<std::complex<double> >* solver = create_linear_solver<std::complex<double> >(matrix_solver_type, matrix, rhs);

    // Assemble the stiffness matrix and right-hand side vector.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the linear system and if successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution<std::complex<double> >::vector_to_solutions(solver->get_sln_vector(), Hermes::vector<Space<std::complex<double> >*>(phi_space, psi_space), Hermes::vector<Solution<std::complex<double> > *>(&phi_prev_time, &psi_prev_time));
    else
      error ("Matrix solver failed.\n");

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Time step %d", ts);
    view.set_title(title);
    view.show(&psi_prev_time);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
