#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This example shows how to set an arbitrary initial guess for the
//  Newton's method, and nonzero Dirichlet boundary conditions.
//
//  PDE: magnetostatics with nonlinear magnetic permeability
//  curl[1/mu curl u] = current_density.
//
//  Domain: unit square (-10,10)^2.
//
//  BC: Dirichlet, see function dir_lift() below.
//
//  The following parameters can be changed:

// Initial polynomial degree.
const int P_INIT = 3;                             
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-10;                  
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 1000;                 
// Number between 0 and 1 to damp Newton's iterations.
const double NEWTON_DAMPING = 1.0;                
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Problem parameters.
double MU_VACUUM = 4. * M_PI * 1e-7;
// Constant initial condition for the magnetic potential.
double INIT_COND = 0.0;                           
// Volume source term.
double CURRENT_DENSITY = 1e9;                     

// Material and boundary markers.
const std::string MAT_AIR = "e2";
const std::string MAT_IRON_1 = "e0";
const std::string MAT_IRON_2 = "e3";
const std::string MAT_COPPER = "e1";
const std::string BDY_DIRICHLET = "bdy";

int main(int argc, char* argv[])
{
  // Define nonlinear magnetic permeability via a cubic spline.
  Hermes::vector<double> mu_inv_pts(0.0,      0.5,      0.9,      1.0,      1.1,      1.2,      1.3,
                                    1.4,      1.6,      1.7,      1.8,      1.9,      3.0,      5.0,     10.0);
  Hermes::vector<double> mu_inv_val(1/1500.0, 1/1480.0, 1/1440.0, 1/1400.0, 1/1300.0, 1/1150.0, 1/950.0,
                                    1/750.0,  1/250.0,  1/180.0,  1/175.0,  1/150.0,  1/20.0,   1/10.0,  1/5.0);

  // Create the cubic spline (and plot it for visual control).
  double bc_left = 0.0;
  double bc_right = 0.0;
  bool first_der_left = false;
  bool first_der_right = false;
  bool extrapolate_der_left = false;
  bool extrapolate_der_right = false;
  CubicSpline mu_inv_iron(mu_inv_pts, mu_inv_val, bc_left, bc_right, first_der_left, first_der_right,
                          extrapolate_der_left, extrapolate_der_right);
  info("Saving cubic spline into a Pylab file spline.dat.");
  // The interval of definition of the spline will be
  // extended by "interval_extension" on both sides.
  double interval_extension = 1.0; 
  bool plot_derivative = false;
  mu_inv_iron.plot("spline.dat", interval_extension, plot_derivative);
  plot_derivative = true;
  mu_inv_iron.plot("spline_der.dat", interval_extension, plot_derivative);

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("actuator.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_essential(BDY_DIRICHLET, 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof: %d", ndof);

  // Initialize the weak formulation
  // This determines the increase of integration order
  // for the axisymmetric term containing 1/r. Default is 3.
  int order_inc = 3; 
  CustomWeakFormMagnetostatics wf(MAT_IRON_1, MAT_IRON_2, &mu_inv_iron, MAT_AIR,
      MAT_COPPER, MU_VACUUM, CURRENT_DENSITY, order_inc);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);
  
  // Initialize the solution.
  ConstantSolution<double> sln(&mesh, INIT_COND);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting to obtain initial vector for the Newton's method.");
  double* coeff_vec = new double[ndof] ;
  OGProjection<double>::project_global(&space, &sln, coeff_vec, matrix_solver);

  // Perform Newton's iteration.
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver);
  bool verbose = true;
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

  // Translate the resulting coefficient vector into the Solution sln.
  Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);

  // Cleanup.
  delete [] coeff_vec;

  // Visualise the solution and mesh.
  ScalarView s_view1("Vector potential", new WinGeom(0, 0, 350, 450));
  FilterVectorPotential vector_potential(Hermes::vector<MeshFunction<double> *>(&sln, &sln),
      Hermes::vector<int>(H2D_FN_VAL, H2D_FN_VAL));
  s_view1.show_mesh(false);
  s_view1.show(&vector_potential);

  ScalarView s_view2("Flux density", new WinGeom(360, 0, 350, 450));
  FilterFluxDensity flux_density(Hermes::vector<MeshFunction<double> *>(&sln, &sln));
  s_view2.show_mesh(false);
  s_view2.show(&flux_density);

  // Output solution in VTK format.
  Linearizer lin;
  bool mode_3D = true;
  lin.save_solution_vtk(&flux_density, "sln.vtk", "Flux density", mode_3D);
  info("Solution in VTK format saved to file %s.", "sln.vtk");

  OrderView o_view("Mesh", new WinGeom(720, 0, 350, 450));
  o_view.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

