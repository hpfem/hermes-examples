#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// This example shows how to model harmonic steady state in parallel plate waveguide.
// The complex-valued Helmholtz equation is solved by decomposing it into two equations 
// for the real and imaginary part of the E field. Two typical boundary conditions used in 
// high-frequency problems are demonstrated.
//
// PDE: Helmholtz equation for electric field
//
//    -Laplace E  - omega^2*mu*epsilon*E + j*omega*sigma*mu*E = 0
//
// BC:                     Gamma_1 (perfect)
//                   ----------------------------
//  Gamma_3 (left)   |                           |  Gamma_4 (impedance)
//                   ----------------------------
//                         Gamma_2 (perfect)
//
//     1) Perfect conductor boundary condition Ex = 0 on Gamma_1 and Gamma_2.
//     2) Essential (Dirichlet) boundary condition on Gamma_3
//          Ex(y) = E_0 * cos(y*M_PI/h), where h is height of the waveguide.
//     3) Impedance boundary condition on Gamma_4
//          dE/dn = j*beta*E.
//
// The following parameters can be changed:

// Initial polynomial degree of all elements.
const int P_INIT = 6;                                  
// Number of initial mesh refinements.
const int INIT_REF_NUM = 3;                            
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;       

// Newton's method.
const double NEWTON_TOL = 1e-8;
const int NEWTON_MAX_ITER = 100;

// Problem parameters.
// Relative permittivity.
const double epsr = 1.0;                    
// Permittivity of vacuum F/m.
const double eps0 = 8.85418782e-12;         
const double eps = epsr * eps0;
// Relative permeablity.
const double mur = 1.0;                     
// Permeability of vacuum H/m.
const double mu0 = 4*M_PI*1e-7;             
const double mu = mur * mu0;
// Frequency MHz.
const double frequency = 3e9;               
// Angular velocity.
const double omega = 2*M_PI * frequency;    
// Conductivity Ohm/m.
const double sigma = 0;                     
// Propagation constant.
const double beta = 54;                     
// Input electric intensity.
const double E0 = 100;                      
// Height of waveguide.
const double h = 0.1;                       

int main(int argc, char* argv[])
{
    // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform uniform mesh refinement.
  // 2 is for vertical split.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(2); 

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc1("Bdy_perfect", 0.0);
  EssentialBCNonConst bc2("Bdy_left");
  DefaultEssentialBCConst<double> bc3("Bdy_left", 0.0);
  EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc1, &bc2));
  EssentialBCs<double> bcs_im(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc1, &bc3));

  // Create an H1 space with default shapeset.
  H1Space<double> e_r_space(&mesh, &bcs, P_INIT);
  H1Space<double> e_i_space(&mesh, &bcs_im, P_INIT);
  int ndof = Space<double>::get_num_dofs(&e_r_space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakFormHelmholtz wf(eps, mu, omega, sigma, beta, E0, h);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<const Space<double>*>(&e_r_space, &e_i_space));

  // Initialize the solutions.
  Solution<double> e_r_sln, e_i_sln;

  // Initial coefficient vector for the Newton's method.  
  ndof = Space<double>::get_num_dofs(Hermes::vector<const Space<double>*>(&e_r_space, &e_i_space));

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

  // Translate the resulting coefficient vector into Solutions.
  Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<const Space<double>*>(&e_r_space, &e_i_space), 
      Hermes::vector<Solution<double>*>(&e_r_sln, &e_i_sln));

  // Visualize the solution.
  ScalarView viewEr("Er [V/m]", new WinGeom(0, 0, 800, 400));
  viewEr.show(&e_r_sln);
  // viewEr.save_screenshot("real_part.bmp");

  ScalarView viewEi("Ei [V/m]", new WinGeom(0, 450, 800, 400));
  viewEi.show(&e_i_sln);
  // viewEi.save_screenshot("imaginary_part.bmp");

  // Wait for the view to be closed.
  View::wait();

  return 0;
}
