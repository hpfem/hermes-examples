#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// This example shows how to model harmonic steady state in parallel plate waveguide.
// The Helmholtz equation is solved and there are demonstrated two typical boundary
// conditions used in high frequency domain.
//
// PDE: Helmholtz equation for electric field
//
//    Delta E  + (omega^2*mu*epsilon - j*omega*sigma*mu)*E = 0
//
// BC:              Gamma_1
//             ----------------------------
//  Gamma_3    |                           |  Gamma_4
//             ----------------------------
//                  Gamma_2
//
//     1) Dirichlet boundary condition Ex = 0 (perfect eletric conductor) on Gamma_1 and Gamma_2.
//     2) Essential (Dirichlet) boundary condition on Gamma_3
//          Ex(y) = E_0 * cos(y*M_PI/h), where h is height of the waveguide ()
//     3) Newton boundary condition (impedance matching) on Gamma_4
//          dE/dn = j*beta*E
//
// The following parameters can be changed:

// Initial polynomial degree of all elements.
const int P_INIT = 6;                                  
// Number of initial mesh refinements.
const int INIT_REF_NUM = 3;                            
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;       

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
  EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc1, &bc2));

  // Create an H1 space with default shapeset.
  H1Space<double> e_r_space(&mesh, &bcs, P_INIT);
  H1Space<double> e_i_space(&mesh, &bcs, P_INIT);
  int ndof = Space<double>::get_num_dofs(&e_r_space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakFormHelmholtz wf(eps, mu, omega, sigma, beta, E0, h);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<Space<double>*>(&e_r_space, &e_i_space));

  // Initialize the solutions.
  Solution<double> e_r_sln, e_i_sln;

  // Initial coefficient vector for the Newton's method.  
  ndof = Space<double>::get_num_dofs(Hermes::vector<Space<double>*>(&e_r_space, &e_i_space));
  double* coeff_vec = new double[ndof];
  memset(coeff_vec, 0, ndof * sizeof(double));

  Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver_type);
  try
  {
    newton.solve(coeff_vec);
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.printMsg();
    error("Newton's iteration failed.");
  };

  // Translate the resulting coefficient vector into the Solution<double> sln.
  Solution<double>::vector_to_solutions(coeff_vec, Hermes::vector<Space<double>*>(&e_r_space, &e_i_space), Hermes::vector<Solution<double>*>(&e_r_sln, &e_i_sln));

  // Visualize the solution.
  ScalarView viewEr("Er [V/m]", new WinGeom(0, 0, 800, 400));
  viewEr.show(&e_r_sln);
  // viewEr.save_screenshot("real_part.bmp");

  ScalarView viewEi("Ei [V/m]", new WinGeom(0, 450, 800, 400));
  viewEi.show(&e_i_sln);
  // viewEi.save_screenshot("imaginary_part.bmp");

  // Wait for the view to be closed.
  View::wait();

  // Clean up.
  delete [] coeff_vec;

  return 0;
}
