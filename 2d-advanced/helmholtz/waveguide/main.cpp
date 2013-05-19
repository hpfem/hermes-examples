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
const int P_INIT = 4;                                  
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
double frequency = 16e9;
// Angular velocity.
double omega = 2*M_PI * frequency;    
// Conductivity Ohm/m.
const double sigma = 0;                     
// Propagation constant.
const double beta = 54;                     
// Input electric intensity.
const double E0 = 100;                      
// Height of waveguide.
const double h = 1.0;

int main(int argc, char* argv[])
{
    // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("initial.xml", mesh);

  mesh->refine_all_elements();

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc1(Hermes::vector<std::string>("2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 0.0);
  EssentialBCNonConst bc2("1");
  DefaultEssentialBCConst<double> bc3("1", 0.0);
  EssentialBCs<double> bcs(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc1, &bc2));
  EssentialBCs<double> bcs_im(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc1, &bc3));

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> e_r_space(new H1Space<double> (mesh, &bcs, P_INIT));
  SpaceSharedPtr<double> e_i_space(new H1Space<double>(mesh, &bcs_im, P_INIT));
  int ndof = Space<double>::get_num_dofs(e_r_space);
  
  // Initialize the weak formulation.
  WeakFormHelmholtz wf(eps, mu, omega, sigma, beta, E0, h);
  wf.set_verbose_output(false);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<SpaceSharedPtr<double> >(e_r_space, e_i_space));

  // Initialize the solutions.
  MeshFunctionSharedPtr<double> e_r_sln(new Solution<double>), e_i_sln(new Solution<double>);

  // Initial coefficient vector for the Newton's method.  
  ndof = Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(e_r_space, e_i_space));

  Hermes::Hermes2D::NewtonSolver<double> newton(&dp);
  newton.set_tolerance(NEWTON_TOL);
  newton.set_max_allowed_iterations(NEWTON_MAX_ITER);
  
  ScalarView viewEr("Er [V/m]", new WinGeom(600, 0, 700, 200));
  viewEr.show_mesh(false);
  viewEr.set_3d_mode(true);
  ScalarView viewEi("Ei [V/m]", new WinGeom(600, 220, 700, 200));
  
  ScalarView viewMagnitude("Magnitude of E [V/m]", new WinGeom(600, 440, 700, 200));
  viewMagnitude.show_mesh(false);
  viewMagnitude.show_contours(50., 0.);
  
  while(true)
  {
    std::cout << "Frequency: " << frequency << " Hz" << std::endl;
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    };

    // Translate the resulting coefficient vector into Solutions.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<SpaceSharedPtr<double> >(e_r_space, e_i_space), 
        Hermes::vector<MeshFunctionSharedPtr<double> >(e_r_sln, e_i_sln));

    // Visualize the solution.
    viewEr.show(e_r_sln, 0.1);
    // viewEr.save_screenshot("real_part.bmp");

    viewEi.show(e_i_sln, 0.1);
    // viewEi.save_screenshot("imaginary_part.bmp");

    MeshFunctionSharedPtr<double> magnitude(new MagFilter<double>(Hermes::vector<MeshFunctionSharedPtr<double> >(e_r_sln, e_i_sln)));
    viewMagnitude.show(magnitude, HERMES_EPS_LOW);

    char* change_state = new char[1000];
    std::cout << "Done?";
    std::cin.getline(change_state, 1);
    if(!strcmp(change_state, "y"))
      break;
    std::cout << std::endl;
    std::cout << "Frequency change [1e9 Hz]: ";
    double frequency_change;
    std::cin >> frequency_change;
    
    frequency += 1e9 * frequency_change;
    omega = 2*M_PI * frequency;  
    wf.set_parameters(eps, mu, omega, sigma, beta, E0, h);
    newton.set_weak_formulation(&wf);
  }

  // Wait for the view to be closed.
  View::wait();

  return 0;
}
