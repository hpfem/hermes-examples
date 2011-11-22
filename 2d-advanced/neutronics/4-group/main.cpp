#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "problem_data.h"

// This example solves a 4-group neutron diffusion equation in the reactor core.
// The eigenproblem is solved using power interations.
//
// The reactor neutronics is given by the following eigenproblem:
//
//  - \nabla \cdot D_g \nabla \phi_g + \Sigma_{Rg}\phi_g - \sum_{g' \neq g} \Sigma_s^{g'\to g} \phi_{g'} =
//  = \frac{\chi_g}{k_{eff}} \sum_{g'} \nu_{g'} \Sigma_{fg'}\phi_{g'}
//
// where 1/k_{eff} is eigenvalue and \phi_g, g = 1,...,4 are eigenvectors (neutron fluxes). The current problem
// is posed in a 3D cylindrical axisymmetric geometry, leading to a 2D problem with r-z as the independent spatial 
// coordinates. The corresponding diffusion operator is given by (r = x, z = y):
//
//	\nabla \cdot D \nabla \phi = \frac{1}{x} (x D \phi_x)_x  + (D \phi_y)_y 
//
// BC:
//
// Homogeneous neumann on symmetry axis,
// d \phi_g / d n = - 0.5 \phi_g   elsewhere
//
// The eigenproblem is numerically solved using common technique known as the power method (power iterations):
//
//  1) Make an initial estimate of \phi_g and k_{eff}
//  2) For n = 1, 2,...
//         solve for \phi_g using previous k_prev
//         solve for new k_{eff}
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{new}
//               k_new =  k_prev -------------------------------------------------------------------------
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{prev}
//  3) Stop iterations when
//
//     |   k_new - k_prev  |
//     | ----------------- |  < epsilon
//     |       k_new       |
//
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Initial polynomial degrees for approximation of group fluxes.
const int P_INIT_1 = 1, P_INIT_2 = 1, P_INIT_3 = 1, P_INIT_4 = 1;                           
// Tolerance for the eigenvalue.
const double ERROR_STOP = 1e-5;                   
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-8;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                  

// Initial eigenvalue approximation.
double k_eff = 1.0;         

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load(mesh_file.c_str(), &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Solution variables.
  Solution<double> sln1, sln2, sln3, sln4;
  Hermes::vector<Solution<double>*> solutions(&sln1, &sln2, &sln3, &sln4);
  
  // Define initial conditions.
  info("Setting initial conditions.");
  ConstantSolution<double> iter1(&mesh, 1.00), iter2(&mesh, 1.00), iter3(&mesh, 1.00), iter4(&mesh, 1.00);

  Hermes::vector<MeshFunction<double>*> iterates(&iter1, &iter2, &iter3, &iter4);

  // Create H1 spaces with default shapesets.
  H1Space<double> space1(&mesh, P_INIT_1);
  H1Space<double> space2(&mesh, P_INIT_2);
  H1Space<double> space3(&mesh, P_INIT_3);
  H1Space<double> space4(&mesh, P_INIT_4);
  Hermes::vector<const Space<double>* > spaces(&space1, &space2, &space3, &space4);
  int ndof = Space<double>::get_num_dofs(spaces);
  info("ndof = %d", ndof);

  // Initialize views.
  ScalarView view1("Neutron flux 1", new WinGeom(0, 0, 320, 600));
  ScalarView view2("Neutron flux 2", new WinGeom(350, 0, 320, 600));
  ScalarView view3("Neutron flux 3", new WinGeom(700, 0, 320, 600));
  ScalarView view4("Neutron flux 4", new WinGeom(1050, 0, 320, 600));
  
  // Do not show meshes, set 3D mode.
  view1.show_mesh(false); view1.set_3d_mode(true);
  view2.show_mesh(false); view2.set_3d_mode(true);
  view3.show_mesh(false); view3.set_3d_mode(true);
  view4.show_mesh(false); view4.set_3d_mode(true);
  
  // Load physical data of the problem for the 4 energy groups.
  Hermes::Hermes2D::WeakFormsNeutronics::Multigroup::MaterialProperties::Diffusion::MaterialPropertyMaps matprop(4);
  matprop.set_D(D);
  matprop.set_Sigma_r(Sr);
  matprop.set_Sigma_s(Ss);
  matprop.set_Sigma_a(Sa);
  matprop.set_Sigma_f(Sf);
  matprop.set_nu(nu);
  matprop.set_chi(chi);
  matprop.validate();
  
  // Printing table of material properties.
  std::cout << matprop;
  
  // Initialize the weak formulation.
  CustomWeakForm wf(matprop, iterates, k_eff, bdy_vacuum);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, spaces);

  // Initialize Newton solver.
  NewtonSolver<double> newton(&dp, matrix_solver);

  // Time measurement.
  TimePeriod cpu_time;
      
  // Main power iteration loop:
  int it = 1; bool done = false;
  do
  {
    info("------------ Power iteration %d:", it);
    
    info("Newton's method (matrix problem solved by %s).", MatrixSolverNames[matrix_solver].c_str());
    
    // Perform Newton's iteration.
    try
    {
      // The problem is linear and we can use NULL 
      // to make Newton start from zero vector.
      newton.solve_keep_jacobian(NULL, NEWTON_TOL, NEWTON_MAX_ITER);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    }
       
    // Debug.
    //printf("\n=================================================\n");
    //for (int d = 0; d < ndof; d++) printf("%g ", newton.get_sln_vector()[d]);

    // Translate the resulting coefficient vector into a Solution.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces, solutions);
    
    // Show intermediate solutions.
    view1.show(&sln1);    
    view2.show(&sln2);
    view3.show(&sln3);    
    view4.show(&sln4);
    
    // Compute eigenvalue.
    SourceFilter source(solutions, &matprop, core);
    SourceFilter source_prev(iterates, &matprop, core);
    
    double k_new = k_eff * (integrate(&source, core) / integrate(&source_prev, core));
    info("Largest eigenvalue: %.8g, rel. difference from previous it.: %g", k_new, fabs((k_eff - k_new) / k_new));
    
    // Stopping criterion.
    if (fabs((k_eff - k_new) / k_new) < ERROR_STOP) done = true;

    // Update eigenvalue.
    k_eff = k_new;
    wf.update_keff(k_eff);
    
    if (!done)
    {
      // Save solutions for the next iteration.
      iter1.copy(&sln1);    
      iter2.copy(&sln2);
      iter3.copy(&sln3);    
      iter4.copy(&sln4);
      
      it++;
    }
  }
  while (!done);
  
  // Time measurement.
  cpu_time.tick();
  
  // Show solutions.
  view1.show(&sln1);
  view2.show(&sln2);
  view3.show(&sln3);    
  view4.show(&sln4);
  
  // Skip visualization time.
  cpu_time.tick(HERMES_SKIP);

  // Print timing information.
  verbose("Total running time: %g s", cpu_time.accumulated());
    
  // Wait for all views to be closed.
  View::wait();
  return 0;
}
