#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// This example solves adaptively the electric field in a simplified microwave oven.
// The waves are generated using a harmonic surface current on the right-most edge.
// (Such small cavity is present in every microwave oven). There is a circular
// load located in the middle of the main cavity, defined through a different
// permittivity -- see function in_load(...). One can either use a mesh that is
// aligned to the load via curvilinear elements (ALIGN_MESH = true), or an unaligned
// mesh (ALIGN_MESH = false). Convergence graphs are saved both wrt. the dof number
// and cpu time.
//
// PDE: time-harmonic Maxwell's equations;
//      there is circular load in the middle of the large cavity, whose permittivity
//      is different from the rest of the domain.
//
// Domain: square cavity with another small square cavity attached from outside
//         on the right.
//
// Meshes: you can either use "oven_load_circle.mesh" containing curved elements
//         aligned with the circular load, or "oven_load_square.mesh" which is not
//         aligned.
//
// BC: perfect conductor on the boundary except for the right-most edge of the small
//     cavity, where a harmonic surface current is prescribed
//
// The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Initial polynomial degree. NOTE: The meaning is different from
// standard continuous elements in the space H1. Here, P_INIT refers
// to the maximum poly order of the tangential component, and polynomials
// of degree P_INIT + 1 are present in element interiors. P_INIT = 0
// is for Whitney elements.
const int P_INIT = 2;                             
// if ALIGN_MESH == true, curvilinear elements aligned with the
// circular load are used, otherwise one uses a non-aligned mesh.
const bool ALIGN_MESH = true;                     
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;                     
// Adapt<std::complex<double> >ive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
const int STRATEGY = 0;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;          
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0. 
const double CONV_EXP = 1.0;                      
// Stopping criterion for adaptivity.
const double ERR_STOP = 1.0;                      
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Newton's method.
const double NEWTON_TOL = 1e-6;
const int NEWTON_MAX_ITER = 100;

// Problem parameters.
const double e_0 = 8.8541878176 * 1e-12;
const double mu_0 = 1.256 * 1e-6;
const double e_r = 1.0;
const double mu_r = 1.0;
const double rho = 3820.0;
const double Cp = 7.531000;
const double freq = 1.0*2450000000.0;
const double omega = 2 * M_PI * freq;
const double c = 1 / std::sqrt(e_0 * mu_0);
const double kappa  = 2 * M_PI * freq * std::sqrt(e_0 * mu_0);
const double J = 0.0000033333;

//  Boundary markers.
const std::string BDY_PERFECT_CONDUCTOR = "b2";
const std::string BDY_CURRENT = "b1";

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  if (ALIGN_MESH) 
    mloader.load("oven_load_circle.mesh", &mesh);
  else 
    mloader.load("oven_load_square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions
  DefaultEssentialBCConst<std::complex<double> > bc_essential(BDY_PERFECT_CONDUCTOR, std::complex<double>(0.0, 0.0));

  EssentialBCs<std::complex<double> > bcs(&bc_essential);

  // Create an Hcurl space with default shapeset.
  HcurlSpace<std::complex<double> > space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakForm wf(e_0, mu_0, mu_r, kappa, omega, J, ALIGN_MESH, &mesh, BDY_CURRENT);

  // Initialize coarse and reference mesh solution.
  Solution<std::complex<double> > sln, ref_sln;

  // Initialize refinements selector.
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  VectorView eview("Electric field", new WinGeom(0, 0, 580, 400));
  OrderView  oview("Polynomial orders", new WinGeom(590, 0, 550, 400));
  
  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space<std::complex<double> >* ref_space = Space<std::complex<double> >::construct_refined_space(&space);
    int ndof_ref = Space<std::complex<double> >::get_num_dofs(ref_space);

    // Initialize reference problem.
    info("Solving on reference mesh.");
    DiscreteProblem<std::complex<double> > dp(&wf, ref_space);

    // Time measurement.
    cpu_time.tick();

    // Perform Newton's iteration.
    Hermes::Hermes2D::NewtonSolver<std::complex<double> > newton(&dp, matrix_solver);
    try
    {
      newton.solve(NULL, NEWTON_TOL, NEWTON_MAX_ITER);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    };
    // Translate the resulting coefficient vector into the Solution<std::complex<double> > sln.
    Hermes::Hermes2D::Solution<std::complex<double> >::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);
  
    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection<std::complex<double> >::project_global(&space, &ref_sln, &sln, matrix_solver); 
   
    // View the coarse mesh solution and polynomial orders.
    ComplexAbsFilter magn(&sln);
    char title[100];
    sprintf(title, "Electric field, adaptivity step %d", as);
    eview.set_title(title);
    eview.set_min_max_range(0.0, 4e3);
    eview.show(&magn);
    sprintf(title, "Polynomial orders, adaptivity step %d", as);
    oview.set_title(title);
    oview.show(&space);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt<std::complex<double> >* adaptivity = new Adapt<std::complex<double> >(&space);

    // Set custom error form and calculate error estimate.
    CustomErrorForm cef(kappa);
    adaptivity->set_error_form(0, 0, &cef);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
      Space<std::complex<double> >::get_num_dofs(&space), 
      Space<std::complex<double> >::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space<std::complex<double> >::get_num_dofs(&space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (space.get_num_dofs() >= NDOF_STOP) done = true;

    delete adaptivity;
    if(done == false) delete ref_space->get_mesh();
    delete ref_space;
    
    // Increase counter.
    as++;
  }
  while (done == false);
  
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  eview.set_title("Fine mesh solution - magnitude");
  RealFilter ref_magn(&ref_sln);
  eview.show(&ref_magn);

  // Output solution in VTK format.
  Linearizer lin;
  bool mode_3D = true;
  lin.save_solution_vtk(&ref_magn, "sln.vtk", "Magnitude of E", mode_3D);
  info("Solution in VTK format saved to file %s.", "sln.vtk");

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

