#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;

//  This is the second in the series of NIST benchmarks with known exact solutions. This benchmark
//  has four different versions, use the global variable PARAM below to switch among them.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: -Laplace u = 0
//
//  Known exact solution: (pow(sqrt(x*x + y*y), alpha) * sin(alpha * atan2(y,x))).
//  See functions fn() and fndd() in "exact_solution.cpp".
//
//  Domain: square (-1, 1)^2, with a section removed from the clockwise side of the positive x-axis.
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

// PARAM determines which parameter values you wish to use for the strength of the singularity in
// the current (nist-2) Reentrant Corner problem.
// PARAM      strength         omega            alpha
// 0:            1             5*Pi/4           4/5
// 1:            2             3*Pi/2           2/3
// 2:            3             7*Pi/4           4/7
// 3:            4             2*Pi             1/2
int PARAM = 1;

// Initial polynomial degree of mesh elements.
const int P_INIT = 3;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;                       
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;                     
// Adaptive strategy:
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
const CandList CAND_LIST = H2D_HP_ANISO_H;          
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
// Stopping criterion for adaptivity (rel. error tolerance between the
// reference mesh and coarse mesh solution in percent).
const double ERR_STOP = 0.1;                      
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;  

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;

  double alpha = 0, omega = 0;
  switch (PARAM) 
  {
    case 0: 
    mloader.load("geom0.mesh", &mesh); 
    omega = ((5.0 * M_PI)/ 4.0);
    alpha = (M_PI/ omega);
    break;
    case 1: 
    mloader.load("geom1.mesh", &mesh); 
    omega = ((3.0 * M_PI)/ 2.0);
    alpha = (M_PI/ omega);
    break;
    case 2: 
    mloader.load("geom2.mesh", &mesh); 
    omega = ((7.0 * M_PI)/ 4.0);
    alpha = (M_PI/ omega);
    break;
    case 3: mloader.load("geom3.mesh", &mesh); 
    omega = (2.0 * M_PI);
    alpha = (M_PI/omega);
    break;
    default: error("Admissible values of PARAM are 0, 1, 2, 3.");
  }

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Set exact solution.
  CustomExactSolution exact_sln(&mesh, alpha);

  // Initialize weak formulation.
  Hermes1DFunction<double> lambda(1.0);
  WeakFormsH1::DefaultWeakFormLaplace<double> wf(HERMES_ANY, &lambda);

  // Initialize boundary conditions
  DefaultEssentialBCNonConst<double> bc_essential("Bdy", &exact_sln);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);

  // Initialize approximate solution.
  Solution<double> sln;
  
  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 440, 350));
  sview.show_mesh(false);
  sview.fix_scale_width(50);
  Views::OrderView  oview("Polynomial orders", new Views::WinGeom(450, 0, 420, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  Hermes::TimePeriod cpu_time;

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    cpu_time.tick();

    // Construct globally refined reference mesh and setup reference space.
    Space<double>* ref_space = Space<double>::construct_refined_space(&space, 1);
    int ndof_ref = ref_space->get_num_dofs();

    info("---- Adaptivity step %d (%d DOF):", as, ndof_ref);
    cpu_time.tick();
    
    info("Solving on reference mesh.");
    
    // Assemble the discrete problem.    
    DiscreteProblem<double> dp(&wf, ref_space);
    
    // Initial coefficient vector for the Newton's method.  
    double* coeff_vec = new double[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(double));
    
    NewtonSolver<double> newton(&dp, matrix_solver_type);
    newton.set_verbose_output(false);
    
    Solution<double> ref_sln;
    try
    {
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    };

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);
    
    cpu_time.tick();
    verbose("Solution: %g s", cpu_time.last());
    
    // Project the fine mesh solution onto the coarse mesh.
    info("Calculating error estimate and exact error.");
    OGProjection<double>::project_global(&space, &ref_sln, &sln, matrix_solver_type);

    // Calculate element errors and total error estimate.
    Adapt<double> adaptivity(&space);
    double err_est_rel = adaptivity.calc_err_est(&sln, &ref_sln) * 100;

    // Calculate exact error.
    double err_exact_rel = Global<double>::calc_rel_error(&sln, &exact_sln, HERMES_H1_NORM) * 100;

    cpu_time.tick();
    verbose("Error calculation: %g s", cpu_time.last());
    
    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d", space.get_num_dofs(), ref_space->get_num_dofs());
    info("err_est_rel: %g%%, err_exact_rel: %g%%", err_est_rel, err_exact_rel);

    // Time measurement.
    cpu_time.tick();
    double accum_time = cpu_time.accumulated();
    
    // View the coarse mesh solution and polynomial orders.
    sview.show(&sln);
    oview.show(&space);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(space.get_num_dofs(), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(accum_time, err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");
    graph_dof_exact.add_values(space.get_num_dofs(), err_exact_rel);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(accum_time, err_exact_rel);
    graph_cpu_exact.save("conv_cpu_exact.dat");
    
    cpu_time.tick(Hermes::HERMES_SKIP);

    // If err_est too large, adapt the mesh. The NDOF test must be here, so that the solution may be visualized
    // after ending due to this criterion.
    if (err_exact_rel < ERR_STOP || space.get_num_dofs() >= NDOF_STOP) 
      done = true;
    else
      done = adaptivity.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
   
    cpu_time.tick();
    verbose("Adaptation: %g s", cpu_time.last());
    
    // Increase the counter of adaptivity steps.
    if (done == false)  
      as++;
/*    else 
    {
      sview.show(&sln);
      oview.show(&space);
      info("err_est_rel: %g%%, err_exact_rel: %g%%", err_est_rel, err_exact_rel);
    }
*/       
    // Clean up.
    delete [] coeff_vec;
    
    if(done == false) 
      delete ref_space->get_mesh();
    delete ref_space;
  }
  while (done == false);
  
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
