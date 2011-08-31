#define HERMES_REPORT_ALL
#include "definitions.h"

//  This is the ninth in the series of NIST benchmarks with known exact solutions. This benchmark
//  has four different versions, use the global variable PROB_PARAM below to switch among them.
//  It differs from 09-wave-front in the mesh adaptation method.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: -Laplace u - f = 0
//
//  Known exact solution; atan(ALPHA * (sqrt(pow(x - X_LOC, 2) + pow(y - Y_LOC, 2)) - R_ZERO));
//  See the class CustomExactSolution.
//
//  Domain: unit square (0, 1) x (0, 1), see the file square.mesh.
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

int PARAM = 3;         // PARAM determines which parameter values you wish to use 
                       //            for the steepness and location of the wave front. 
                       // #| name   |   ALPHA | X_LOC	| Y_LOC | R_ZERO
                       // 0: mild		    20      -0.05  -0.05    0.7
                       // 1: steep      1000    -0.05  -0.05    0.7
                       // 2: asymmetric 1000     1.5    0.25    0.92
                       // 3: well       50       0.5    0.5     0.25

const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double ERR_STOP = 0.5;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const bool USE_RESIDUAL_ESTIMATOR = true;         // Add also the norm of the residual to the error estimate of each element.
const bool USE_ENERGY_NORM_NORMALIZATION = false; // Use energy norm for error estimate normalization and measuring of exact error.
const bool TEST_ELEMENT_BASED_KELLY = false;      // Test if two possible approaches to interface error estimators accumulation give same results.

int main(int argc, char* argv[])
{
  // Define problem parameters: (x_loc, y_loc) is the center of the circular wave front, R_ZERO is the distance from the 
  // wave front to the center of the circle, and alpha gives the steepness of the wave front.
  double alpha, x_loc, y_loc, r_zero;
  switch(PARAM) {
  case 0:
    alpha = 20;
    x_loc = -0.05;
    y_loc = -0.05;
    r_zero = 0.7;
    break;
  case 1:
    alpha = 1000;
    x_loc = -0.05;
    y_loc = -0.05;
    r_zero = 0.7;
    break;
  case 2:
    alpha = 1000;
    x_loc = 1.5;
    y_loc = 0.25;
    r_zero = 0.92;
    break;
  case 3:
    alpha = 50;
    x_loc = 0.5;
    y_loc = 0.5;
    r_zero = 0.25;
    break;
  default:   // The same as 0.
    alpha = 20;
    x_loc = -0.05;
    y_loc = -0.05;
    r_zero = 0.7;
    break;
  }

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square_quad.mesh", &mesh);     // quadrilaterals
  // mloader.load("square_tri.mesh", &mesh);   // triangles

  // Perform initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  
  // Set exact solution.
  CustomExactSolution exact(&mesh, alpha, x_loc, y_loc, r_zero);

  // Define right-hand side.
  CustomRightHandSide rhs(alpha, x_loc, y_loc, r_zero);

  // Initialize the weak formulation.
  CustomWeakForm wf(&rhs);
  // Equivalent, but slower:
  // DefaultWeakFormPoisson<double> wf(Hermes::HERMES_ANY, HERMES_ONE, &rhs);

  // Initialize boundary conditions.
  DefaultEssentialBCNonConst<double> bc("Bdy", &exact);
  EssentialBCs<double> bcs(&bc);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  
  // Initialize approximate solution.
  Solution<double> sln;

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 440, 350));
  sview.show_mesh(false);
  sview.fix_scale_width(50);
  Views::OrderView<double>  oview("Polynomial orders", new Views::WinGeom(450, 0, 420, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  Hermes::TimePeriod cpu_time;

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d (%d DOF):", as, space.get_num_dofs());
    cpu_time.tick();

    // Assemble the discrete problem.    
    DiscreteProblem<double> dp(&wf, &space);

    // Initial coefficient vector for the Newton's method.  
    int ndof = space.get_num_dofs();
    double* coeff_vec = new double[ndof];
    memset(coeff_vec, 0, ndof * sizeof(double));
    
    NewtonSolver<double> newton(&dp, matrix_solver_type);
    newton.set_verbose_output(false);

    if (!newton.solve(coeff_vec)) 
      error("Newton's iteration failed.");
    else
      Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);
    
    cpu_time.tick();
    verbose("Solution: %g s", cpu_time.last());
    
    // Calculate element errors and total error estimate.
    BasicKellyAdapt<double> adaptivity(&space);
    
    if (USE_ENERGY_NORM_NORMALIZATION)
      adaptivity.set_error_form(new EnergyErrorForm(&wf));

    if (USE_RESIDUAL_ESTIMATOR) 
      adaptivity.add_error_estimator_vol(new ResidualErrorForm(&rhs));
    
    double err_est_rel = adaptivity.calc_err_est(&sln) * 100;  
    double err_exact_rel = Global<double>::calc_rel_error(&sln, &exact, HERMES_H1_NORM) * 100;
    
    cpu_time.tick();
    verbose("Error calculation: %g s", cpu_time.last());
    
    info("err_est_rel: %g%%, err_exact_rel: %g%%", err_est_rel, err_exact_rel);

    if (TEST_ELEMENT_BASED_KELLY)
    {
      cpu_time.tick();
      
      // Ensure that the two possible approaches to interface error estimators accumulation give same results.
      KellyTypeAdapt<double> adaptivity2(&space, false);
      adaptivity2.disable_aposteriori_interface_scaling();
      adaptivity2.add_error_estimator_surf(new InterfaceErrorForm);
      
      if (USE_ENERGY_NORM_NORMALIZATION)
        adaptivity2.set_error_form(new EnergyErrorForm(&wf));
      
      if (USE_RESIDUAL_ESTIMATOR) 
      {
        adaptivity2.add_error_estimator_vol(new ResidualErrorForm(&rhs));
        adaptivity2.set_volumetric_scaling_const(1./24.);
      }
      
      double err_est_rel2 = adaptivity2.calc_err_est(&sln) * 100;  
      double err_exact_rel2 = adaptivity2.calc_err_exact(&sln, &exact, false) * 100;
      
      info("err_est_rel_2: %g%%, err_exact_rel_2: %g%%", err_est_rel2, err_exact_rel2);
      
      if (fabs(err_est_rel2 - err_est_rel) >= 1e-13 || fabs(err_exact_rel2 - err_exact_rel) >= 1e-13)
      {
        info("err_est_rel diff: %1.15g, err_exact_rel diff: %1.15g", 
        std::abs(err_est_rel2 - err_est_rel), std::abs(err_exact_rel2 - err_exact_rel));
        err_est_rel = err_exact_rel = 0; // Exit the adaptivity loop.
      }
      
      cpu_time.tick(Hermes::HERMES_SKIP);
    }
    
    // Report results.
    
    cpu_time.tick();
    double accum_time = cpu_time.accumulated();
    
    // View the approximate solution and polynomial orders.
    sview.show(&sln);
    oview.show(&space);
       
    // Add entry to DOF and CPU convergence graphs.    
    graph_dof.add_values(space.get_num_dofs(), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(accum_time, err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
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
      done = adaptivity.adapt(THRESHOLD, STRATEGY, MESH_REGULARITY);
    
    cpu_time.tick();
    verbose("Adaptation: %g s", cpu_time.last());
    
    // Increase the counter of performed adaptivity steps.
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
  }
  while (done == false);

  cpu_time.tick();
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
