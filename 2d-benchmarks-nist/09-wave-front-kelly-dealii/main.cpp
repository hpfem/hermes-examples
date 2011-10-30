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

// PARAM determines which parameter values you wish to use 
//            for the steepness and location of the wave front. 
// #| name   |   ALPHA | X_LOC	| Y_LOC | R_ZERO
// 0: mild		    20      -0.05  -0.05    0.7
// 1: steep      1000    -0.05  -0.05    0.7
// 2: asymmetric 1000     1.5    0.25    0.92
// 3: well       50       0.5    0.5     0.25
int PARAM = 3;         

// Initial polynomial degree of mesh elements.
const int P_INIT = 2;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;                       
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// Stopping criterion for adaptivity (rel. error tolerance between the
// reference mesh and coarse mesh solution in percent).
const double ERR_STOP = 0.1;                      
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  

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
  default:   
    // The same as 0.
    alpha = 20;
    x_loc = -0.05;
    y_loc = -0.05;
    r_zero = 0.7;
    break;
  }
  
  // Set adaptivity options according to the command line argument.
  int refinement_mode;
  if (argc != 2)
    refinement_mode = 0;
  else
  {
    refinement_mode = atoi(argv[1]);
    if (refinement_mode < 1 || refinement_mode > 12)
      error("Invalid run case: %d (valid range is [1,12])", refinement_mode);
  }
  
  double threshold = 0.3;                     
  // strategy = 0 ... Refine elements until sqrt(threshold) times total error is processed. 
  //                  If more elements have similar errors, refine all to keep the mesh symmetric.
  int strategy = 0;       
  // strategy = 1 ... Refine all elements whose error is larger than threshold times max. element error.
  // Add also the norm of the residual to the error estimate of each element.
  bool use_residual_estimator = false;        
  // Use energy norm for error estimate normalization and measuring of exact error.
  bool use_energy_norm_normalization = false; 
  
  switch (refinement_mode)
  {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    {
      strategy = 0;
      break;
    }
    
    case 7 :
    case 8 :
    case 9 :
    case 10:
    case 11:
    case 12:
    {
      strategy = 1;
      break;
    }
  }
  
  switch (refinement_mode)
  {
    case 1:
    case 2:
    case 3:
    case 7:
    case 8:
    case 9:
    {
      threshold = 0.3;
      break;
    }
    
    case 4:
    case 5:
    case 6:
    case 10:
    case 11:
    case 12:
    {
      threshold = 0.3*0.3;
      break;
    }
  }
  
  switch (refinement_mode)
  {
    case 2:
    case 3:
    case 5:
    case 6:
    case 8:
    case 9:
    case 11:
    case 12:
    {
      use_residual_estimator = true;
      break;
    }
  }
  
  switch (refinement_mode)
  {
    case 3:
    case 6:
    case 9:
    case 12:
    {
      use_energy_norm_normalization = true;
      break;
    }
  }
  
  double setup_time = 0, assemble_time = 0, solve_time = 0, adapt_time = 0;
  
  // Time measurement.
  Hermes::TimePeriod wall_clock;
  // Stop counting time for adaptation.
  wall_clock.tick(); 

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square_quad.mesh", &mesh);    

  // Perform initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  
  // Stop counting time for adaptation.
  wall_clock.tick(); 
  adapt_time += wall_clock.last();
  
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
  Views::ScalarView sview("Solution", new Views::WinGeom(800, 0, 400, 400));
  sview.show_mesh(false);
  sview.set_3d_mode();
  sview.set_palette(Views::H2DV_PT_HUESCALE);
  sview.fix_scale_width(50);
  Views::OrderView  oview("Mesh", new Views::WinGeom(0, 0, 800, 800));
  oview.set_palette(Views::H2DV_PT_INVGRAYSCALE);

  // DOF and CPU convergence graphs.
  ConvergenceTable conv_table;
  conv_table.add_column(" cycle ", "%6d ");
  conv_table.add_column("  H1 error ", " %.4e ");
  conv_table.add_column("   ndof  ", " %6d ");
  conv_table.add_column(" total time ", " %8.3f  ");
  conv_table.add_column(" setup time ", "  %8.3f  ");
  conv_table.add_column(" assem time ", "  %8.3f  ");
  conv_table.add_column(" solve time ", "  %8.3f  ");
  conv_table.add_column(" adapt time ", "  %8.3f  ");
  
  wall_clock.tick(Hermes::HERMES_SKIP);

  // Adaptivity loop:
  int as = 0; bool done = false;
  do
  {
    // Start counting setup time.
    wall_clock.tick(); 
    
    // Assemble the discrete problem.    
    DiscreteProblem<double> dp(&wf, &space);

    // Actual ndof.
    int ndof = space.get_num_dofs();
    
    NewtonSolver<double> newton(&dp, matrix_solver);
    newton.set_verbose_output(false);
    newton.attach_timer(&wall_clock);

    // Setup time continues in NewtonSolver::solve().
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    };

    setup_time += newton.get_setup_time();
    assemble_time += newton.get_assemble_time();
    solve_time += newton.get_solve_time();
     
    // Start counting time for adaptation.
    wall_clock.tick();  
    Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);
    
    double err_exact = Global<double>::calc_abs_error(&sln, &exact, HERMES_H1_NORM);
   
    // Report results.
    info(" Cycle %d:", as);
    info("    Number of degrees of freedom: %d", ndof);
    info("    H1 error w.r.t. exact soln.:  %g", err_exact);
    
    // Stop counting time for adaptation.
    wall_clock.tick(); 
    double accum_time = wall_clock.accumulated();
    adapt_time += wall_clock.last();
    
    // View the approximate solution and polynomial orders.
    //sview.show(&sln);
    //oview.show(&space);
    //Views::View::wait(Views::HERMES_WAIT_KEYPRESS);
    
    conv_table.add_value(0, as);
    conv_table.add_value(1, err_exact);
    conv_table.add_value(2, ndof);
    conv_table.add_value(3, accum_time);
    conv_table.add_value(4, setup_time);
    conv_table.add_value(5, assemble_time);
    conv_table.add_value(6, solve_time);
    conv_table.add_value(7, adapt_time);
    
    // Start counting time for adaptation.
    wall_clock.tick(); 
    
    if (err_exact < ERR_STOP) 
      done = true;
    else
    {
      // Calculate element errors and total error estimate.
      Hermes::Hermes2D::BasicKellyAdapt<double> adaptivity(&space);
      
      unsigned int error_flags = HERMES_TOTAL_ERROR_ABS;
      
      if (use_energy_norm_normalization)
      {
        error_flags |= HERMES_ELEMENT_ERROR_REL;
        adaptivity.set_error_form(new EnergyErrorForm(&wf));
      }
      else
        error_flags |= HERMES_ELEMENT_ERROR_ABS;
      
      if (use_residual_estimator) 
        adaptivity.add_error_estimator_vol(new ResidualErrorForm(&rhs));
      
      double err_est_rel = adaptivity.calc_err_est(&sln, error_flags);  
      
      done = adaptivity.adapt(threshold, strategy, MESH_REGULARITY);
      
      // Stop counting time for adaptation.
      wall_clock.tick(); 
      adapt_time += wall_clock.last();
    }
        
    // Increase the counter of performed adaptivity steps.
    if (done == false)  
      as++;
    else
    {
      info("Total running time:  %g s", wall_clock.accumulated());
      info("   Setup:            %g s", setup_time);
      info("   Assemble:         %g s", assemble_time);
      info("   Solve:            %g s", solve_time);
      info("   Adapt:            %g s", adapt_time);
      
      //sview.show(&sln);
      oview.show(&space);
      oview.save_screenshot(("final_mesh-"+itos(refinement_mode)+".bmp").c_str());
      oview.close();
      conv_table.save(("conv_table-"+itos(refinement_mode)+".dat").c_str());
    }
  }
  while (done == false);

  // Wait for all views to be closed.
  //Views::View::wait();
  return 0;
}
