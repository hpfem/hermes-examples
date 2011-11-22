#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This example uses adaptivity with dynamical meshes to solve
//  the time-dependent Richard's equation. The time discretization 
//  is backward Euler or Crank-Nicolson, and the nonlinear solver 
//  in each time step is either Newton or Picard. 
//
//  PDE: C(h)dh/dt - div(K(h)grad(h)) - (dK/dh)*(dh/dy) = 0
//  where K(h) = K_S*exp(alpha*h)                          for h < 0,
//        K(h) = K_S                                       for h >= 0,
//        C(h) = alpha*(theta_s - theta_r)*exp(alpha*h)    for h < 0,
//        C(h) = alpha*(theta_s - theta_r)                 for h >= 0.
//
//  Picard's linearization: C(h^k)dh^{k+1}/dt - div(K(h^k)grad(h^{k+1})) - (dK/dh(h^k))*(dh^{k+1}/dy) = 0
//                          Note: the index 'k' does not refer to time stepping.
//  Newton's method is more involved, see the file definitions.cpp.
//
//  Domain: rectangle (0, 8) x (0, 6.5).
//  Units: length: cm
//         time: days
//
//  BC: Dirichlet, given by the initial condition.
//
//  The following parameters can be changed:

// Choose full domain or half domain.
// const char* mesh_file = "domain-full.mesh";
std::string mesh_file = "domain-half.mesh";

// Methods.
// 1 = Newton, 2 = Picard.
const int ITERATIVE_METHOD = 1;		          
// 1 = implicit Euler, 2 = Crank-Nicolson.
const int TIME_INTEGRATION = 1;                   

// Adaptive time stepping.
// Time step (in days).
double time_step = 0.5;                           
// Timestep decrease ratio after unsuccessful nonlinear solve.
double time_step_dec = 0.5;                       
// Timestep increase ratio after successful nonlinear solve.
double time_step_inc = 1.1;                       
// Computation will stop if time step drops below this value. 
double time_step_min = 1e-8; 			                
// Maximal time step.
double time_step_max = 1.0;                       

// Elements orders and initial refinements.
// Initial polynomial degree of all mesh elements.
const int P_INIT = 1;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Number of initial mesh refinements towards the top edge.
const int INIT_REF_NUM_BDY_TOP = 0;               

// Adaptivity.
// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 1;                         
// 1... mesh reset to basemesh and poly degrees to P_INIT.   
// 2... one ref. layer shaved off, poly degrees reset to P_INIT.
// 3... one ref. layer shaved off, poly degrees decreased by one. 
const int UNREF_METHOD = 3;                       
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

// Constitutive relations.
enum CONSTITUTIVE_RELATIONS {
    CONSTITUTIVE_GENUCHTEN,    // Van Genuchten.
    CONSTITUTIVE_GARDNER       // Gardner.
};
// Use van Genuchten's constitutive relations, or Gardner's.
CONSTITUTIVE_RELATIONS constitutive_relations_type = CONSTITUTIVE_GENUCHTEN;

// Newton's and Picard's methods.
// Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL = 1e-5;                   
// Maximum allowed number of Newton iterations.
int NEWTON_MAX_ITER = 10;                         
// Stopping criterion for Picard on fine mesh.
const double PICARD_TOL = 1e-2;                   
// Maximum allowed number of Picard iterations.
int PICARD_MAX_ITER = 23;                         

// Times.
// Start-up time for time-dependent Dirichlet boundary condition.
const double STARTUP_TIME = 5.0;                  
// Time interval length.
const double T_FINAL = 1000.0;                    
// Time interval of the top layer infiltration.
const double PULSE_END_TIME = 1000.0;             
// Global time variable.
double current_time = time_step;                          

// Problem parameters.
// Initial pressure head.
double H_INIT = -15.0;                            
// Top constant pressure head -- an infiltration experiment.
double H_ELEVATION = 10.0;                        
double K_S_vals[4] = {350.2, 712.8, 1.68, 18.64}; 
double ALPHA_vals[4] = {0.01, 1.0, 0.01, 0.01};
double N_vals[4] = {2.5, 2.0, 1.23, 2.5};
double M_vals[4] = {0.864, 0.626, 0.187, 0.864};

double THETA_R_vals[4] = {0.064, 0.0, 0.089, 0.064};
double THETA_S_vals[4] = {0.14, 0.43, 0.43, 0.24};
double STORATIVITY_vals[4] = {0.1, 0.1, 0.1, 0.1};

// Precalculation of constitutive tables.
const int MATERIAL_COUNT = 4;
// 0 - constitutive functions are evaluated directly (slow performance).
// 1 - constitutive functions are linearly approximated on interval 
//     <TABLE_LIMIT; LOW_LIMIT> (very efficient CPU utilization less 
//     efficient memory consumption (depending on TABLE_PRECISION)).
// 2 - constitutive functions are aproximated by quintic splines.
const int CONSTITUTIVE_TABLE_METHOD = 2;
						  
/* Use only if CONSTITUTIVE_TABLE_METHOD == 2 */					  
// Number of intervals.        
const int NUM_OF_INTERVALS = 16;                                
// Low limits of intervals approximated by quintic splines.
double INTERVALS_4_APPROX[16] = 
      {-1.0, -2.0, -3.0, -4.0, -5.0, -8.0, -10.0, -12.0, 
      -15.0, -20.0, -30.0, -50.0, -75.0, -100.0,-300.0, -1000.0}; 
// This array contains for each integer of h function appropriate polynomial ID.
                      
// First DIM is the interval ID, second DIM is the material ID, 
// third DIM is the derivative degree, fourth DIM are the coefficients.

/* END OF Use only if CONSTITUTIVE_TABLE_METHOD == 2 */					  

/* Use only if CONSTITUTIVE_TABLE_METHOD == 1 */
// Limit of precalculated functions (should be always negative value lower 
// then the lowest expect value of the solution (consider DMP!!)
double TABLE_LIMIT = -1000.0; 		          
// Precision of precalculated table use 1.0, 0,1, 0.01, etc.....
const double TABLE_PRECISION = 0.1;               

bool CONSTITUTIVE_TABLES_READY = false;
// Polynomial approximation of the K(h) function close to saturation.
// This function has singularity in its second derivative.
// First dimension is material ID
// Second dimension is the polynomial derivative.
// Third dimension are the polynomial's coefficients.
double*** POLYNOMIALS;                            	  
// Lower bound of K(h) function approximated by polynomials.						  
const double LOW_LIMIT = -1.0;                    
const int NUM_OF_INSIDE_PTS = 0;
/* END OF Use only if CONSTITUTIVE_TABLE_METHOD == 1 */

// Boundary markers.
const std::string BDY_TOP = "1";
const std::string BDY_RIGHT = "2";
const std::string BDY_BOTTOM = "3";
const std::string BDY_LEFT = "4";

// Main function.
int main(int argc, char* argv[])
{
  ConstitutiveRelationsGenuchtenWithLayer constitutive_relations(CONSTITUTIVE_TABLE_METHOD, NUM_OF_INSIDE_PTS, LOW_LIMIT, TABLE_PRECISION, TABLE_LIMIT, K_S_vals, ALPHA_vals, N_vals, M_vals, THETA_R_vals, THETA_S_vals, STORATIVITY_vals);

  // Either use exact constitutive relations (slow) (method 0) or precalculate 
  // their linear approximations (faster) (method 1) or
  // precalculate their quintic polynomial approximations (method 2) -- managed by 
  // the following loop "Initializing polynomial approximation".
  if (CONSTITUTIVE_TABLE_METHOD == 1)
    constitutive_relations.constitutive_tables_ready = get_constitutive_tables(1, &constitutive_relations, MATERIAL_COUNT);  // 1 stands for the Newton's method.


  // The van Genuchten + Mualem K(h) function is approximated by polynomials close 
  // to zero in case of CONSTITUTIVE_TABLE_METHOD==1.
  // In case of CONSTITUTIVE_TABLE_METHOD==2, all constitutive functions are approximated by polynomials.
  info("Initializing polynomial approximations.");
  for (int i=0; i < MATERIAL_COUNT; i++)
  {
    // Points to be used for polynomial approximation of K(h).
    double* points = new double[NUM_OF_INSIDE_PTS];

    init_polynomials(6 + NUM_OF_INSIDE_PTS, LOW_LIMIT, points, NUM_OF_INSIDE_PTS, i, &constitutive_relations, MATERIAL_COUNT, NUM_OF_INTERVALS, INTERVALS_4_APPROX);
  }
  
  constitutive_relations.polynomials_ready = true;
  if (CONSTITUTIVE_TABLE_METHOD == 2)
  {
    constitutive_relations.constitutive_tables_ready = true;
    //Assign table limit to global definition.
    constitutive_relations.table_limit = INTERVALS_4_APPROX[NUM_OF_INTERVALS-1];
  }

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh, basemesh;
  MeshReaderH2D mloader;
  mloader.load(mesh_file.c_str(), &basemesh);
  
  // Perform initial mesh refinements.
  mesh.copy(&basemesh);
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_TOP, INIT_REF_NUM_BDY_TOP);

  // Initialize boundary conditions.
  RichardsEssentialBC bc_essential(BDY_TOP, H_ELEVATION, PULSE_END_TIME, H_INIT, STARTUP_TIME);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = Space<double>::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Solutions for the time stepping and the Newton's method.
  Solution<double> sln, ref_sln;
  InitialSolutionRichards sln_prev_time(&mesh, H_INIT);
  InitialSolutionRichards sln_prev_iter(&mesh, H_INIT);

  // Initialize the weak formulation.
  WeakForm<double>* wf;
  if (ITERATIVE_METHOD == 1) {
    if (TIME_INTEGRATION == 1) {
      info("Creating weak formulation for the Newton's method (implicit Euler in time).");
      wf = new WeakFormRichardsNewtonEuler(&constitutive_relations, time_step, &sln_prev_time, &mesh);
    }
    else {
      info("Creating weak formulation for the Newton's method (Crank-Nicolson in time).");
      wf = new WeakFormRichardsNewtonCrankNicolson(&constitutive_relations, time_step, &sln_prev_time, &mesh);
    }
  }
  else {
    if (TIME_INTEGRATION == 1) {
      info("Creating weak formulation for the Picard's method (implicit Euler in time).");
      wf = new WeakFormRichardsPicardEuler(&constitutive_relations, time_step, &sln_prev_iter, &sln_prev_time, &mesh);
    }
    else {
      info("Creating weak formulation for the Picard's method (Crank-Nicolson in time).");
      error("Not implemented yet.");
    }
  }

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err_est, graph_time_err_exact, 
    graph_time_dof, graph_time_cpu, graph_time_step;
 
  // Visualize the projection and mesh.
  ScalarView view("Initial condition", new WinGeom(0, 0, 630, 350));
  view.fix_scale_width(50);
  OrderView ordview("Initial mesh", new WinGeom(640, 0, 600, 350));
  view.show(&sln_prev_time, HERMES_EPS_VERYHIGH);
  ordview.show(&space);
  //MeshView mview("Mesh", new WinGeom(840, 0, 600, 350));
  //mview.show(&mesh);
  //View::wait();

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/time_step + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    info("---- Time step %d:", ts);

    // Time measurement.
    cpu_time.tick();

    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      switch (UNREF_METHOD) {
        case 1: mesh.copy(&basemesh);
                space.set_uniform_order(P_INIT);
                break;
        case 2: mesh.unrefine_all_elements();
                space.set_uniform_order(P_INIT);
                break;
        case 3: mesh.unrefine_all_elements();
                //space.adjust_element_order(-1, P_INIT);
                space.adjust_element_order(-1, -1, P_INIT, P_INIT);
                break;
        default: error("Wrong global derefinement method.");
      }

      ndof = Space<double>::get_num_dofs(&space);
    }

    // Spatial adaptivity loop. Note: sln_prev_time must not be touched during adaptivity.
    bool done = false;
    int as = 1;
    double err_est_rel;
    do
    {
      info("---- Time step %d, time step lenght %g, time %g (days), adaptivity step %d:", ts, time_step, current_time, as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Space<double>* ref_space = Space<double>::construct_refined_space(&space);
      ndof = Space<double>::get_num_dofs(ref_space);

      // Next we need to calculate the reference solution.
      // Newton's method:
      if(ITERATIVE_METHOD == 1) {
        double* coeff_vec = new double[ref_space->get_num_dofs()];
     
        // Calculate initial coefficient vector for Newton on the fine mesh.
        if (as == 1 && ts == 1) {
          info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
          OGProjection<double>::project_global(ref_space, &sln_prev_time, coeff_vec, matrix_solver);
        }
        else {
          info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
          OGProjection<double>::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver);
          delete ref_sln.get_space();
          delete ref_sln.get_mesh();
        }

        // Initialize the FE problem.
        DiscreteProblem<double> dp(wf, ref_space);

        // Perform Newton's iteration on the reference mesh. If necessary, 
        // reduce time step to make it converge, but then restore time step 
        // size to its original value.
        info("Performing Newton's iteration (tau = %g days):", time_step);
        bool success, verbose = true;
        double* save_coeff_vec = new double[ndof];
        
        // Save coefficient vector.
        for (int i=0; i < ndof; i++) 
          save_coeff_vec[i] = coeff_vec[i];

        bc_essential.set_current_time(current_time);
        
        // Perform Newton's iteration.
        info("Solving nonlinear problem:");
        Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver);
        bool newton_converged = false;
        while(!newton_converged)
        {
          try
          {
            newton.solve(coeff_vec, NEWTON_TOL, NEWTON_MAX_ITER);
            newton_converged = true;
          }
          catch(Hermes::Exceptions::Exception e)
          {
            e.printMsg();
            newton_converged = false;
          };
        
          if(!newton_converged)
          {
            // Restore solution from the beginning of time step.
            for (int i=0; i < ndof; i++) coeff_vec[i] = save_coeff_vec[i];
            // Reducing time step to 50%.
            info("Reducing time step size from %g to %g days for the rest of this time step.", 
                 time_step, time_step * time_step_dec);
            time_step *= time_step_dec;
            // If time_step less than the prescribed minimum, stop.
            if (time_step < time_step_min) error("Time step dropped below prescribed minimum value.");
          }
        }
        // Delete the saved coefficient vector.
        delete [] save_coeff_vec;

        // Translate the resulting coefficient vector 
        // into the desired reference solution. 
        Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);

        // Cleanup.
        delete [] coeff_vec;
      }
      else {
        // Calculate initial condition for Picard on the fine mesh.
        if (as == 1 && ts == 1) {
          info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
          OGProjection<double>::project_global(ref_space, &sln_prev_time, &sln_prev_iter, matrix_solver);
        }
        else {
          info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
          OGProjection<double>::project_global(ref_space, &ref_sln, &sln_prev_iter, matrix_solver);
        }

        // Perform Picard iteration on the reference mesh. If necessary, 
        // reduce time step to make it converge, but then restore time step 
        // size to its original value.
        info("Performing Picard's iteration (tau = %g days):", time_step);
        bool verbose = true;

        bc_essential.set_current_time(current_time);

        DiscreteProblem<double> dp(wf, ref_space);
        PicardSolver<double> picard(&dp, &sln_prev_iter, matrix_solver);
        picard.set_verbose_output(verbose);
        while(!picard.solve(PICARD_TOL, PICARD_MAX_ITER)) 
          {
          // Restore solution from the beginning of time step.
          sln_prev_iter.copy(&sln_prev_time);
          // Reducing time step to 50%.
          info("Reducing time step size from %g to %g days for the rest of this time step", time_step, time_step * time_step_inc);
          time_step *= time_step_dec;
          // If time_step less than the prescribed minimum, stop.
          if (time_step < time_step_min) error("Time step dropped below prescribed minimum value.");
        }	

        ref_sln.copy(&sln_prev_iter);
      }


      /*** ADAPTIVITY ***/

      // Project the fine mesh solution on the coarse mesh.
      info("Projecting fine mesh solution on coarse mesh for error calculation.");
      if(space.get_mesh() == NULL) error("it is NULL");
      OGProjection<double>::project_global(&space, &ref_sln, &sln, matrix_solver);

      // Calculate element errors.
      info("Calculating error estimate."); 
      Adapt<double>* adaptivity = new Adapt<double>(&space);
      
      // Calculate error estimate wrt. fine mesh solution.
      err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, space_err_est_rel: %g%%", 
        Space<double>::get_num_dofs(&space), Space<double>::get_num_dofs(ref_space), err_est_rel);

      // If space_err_est too large, adapt the mesh.
      if (err_est_rel < ERR_STOP) done = true;
      else {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (Space<double>::get_num_dofs(&space) >= NDOF_STOP) {
          done = true;
          break;
        }
        as++;
      }

      delete adaptivity;
    }
    while (!done);

    // Add entries to graphs.
    graph_time_err_est.add_values(current_time, err_est_rel);
    graph_time_err_est.save("time_error_est.dat");
    graph_time_dof.add_values(current_time, Space<double>::get_num_dofs(&space));
    graph_time_dof.save("time_dof.dat");
    graph_time_cpu.add_values(current_time, cpu_time.accumulated());
    graph_time_cpu.save("time_cpu.dat");
    graph_time_step.add_values(current_time, time_step);
    graph_time_step.save("time_step_history.dat");

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time %g days", current_time);
    view.set_title(title);
    view.show(&sln, HERMES_EPS_VERYHIGH);
    sprintf(title, "Mesh, time %g days", current_time);
    ordview.set_title(title);
    ordview.show(&space);
    
    // Save complete Solution.
    char* filename = new char[100];
    sprintf(filename, "outputs/tsln_%f.dat", current_time);
    sln.save(filename);
    info("Solution at time %g saved to file %s.", current_time, filename);

    // Copy new reference level solution into sln_prev_time.
    // This starts new time step.
    sln_prev_time.copy(&ref_sln);

    // Updating time step. Note that time_step might have been reduced during adaptivity.
    current_time += time_step;

    // Increase time step.
    if (time_step*time_step_inc < time_step_max) {
      info("Increasing time step from %g to %g days.", time_step, time_step * time_step_inc);
      time_step *= time_step_inc;
    }
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
