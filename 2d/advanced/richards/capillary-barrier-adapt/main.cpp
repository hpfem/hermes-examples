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
//
//  BC: Dirichlet, given by the initial condition.
//
//  The following parameters can be changed:

// Choose either CONSTITUTIVE_GARDNER or CONSTITUTIVE_GENUCHTEN.
CONSTITUTIVE_RELATIONS constitutive_relations = CONSTITUTIVE_GENUCHTEN;

// Choose full domain or half domain.
// const char* mesh_file = "domain-full.mesh";
std::string mesh_file = "domain-half.mesh";

// Methods.
// 1 = Newton, 2 = Picard.
const int ITERATIVE_METHOD = 2;		          
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
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  

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
// Global time variable initialized with first time step.
double current_time = time_step;                  

// Boundary markers.
const std::string BDY_TOP = "1";
const std::string BDY_RIGHT = "2";
const std::string BDY_BOTTOM = "3";
const std::string BDY_LEFT = "4";

// Main function.
int main(int argc, char* argv[])
{
  ConstitutiveRelations* relations;
  if(constitutive_relations == CONSTITUTIVE_GENUCHTEN)
    relations = new ConstitutiveGenuchten(LOW_LIMIT, POLYNOMIALS_READY, CONSTITUTIVE_TABLE_METHOD, 
                NUM_OF_INSIDE_PTS, TABLE_LIMIT, ALPHA_vals, N_vals, M_vals, K_S_vals, THETA_R_vals,
    THETA_S_vals, STORATIVITY_vals, TABLE_PRECISION, MATERIAL_COUNT, POLYNOMIALS_ALLOCATED, 
                NUM_OF_INTERVALS, ITERATIVE_METHOD);
  else
    relations = new ConstitutiveGardner(K_S_vals[0], ALPHA_vals[0], THETA_S_vals[0], THETA_R_vals[0], 
                CONSTITUTIVE_TABLE_METHOD);

  // Points to be used for polynomial approximation of K(h).
  double* points = new double[NUM_OF_INSIDE_PTS];

  // The van Genuchten + Mualem K(h) function is approximated by polynomials close 
  // to zero in case of CONSTITUTIVE_TABLE_METHOD==1.
  // In case of CONSTITUTIVE_TABLE_METHOD==2, all constitutive functions are approximated by polynomials.
  info("Initializing polynomial approximations.");
  for (int i=0; i < MATERIAL_COUNT; i++)
    relations->init_polynomials(6 + NUM_OF_INSIDE_PTS, LOW_LIMIT, points, NUM_OF_INSIDE_PTS, i);
  relations->polynomials_ready = true;
  if (CONSTITUTIVE_TABLE_METHOD == 2) {
    relations->constitutive_tables_ready = true ;
    //Assign table limit to global definition.
    relations->table_limit = INTERVALS_4_APPROX[NUM_OF_INTERVALS-1];
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
      wf = new WeakFormRichardsNewtonEuler(relations, time_step, &sln_prev_time, &mesh);
    }
    else {
      info("Creating weak formulation for the Newton's method (Crank-Nicolson in time).");
      wf = new WeakFormRichardsNewtonCrankNicolson(relations, time_step, &sln_prev_time, &mesh);
    }
  }
  else {
    if (TIME_INTEGRATION == 1) {
      info("Creating weak formulation for the Picard's method (implicit Euler in time).");
      wf = new WeakFormRichardsPicardEuler(relations, time_step, &sln_prev_iter, &sln_prev_time, &mesh);
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
        double* coeff_vec = new double[Space<double>::get_num_dofs(ref_space)];
     
        // Calculate initial coefficient vector for Newton on the fine mesh.
        if (as == 1 && ts == 1) {
          info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
          OGProjection<double>::project_global(ref_space, &sln_prev_time, coeff_vec, matrix_solver_type);
        }
        else {
          info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
          OGProjection<double>::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver_type);
          delete ref_sln.get_mesh();
        }

        // Initialize the FE problem.
        DiscreteProblem<double> dp(wf, ref_space);

        // Set up the solver, matrix, and rhs according to the solver selection.
        SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver_type);
        Vector<double>* rhs = create_vector<double>(matrix_solver_type);
        LinearSolver<double>* solver = create_linear_solver<double>(matrix_solver_type, matrix, rhs);

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
        Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver_type);
        bool Newton_converged = false;
        while(!Newton_converged)
        {
          try
          {
            newton.solve(coeff_vec, NEWTON_TOL, NEWTON_MAX_ITER);
            Newton_converged = true;
          }
          catch(Hermes::Exceptions::Exception e)
          {
            e.printMsg();
            Newton_converged = false;
          };
        
          if(!Newton_converged)
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
        Solution<double>::vector_to_solution(coeff_vec, ref_space, &ref_sln);

        // Cleanup.
        delete [] coeff_vec;
        delete solver;
        delete matrix;
        delete rhs;
      }
      else {
        // Calculate initial condition for Picard on the fine mesh.
        if (as == 1 && ts == 1) {
          info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
          OGProjection<double>::project_global(ref_space, &sln_prev_time, &sln_prev_iter, matrix_solver_type);
        }
        else {
          info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
          OGProjection<double>::project_global(ref_space, &ref_sln, &sln_prev_iter, matrix_solver_type);
        }

        // Perform Picard iteration on the reference mesh. If necessary, 
        // reduce time step to make it converge, but then restore time step 
        // size to its original value.
        info("Performing Picard's iteration (tau = %g days):", time_step);
        bool verbose = true;

        bc_essential.set_current_time(current_time);

        DiscreteProblem<double> dp(wf, ref_space);
        PicardSolver<double> picard(&dp, &sln_prev_iter, matrix_solver_type);
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
      OGProjection<double>::project_global(&space, &ref_sln, &sln, matrix_solver_type);

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
      delete ref_space;
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
