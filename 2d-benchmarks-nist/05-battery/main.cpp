#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;

//  This is the fifth in the series of NIST benchmarks with unknown exact solution.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: -\frac{\partial }{\partial x}\left(p(x, y)\frac{\partial u}{\partial x}\right)
//       -\frac{\partial }{\partial y}\left(q(x, y)\frac{\partial u}{\partial y}\right) - f = 0.
//
//  Exact solution: unknown.
//
//  Domain: square (0, 8.4) x (0, 24), see the file "battery.mesh".
//
//  BC: Zero Neumann on left edge, Newton on the rest of the boundary:
//
//  The following parameters can be changed:

// Set to "false" to suppress Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;           
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;             
// Initial polynomial degree of mesh elements.
const int P_INIT = 2;
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
// More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int STRATEGY = 0;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
// H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;        
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// Stopping criterion for adaptivity.
const double ERR_STOP = 1.0;                      
// This parameter influences the selection of candidates in hp-adaptivity. 
// Default value is 1.0.
const double CONV_EXP = 0.3;                      
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK; 

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("battery.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, P_INIT);

  // Initialize weak formulation.
  CustomWeakFormPoisson wf("e1", "e2", "e3", "e4", "e5", 
                           "Bdy_left", "Bdy_top", "Bdy_right", "Bdy_bottom", &mesh);

  // Initialize coarse and fine mesh solution.
  Solution<double> sln, ref_sln;
  
  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 320, 600));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  OrderView  oview("Polynomial orders", new WinGeom(330, 0, 300, 600));
  
  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;
  
  // Time measurement.
  TimePeriod cpu_time;

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    
    // Time measurement.
    cpu_time.tick();

    // Construct globally refined mesh and setup fine mesh space.
    Space<double>* ref_space = Space<double>::construct_refined_space(&space);
    int ndof_ref = ref_space->get_num_dofs();

    // Initialize fine mesh problem.
    info("Solving on fine mesh.");
    DiscreteProblem<double> dp(&wf, ref_space);
    
    NewtonSolver<double> newton(&dp, matrix_solver_type);
    newton.set_verbose_output(false);

    // Initial coefficient vector for the Newton's method.  
    double* coeff_vec = new double[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(double));

    // Perform Newton's iteration.
    try
    {
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    }

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);
    
    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting fine mesh solution on coarse mesh.");
    OGProjection<double>::project_global(&space, &ref_sln, &sln, matrix_solver_type);

    // Time measurement.
    cpu_time.tick();

    // VTK output.
    if (VTK_VISUALIZATION) 
    {
      // Output solution in VTK format.
      Views::Linearizer lin;
      char* title = new char[100];
      sprintf(title, "sln-%d.vtk", as);
      lin.save_solution_vtk(&sln, title, "Potential", false);
      info("Solution in VTK format saved to file %s.", title);

      // Output mesh and element orders in VTK format.
      Views::Orderizer ord;
      sprintf(title, "ord-%d.vtk", as);
      ord.save_orders_vtk(&space, title);
      info("Element orders in VTK format saved to file %s.", title);
    }

    // View the coarse mesh solution and polynomial orders.
    if (HERMES_VISUALIZATION) 
    {
      sview.show(&sln);
      oview.show(&space);
    }

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate.");
    Adapt<double> adaptivity(&space);
    bool solutions_for_adapt = true;
    // In the following function, the Boolean parameter "solutions_for_adapt" determines whether
    // the calculated errors are intended for use with adaptivity (this may not be the case, for example,
    // when error wrt. an exact solution is calculated). The default value is solutions_for_adapt = true,
    // The last parameter "error_flags" determine whether the total and element errors are treated as
    // absolute or relative. Its default value is error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL.
    // In subsequent examples and benchmarks, these two parameters will be often used with
    // their default values, and thus they will not be present in the code explicitly.
    double err_est_rel = adaptivity.calc_err_est(&sln, &ref_sln, solutions_for_adapt,
                         HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      space.get_num_dofs(), ref_space->get_num_dofs(), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    cpu_time.tick();    
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    graph_dof.add_values(space.get_num_dofs(), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    
    // Skip the time spent to save the convergence graphs.
    cpu_time.tick(HERMES_SKIP);

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) 
      done = true;
    else
    {
      info("Adapting coarse mesh.");
      done = adaptivity.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of performed adaptivity steps.
      if (done == false)  
        as++;
    }
    if (space.get_num_dofs() >= NDOF_STOP) 
      done = true;

    // Clean up.
    delete [] coeff_vec;
    // Keep the mesh from final step to allow further work with the final fine mesh solution.
    if(done == false) 
      delete ref_space->get_mesh(); 
    delete ref_space;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the fine mesh solution - final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(&ref_sln);

  // Wait for all views to be closed.
  Views::View::wait();

  return 0;
}

