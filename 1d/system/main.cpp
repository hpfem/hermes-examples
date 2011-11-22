#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace RefinementSelectors;

// This example is an analogy to the 2D tutorial example P04-adapt/03-system.
// It explains how to use the multimesh adaptive hp-FEM for 1D problems whose
// solution components behave very differently, so that using the same mesh 
// for both is a waste of DOFs.
//
// PDE: Linearized FitzHugh-Nagumo equation
//      -d_u^2 u'' - f(u) + \sigma v - g1 = 0,
//      -d_v^2 v'' - u + v - g2 = 0.
// In the original equation, f(u) = \lambda u - u^3 - \kappa. For
// simplicity, here we just take f(u) = u.
//
// Domain: Interval (-1, 1).
//
// BC: Both solution components are zero on the boundary.
//
// Exact solution: The functions g1 and g2 were calculated so that
//                 the exact solution is:
//        U(x) = cos(M_PI*x/2)
//        V(x) = 1 - (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K))
// Note: V(x) is the exact solution of the 1D singularly perturbed equation
//       -u'' + K*K*u = K*K in (-1, 1) with zero Dirichlet BC.
//
// The following parameters can be changed: In particular, compare hp- and
// h-adaptivity via the CAND_LIST option, and compare the multi-mesh vs.
// single-mesh using the MULTI parameter.

// Should exact solution be used and exact error calculated?
#define WITH_EXACT_SOLUTION                       
// Initial polynomial degree for u.
const int P_INIT_U = 1;                           
// Initial polynomial degree for v.
const int P_INIT_V = 1;                           
// Number of initial boundary refinements
const int INIT_REF_BDY = 5;                       
// MULTI = true  ... use multi-mesh,
// MULTI = false ... use single-mesh.
// Note: In the single mesh option, the meshes are
// forced to be geometrically the same but the
// polynomial degrees can still vary.
const bool MULTI = true;                          
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
const int STRATEGY = 1;                           
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
const double CONV_EXP = 1;                        
// Stopping criterion for adaptivity.
const double ERR_STOP = 0.1;                      
// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 1000;                       
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  

// Problem parameters.
const double D_u = 1;
const double D_v = 1;
const double SIGMA = 1;
const double LAMBDA = 1;
const double KAPPA = 1;
const double K = 100.;

int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh u_mesh, v_mesh;
  MeshReaderH1DXML mloader;
  mloader.load("domain.xml", &u_mesh);
  if (MULTI == false) 
  {
    u_mesh.refine_towards_boundary("Left", INIT_REF_BDY);
    // Minus one for the sake of mesh symmetry.     
    if (INIT_REF_BDY > 1) u_mesh.refine_towards_boundary("Right", INIT_REF_BDY - 1);  
  }

  // Create initial mesh (master mesh).
  v_mesh.copy(&u_mesh);

  // Initial mesh refinements in the v_mesh towards the boundary.
  if (MULTI == true) 
  {
    v_mesh.refine_towards_boundary("Left", INIT_REF_BDY);
    // Minus one for the sake of mesh symmetry. 
    if (INIT_REF_BDY > 1) v_mesh.refine_towards_boundary("Right", INIT_REF_BDY - 1);  
  }

#ifdef WITH_EXACT_SOLUTION
  // Set exact solutions.
  ExactSolutionFitzHughNagumo1 exact_u(&u_mesh);
  ExactSolutionFitzHughNagumo2 exact_v(&v_mesh, K);
#endif

  // Define right-hand sides.
  CustomRightHandSide1 g1(K, D_u, SIGMA);
  CustomRightHandSide2 g2(K, D_v);

  // Initialize the weak formulation.
  CustomWeakForm wf(&g1, &g2);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_u(Hermes::vector<std::string>("Left", "Right"), 0.0);
  EssentialBCs<double> bcs_u(&bc_u);
  DefaultEssentialBCConst<double> bc_v(Hermes::vector<std::string>("Left", "Right"), 0.0);
  EssentialBCs<double> bcs_v(&bc_v);

  // Create H1 spaces with default shapeset for both displacement components.
  H1Space<double> u_space(&u_mesh, &bcs_u, P_INIT_U);
  H1Space<double> v_space(MULTI ? &v_mesh : &u_mesh, &bcs_v, P_INIT_V);

  // Initialize coarse and reference mesh solutions.
  Solution<double> u_sln, v_sln, u_ref_sln, v_ref_sln;

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  Views::ScalarView s_view_0("Solution[0]", new Views::WinGeom(0, 0, 440, 350));
  s_view_0.show_mesh(false);
  Views::OrderView  o_view_0("Mesh[0]", new Views::WinGeom(450, 0, 420, 350));
  Views::ScalarView s_view_1("Solution[1]", new Views::WinGeom(880, 0, 440, 350));
  s_view_1.show_mesh(false);
  Views::OrderView o_view_1("Mesh[1]", new Views::WinGeom(1330, 0, 420, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est; 
#ifdef WITH_EXACT_SOLUTION
  SimpleGraph graph_dof_exact, graph_cpu_exact;
#endif

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    // FIXME: This should be increase in the x-direction only.
    int order_increase = 1;          
    // FIXME: This should be '2' but that leads to a segfault.
    int refinement_type = 0;         
    Hermes::vector<Space<double> *>* ref_spaces = 
        Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&u_space, &v_space),
        order_increase, refinement_type);
    Space<double>* u_ref_space = (*ref_spaces)[0];
    Space<double>* v_ref_space = (*ref_spaces)[1];
    Hermes::vector<const Space<double> *> ref_spaces_const((*ref_spaces)[0], (*ref_spaces)[1]);
    int ndof_ref = Space<double>::get_num_dofs(ref_spaces_const);

    // Initialize reference problem.
    info("Solving on reference mesh.");
    DiscreteProblem<double> dp(&wf, ref_spaces_const);

    NewtonSolver<double> newton(&dp, matrix_solver);
    newton.set_verbose_output(true);

    // Time measurement.
    cpu_time.tick();
    
    // Perform Newton's iteration.
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    };

    // Translate the resulting coefficient vector into the Solution<double> sln.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), ref_spaces_const, 
        Hermes::vector<Solution<double> *>(&u_ref_sln, &v_ref_sln));

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solutions on coarse meshes.");
    OGProjection<double>::project_global(Hermes::vector<const Space<double> *>(&u_space, &v_space), 
        Hermes::vector<Solution<double> *>(&u_ref_sln, &v_ref_sln), 
        Hermes::vector<Solution<double> *>(&u_sln, &v_sln), matrix_solver); 
   
    cpu_time.tick();

    // View the coarse mesh solution and polynomial orders.
    s_view_0.show(&u_sln); 
    o_view_0.show(&u_space);
    s_view_1.show(&v_sln); 
    o_view_1.show(&v_space);

    // Calculate element errors.
    info("Calculating error estimate and exact error."); 
    Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<Space<double> *>(&u_space, &v_space));
    
    // Calculate error estimate for each solution component and the total error estimate.
    Hermes::vector<double> err_est_rel;
    double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution<double> *>(&u_sln, &v_sln), 
        Hermes::vector<Solution<double> *>(&u_ref_sln, &v_ref_sln), &err_est_rel) * 100;

#ifdef WITH_EXACT_SOLUTION
    // Calculate exact error for each solution component and the total exact error.
    Hermes::vector<double> err_exact_rel;
    bool solutions_for_adapt = false;
    double err_exact_rel_total = adaptivity->calc_err_exact(Hermes::vector<Solution<double> *>(&u_sln, &v_sln), 
        Hermes::vector<Solution<double> *>(&exact_u, &exact_v), &err_exact_rel, solutions_for_adapt) * 100;
#endif

    // Time measurement.
    cpu_time.tick();

    // Report results.
#ifdef WITH_EXACT_SOLUTION
    info("ndof_coarse[0]: %d, ndof_fine[0]: %d",
         u_space.get_num_dofs(), u_ref_space->get_num_dofs());
    info("err_est_rel[0]: %g%%, err_exact_rel[0]: %g%%", err_est_rel[0]*100, err_exact_rel[0]*100);
    info("ndof_coarse[1]: %d, ndof_fine[1]: %d",
         v_space.get_num_dofs(), v_ref_space->get_num_dofs());
    info("err_est_rel[1]: %g%%, err_exact_rel[1]: %g%%", err_est_rel[1]*100, err_exact_rel[1]*100);
    info("ndof_coarse_total: %d, ndof_fine_total: %d",
         Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u_space, &v_space)), 
         Space<double>::get_num_dofs(ref_spaces_const));
    info("err_est_rel_total: %g%%, err_est_exact_total: %g%%", err_est_rel_total, err_exact_rel_total);
#else
    info("ndof_coarse[0]: %d, ndof_fine[0]: %d", u_space.get_num_dofs(), u_ref_space->get_num_dofs());
    info("err_est_rel[0]: %g%%", err_est_rel[0]*100);
    info("ndof_coarse[1]: %d, ndof_fine[1]: %d", v_space.get_num_dofs(), v_ref_space->get_num_dofs());
    info("err_est_rel[1]: %g%%", err_est_rel[1]*100);
    info("ndof_coarse_total: %d, ndof_fine_total: %d",
         Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u_space, &v_space)), 
         Space<double>::get_num_dofs(*ref_spaces));
    info("err_est_rel_total: %g%%", err_est_rel_total);
#endif

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u_space, &v_space)), 
        err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");
#ifdef WITH_EXACT_SOLUTION
    graph_dof_exact.add_values(Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u_space, &v_space)), 
        err_exact_rel_total);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel_total);
    graph_cpu_exact.save("conv_cpu_exact.dat");
#endif

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP) 
      done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector), 
          THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u_space, &v_space)) 
        >= NDOF_STOP) done = true;

    // Clean up.
    delete adaptivity;
    for(unsigned int i = 0; i < ref_spaces->size(); i++)
      delete (*ref_spaces)[i]->get_mesh();
    delete ref_spaces;
    
    // Increase counter.
    as++;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
