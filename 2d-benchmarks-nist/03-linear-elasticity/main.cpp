#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;

//  This is the third in the series of NIST benchmarks with known exact solutions.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: Linear elasticity coupled system of two equations given below
//
//  -E \frac{1-nu^2}{1-2*nu} \frac{\partial^{2} u}{\partial x^{2}} - E\frac{1-nu^2}{2-2*nu} \frac{\partial^{2} u}{\partial y^{2}} 
//  -E \frac{1-nu^2}{(1-2*nu)(2-2*nu)} \frac{\partial^{2} v}{\partial x \partial y} - F_{x} = 0
//
//  -E \frac{1-nu^2}{2-2*nu} \frac{\partial^{2} v}{\partial x^{2}} - E\frac{1-nu^2}{1-2*nu} \frac{\partial^{2} v}{\partial y^{2}} 
//  -E \frac{1-nu^2}{(1-2*nu)(2-2*nu)} \frac{\partial^{2} u}{\partial x \partial y} - F_{y} = 0
//
//  where F_{x} = F_{y} = 0.
//
//  Known exact solution for mode 1: 
//  u(x, y) = \frac{1}{2G} r^{\lambda}[(k - Q(\lambda + 1))cos(\lambda \theta) - \lambda cos((\lambda - 2) \theta)]
//  v(x, y) = \frac{1}{2G} r^{\lambda}[(k + Q(\lambda + 1))sin(\lambda \theta) + \lambda sin((\lambda - 2) \theta)]
//  here \lambda = 0.5444837367825, Q = 0.5430755788367.
//
//  Known exact solution for mode 2: 
//  u(x, y) =  \frac{1}{2G} r^{\lambda}[(k - Q(\lambda + 1))sin(\lambda \theta) - \lambda sin((\lambda - 2) \theta)]
//  v(x, y) = -\frac{1}{2G} r^{\lambda}[(k + Q(\lambda + 1))cos(\lambda \theta) + \lambda cos((\lambda - 2) \theta)]
//  here \lambda = 0.9085291898461, Q = -0.2189232362488.
//
//  Domain: square (-1, 1)^2, with a slit from (0, 0) to (1, 0).
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

// If this is defined, mode-1 solution is calculated, otherwise 
// mode-2 solution. See Mitchell's paper for details.
#define MODE_1                                    

// Initial polynomial degree for u.
const int P_INIT_U = 2;                           
// Initial polynomial degree for v.
const int P_INIT_V = 2;                                                     
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
// Stopping criterion for adaptivity (rel. error tolerance between the
// reference mesh and coarse mesh solution in percent).
const double ERR_STOP = 1.0;                      
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
const double NEWTON_TOLERANCE = 1E-3;           
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  

// Problem parameters.

// Young modulus.
const double E = 1.0;
// Poisson ratio.
const double nu = 0.3;
#ifdef MODE_1
// lambda for mode-1 solution.
const double lambda = 0.5444837367825;
// mu for mode-1 solution (mu is the same as G)
const double mu = E / (2 * (1 + nu));
// Q for mode-1 solution.
const double Q = 0.5430755788367;
#else
// lambda for mode-2 solution.
const double lambda = 0.9085291898461;
// mu for mode-2 solution (mu is the same as G).
const double mu = E / (2 * (1 + nu));
// Q for mode-2 solution.
const double Q = -0.2189232362488;              
#endif

int main(int argc, char* argv[])
{
  // Load the mesh->
  MeshSharedPtr u_mesh(new Mesh), v_mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("elasticity.mesh", u_mesh);

  // Create initial mesh (master mesh).
  v_mesh->copy(u_mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) {
    u_mesh->refine_all_elements();
    v_mesh->refine_all_elements();
  }

  // Set exact solution for each displacement component.
  MeshFunctionSharedPtr<double> exact_u(new CustomExactSolutionU(u_mesh, E, nu, lambda, Q));
  MeshFunctionSharedPtr<double> exact_v(new CustomExactSolutionV(v_mesh, E, nu, lambda, Q));

  // Initialize the weak formulation.
  // NOTE: We pass all four parameters (temporarily) 
  // since in Mitchell's paper (NIST benchmarks) they 
  // are mutually inconsistent.
  CustomWeakFormElasticityNIST wf(E, nu, mu, lambda);

  // Initialize boundary conditions.
  DefaultEssentialBCNonConst<double> bc_u("Bdy", exact_u);
  EssentialBCs<double> bcs_u(&bc_u);
  DefaultEssentialBCNonConst<double> bc_v("Bdy", exact_v);
  EssentialBCs<double> bcs_v(&bc_v);

  // Create H1 spaces with default shapeset for both displacement components.
  SpaceSharedPtr<double> u_space(new H1Space<double>(u_mesh, &bcs_u, P_INIT_U));
  SpaceSharedPtr<double> v_space(new H1Space<double>(v_mesh, &bcs_v, P_INIT_V));

  // Initialize approximate solution.
  MeshFunctionSharedPtr<double> u_sln(new Solution<double>());
  MeshFunctionSharedPtr<double> u_ref_sln(new Solution<double>());
  MeshFunctionSharedPtr<double> v_sln(new Solution<double>());
  MeshFunctionSharedPtr<double> v_ref_sln(new Solution<double>());

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  Views::ScalarView s_view_u("Solution for u", new WinGeom(0, 0, 440, 350));
  s_view_u.show_mesh(false);
  Views::OrderView  o_view_u("Mesh for u", new WinGeom(450, 0, 420, 350));
  Views::ScalarView s_view_v("Solution for v", new WinGeom(0, 405, 440, 350));
  s_view_v.show_mesh(false);
  Views::OrderView  o_view_v("Mesh for v", new WinGeom(450, 405, 420, 350));
  Views::ScalarView mises_view("Von Mises stress [Pa]", new WinGeom(880, 0, 440, 350));
  mises_view.fix_scale_width(50);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    cpu_time.tick();

    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator refMeshCreatorU(u_mesh);
    MeshSharedPtr ref_u_mesh = refMeshCreatorU.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreatorU(u_space, ref_u_mesh);
    SpaceSharedPtr<double> ref_u_space = refSpaceCreatorU.create_ref_space();

    Mesh::ReferenceMeshCreator refMeshCreatorV(u_mesh);
    MeshSharedPtr ref_v_mesh = refMeshCreatorV.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreatorV(v_space, ref_v_mesh);
    SpaceSharedPtr<double> ref_v_space = refSpaceCreatorV.create_ref_space();

    Hermes::vector<SpaceSharedPtr<double> > ref_spaces(ref_u_space, ref_v_space);

    int ndof_ref = Space<double>::get_num_dofs(ref_spaces);

    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d (%d DOF):", as, ndof_ref);
    cpu_time.tick();

    Hermes::Mixins::Loggable::Static::info("Solving on reference mesh->");

    // Assemble the discrete problem.    
    DiscreteProblem<double> dp(&wf, ref_spaces);

    NewtonSolver<double> newton(&dp);
    //newton.set_verbose_output(false);
    newton.set_tolerance(NEWTON_TOLERANCE);

    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    };

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), ref_spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(u_ref_sln, v_ref_sln));

    cpu_time.tick();
    Hermes::Mixins::Loggable::Static::info("Solution: %g s", cpu_time.last());

    // Project the fine mesh solution onto the coarse mesh->
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate and exact error.");
    OGProjection<double> ogProjection; ogProjection.project_global(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space), 
      Hermes::vector<MeshFunctionSharedPtr<double> >(u_ref_sln, v_ref_sln), 
      Hermes::vector<MeshFunctionSharedPtr<double> >(u_sln, v_sln));

    // Calculate element errors and total error estimate.
    Hermes::vector<double> err_est_rel;
    Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space));
    double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<MeshFunctionSharedPtr<double> >(u_sln, v_sln), 
      Hermes::vector<MeshFunctionSharedPtr<double> >(u_ref_sln, v_ref_sln), &err_est_rel) * 100.;

    // Calculate exact error for each solution component and the total exact error.
    Hermes::vector<double> err_exact_rel;
    bool solutions_for_adapt = false;
    double err_exact_rel_total = adaptivity->calc_err_exact(Hermes::vector<MeshFunctionSharedPtr<double> >(u_sln, v_sln), 
      Hermes::vector<MeshFunctionSharedPtr<double> >(exact_u, exact_v), 
      &err_exact_rel, solutions_for_adapt) * 100.;

    cpu_time.tick();
    Hermes::Mixins::Loggable::Static::info("Error calculation: %g s", cpu_time.last());

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse[u]: %d, ndof_fine[u]: %d",
      u_space->Space<double>::get_num_dofs(), Space<double>::get_num_dofs(ref_u_space));
    Hermes::Mixins::Loggable::Static::info("err_est_rel[u]: %g%%, err_exact_rel[u]: %g%%", err_est_rel[0]*100, err_exact_rel[0]*100);
    Hermes::Mixins::Loggable::Static::info("ndof_coarse[v]: %d, ndof_fine[v]: %d",
      v_space->Space<double>::get_num_dofs(), Space<double>::get_num_dofs(ref_v_space));
    Hermes::Mixins::Loggable::Static::info("err_est_rel[v]: %g%%, err_exact_rel[v]: %g%%", err_est_rel[1]*100, err_exact_rel[1]*100);
    Hermes::Mixins::Loggable::Static::info("ndof_coarse_total: %d, ndof_fine_total: %d",
      Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space)), Space<double>::get_num_dofs(ref_spaces));
    Hermes::Mixins::Loggable::Static::info("err_est_rel_total: %g%%, err_est_exact_total: %g%%", err_est_rel_total, err_exact_rel_total);

    // Time measurement.
    cpu_time.tick();
    double accum_time = cpu_time.accumulated();

    // View the coarse mesh solution and polynomial orders.
    s_view_u.show(u_sln); 
    o_view_u.show(u_space);
    s_view_v.show(v_sln); 
    o_view_v.show(v_space);
    MeshFunctionSharedPtr<double> stress(new VonMisesFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(u_sln, v_sln), lambda, mu));
    mises_view.show(stress, HERMES_EPS_NORMAL, H2D_FN_VAL_0, u_sln, v_sln, 0.03);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space)), err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");
    graph_dof_exact.add_values(Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space)), err_exact_rel_total);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel_total);
    graph_cpu_exact.save("conv_cpu_exact.dat");

    cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

    // If err_est too large, adapt the mesh->
    if (err_est_rel_total < ERR_STOP) 
      done = true;
    else 
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh->");
      done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector), 
        THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space)) >= NDOF_STOP) done = true;

    cpu_time.tick();
    Hermes::Mixins::Loggable::Static::info("Adaptation: %g s", cpu_time.last());

    // Increase the counter of adaptivity steps.
    if (done == false)  
      as++;

    delete adaptivity;
  }
  while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}
