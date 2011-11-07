#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// This example uses adaptive multimesh hp-FEM to solve a simple problem
// of linear elasticity. Note that since both displacement components
// have similar qualitative behavior, the advantage of the multimesh 
// discretization is less striking than for example in the tutorial 
// example P04-linear-adapt/02-system-adapt.
//
// PDE: Lame equations of linear elasticity. No external forces, the 
//      object is loaded with its own weight.
//
// BC: u_1 = u_2 = 0 on Gamma_1 (left edge)
//     du_1/dn = du_2/dn = 0 elsewhere, including two horizontal
//               cracks inside the domain. The width of the cracks
//               is currently zero, it can be set in the mesh file
//               via the parameter 'w'.
//
// The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Initial polynomial degree of mesh elements (u-displacement).
const int P_INIT_U1 = 2;                          
// Initial polynomial degree of mesh elements (v-displacement).
const int P_INIT_U2 = 2;                          
// true = use multi-mesh, false = use single-mesh.
// Note: in the single mesh option, the meshes are
// forced to be geometrically the same but the
// polynomial degrees can still vary.
const bool MULTI = true;                          
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD_MULTI = 0.35;              
const double THRESHOLD_SINGLE = 0.7;              
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
const double ERR_STOP = 0.1;                      
// Adaptivity process stops when the number of degrees of freedom grows.
const int NDOF_STOP = 60000;                      
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-6;                            
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                           
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Problem parameters.
// Young modulus for steel: 200 GPa.
const double E  = 200e9;                          
// Poisson ratio.
const double nu = 0.3;                            
// Gravitational acceleration.
const double g1 = -9.81;                          
// Material density in kg / m^3. 
const double rho = 8000;                          
// Top surface force in x-direction.
const double f0  = 0;                             
// Top surface force in y-direction.
const double f1  = 0;

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh u1_mesh, u2_mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &u1_mesh);

  // Perform initial uniform mesh refinement.
  for (int i=0; i < INIT_REF_NUM; i++) u1_mesh.refine_all_elements();

  // Create initial mesh for the vertical displacement component.
  // This also initializes the multimesh hp-FEM.
  u2_mesh.copy(&u1_mesh);

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> zero_disp("bdy left", 0.0);
  EssentialBCs<double> bcs(&zero_disp);

  // Create x- and y- displacement space using the default H1 shapeset.
  H1Space<double> u1_space(&u1_mesh, &bcs, P_INIT_U1);
  H1Space<double> u2_space(&u2_mesh, &bcs, P_INIT_U2);
  info("ndof = %d.", Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u1_space, &u2_space)));

  // Initialize the weak formulation.
  // NOTE; These weak forms are identical to those in example P01-linear/08-system.
  CustomWeakFormLinearElasticity wf(E, nu, rho*g1, "bdy_top", f0, f1);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, Hermes::vector<const Space<double> *>(&u1_space, &u2_space));

  // Initialize coarse and reference mesh solutions.
  Solution<double> u1_sln, u2_sln;
  Solution<double> u1_sln_ref, u2_sln_ref;

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView s_view_0("Solution (x-displacement)", new WinGeom(0, 0, 700, 350));
  s_view_0.show_mesh(false);
  ScalarView s_view_1("Solution (y-displacement)", new WinGeom(760, 0, 700, 350));
  s_view_1.show_mesh(false);
  OrderView  o_view_0("Mesh (x-displacement)", new WinGeom(410, 0, 700, 350));
  OrderView  o_view_1("Mesh (y-displacement)", new WinGeom(1170, 0, 700, 350));
  ScalarView mises_view("Von Mises stress [Pa]", new WinGeom(0, 405, 700, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Hermes::vector<Space<double> *>* ref_spaces = 
        Space<double>::construct_refined_spaces(Hermes::vector<Space<double> *>(&u1_space, &u2_space));
    Hermes::vector<const Space<double>*> ref_spaces_const((*ref_spaces)[0], (*ref_spaces)[1]);
    Space<double>* u_ref_space = (*ref_spaces)[0];
    Space<double>* v_ref_space = (*ref_spaces)[1];
    
    int ndof_ref = Space<double>::get_num_dofs(ref_spaces_const);

    // Initialize the FE problem.
    DiscreteProblem<double> dp(&wf, ref_spaces_const);

    // Initialize Newton solver.
    NewtonSolver<double> newton(&dp, matrix_solver);
    newton.set_verbose_output(true);

    // Time measurement.
    cpu_time.tick();

    // Perform Newton's iteration.
    info("Solving on reference mesh.");
    try
    {
      newton.solve(NULL, NEWTON_TOL, NEWTON_MAX_ITER);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    }
  
    // Time measurement.
    cpu_time.tick();

    // Translate the resulting coefficient vector into the Solution sln.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), ref_spaces_const, 
        Hermes::vector<Solution<double> *>(&u1_sln_ref, &u2_sln_ref));

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection<double>::project_global(Hermes::vector<const Space<double> *>(&u1_space, &u2_space), 
        Hermes::vector<Solution<double> *>(&u1_sln_ref, &u2_sln_ref), 
        Hermes::vector<Solution<double> *>(&u1_sln, &u2_sln), matrix_solver); 
   
    // View the coarse mesh solution and polynomial orders.
    s_view_0.show(&u1_sln); 
    o_view_0.show(&u1_space);
    s_view_1.show(&u2_sln); 
    o_view_1.show(&u2_space);
    // For von Mises stress Filter.
    double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
    double mu = E / (2*(1 + nu));
    VonMisesFilter stress(Hermes::vector<MeshFunction<double> *>(&u1_sln, &u2_sln), lambda, mu);
    mises_view.show(&stress, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u1_sln, &u2_sln, 1e3);

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Initialize adaptivity.
    Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<Space<double> *>(&u1_space, &u2_space));

    /* 
    // Register custom forms for error calculation.
    adaptivity->set_error_form(0, 0, bilinear_form_0_0<double, double>, bilinear_form_0_0<Ord, Ord>);
    adaptivity->set_error_form(0, 1, bilinear_form_0_1<double, double>, bilinear_form_0_1<Ord, Ord>);
    adaptivity->set_error_form(1, 0, bilinear_form_1_0<double, double>, bilinear_form_1_0<Ord, Ord>);
    adaptivity->set_error_form(1, 1, bilinear_form_1_1<double, double>, bilinear_form_1_1<Ord, Ord>);
    */

    // Calculate error estimate for each solution component and the total error estimate.
    info("Calculating error estimate and exact error."); 
    Hermes::vector<double> err_est_rel;
    double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution<double> *>(&u1_sln, &u2_sln), 
                               Hermes::vector<Solution<double> *>(&u1_sln_ref, &u2_sln_ref), &err_est_rel) * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    info("ndof_coarse[0]: %d, ndof_fine[0]: %d, err_est_rel[0]: %g%%", 
         u1_space.Space<double>::get_num_dofs(), Space<double>::get_num_dofs((*ref_spaces)[0]), err_est_rel[0]*100);
    info("ndof_coarse[1]: %d, ndof_fine[1]: %d, err_est_rel[1]: %g%%",
         u2_space.Space<double>::get_num_dofs(), Space<double>::get_num_dofs((*ref_spaces)[1]), err_est_rel[1]*100);
    info("ndof_coarse_total: %d, ndof_fine_total: %d, err_est_rel_total: %g%%",
         Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u1_space, &u2_space)), 
         Space<double>::get_num_dofs(ref_spaces_const), err_est_rel_total);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u1_space, &u2_space)), 
                             err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP) 
      done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector), 
                               MULTI == true ? THRESHOLD_MULTI : THRESHOLD_SINGLE, STRATEGY, MESH_REGULARITY);
    }
    if (Space<double>::get_num_dofs(Hermes::vector<const Space<double> *>(&u1_space, &u2_space)) >= NDOF_STOP) done = true;

    // Clean up.
    delete adaptivity;
    if(done == false)
      for(unsigned int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i]->get_mesh();
    delete ref_spaces;
    
    // Increase counter.
    as++;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  s_view_0.set_title("Fine mesh solution (x-displacement)");
  s_view_0.show(&u1_sln_ref);
  s_view_1.set_title("Fine mesh solution (y-displacement)");
  s_view_1.show(&u2_sln_ref);
  // For von Mises stress Filter.
  double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
  double mu = E / (2*(1 + nu));
  VonMisesFilter stress(Hermes::vector<MeshFunction<double> *>(&u1_sln_ref, &u2_sln_ref), lambda, mu);
  mises_view.show(&stress, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u1_sln_ref, &u2_sln_ref, 1e3);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

