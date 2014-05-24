#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "function/function.h"

// Read original or XML mesh file.
const bool USE_XML_FORMAT = true;
// Visualisation of the distortion and stress
const bool VIEW_DIST = false;
const bool VIEW_STRESS = false;
// Parameter influencing the candidate selection.
const double THRESHOLD = .8;
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_H_ISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1.0;
// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 3);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Selector.
H1ProjBasedSelector<double> selector(CAND_LIST);

// Initial polynomial degree of all elements.
const int P_INIT_U = 2;
const int P_INIT_V = 2;
const int P_INIT_P = 1;

// Newton's method.
const double NEWTON_TOL_FINE = 1e-0;
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-3;
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 25;

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 3;
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;

// Problem parameters.
// Young modulus (cautchuque).
const double E  = 1000;
// Poisson ratio.
const double nu = 0.49;
// Density.
const double rho = 1010.0;
// Gravitational acceleration.
const double g1 = -9.81;
// Surface force in x-direction.
const double f0  = 0;
// Surface force in y-direction.
const double f1  = 8e4;
// Components of distortion
const double Eo11 = 0.0;
const double Eo12 = 0.0;
const double Eo22 = 0.1;
const double Eo33 = 0.0;
CustomEo11 Eo11d(Eo11);
CustomEo12 Eo12d(Eo12);
CustomEo22 Eo22d(Eo22);
CustomEo33 Eo33d(Eo33);

int main(int argc, char* argv[])
{
  // Time measurement
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  if (USE_XML_FORMAT == true)
  {
    MeshReaderH2DXML mloader;
    Hermes::Mixins::Loggable::Static::info("Reading mesh in XML format.");
    mloader.load("rect_1-1_Q.xml", mesh);
  }
  else
  {
    MeshReaderH2D mloader;
    Hermes::Mixins::Loggable::Static::info("Reading mesh in original format.");
    mloader.load("rect_1-1_Q.mesh", mesh);
  }

  // Perform uniform mesh refinement.
  int refinement_type = 0;
  for (int i = 0; i < INIT_REF_NUM; i++)
  {
    mesh->refine_all_elements(refinement_type);
  }

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> disp_bot_top_x(Hermes::vector<std::string>("Bottom","Top"), 0.0);
  DefaultEssentialBCConst<double> disp_bot_y("Bottom", 0.0);
  DefaultEssentialBCConst<double> disp_top_y("Top", 0.0);
  EssentialBCs<double> bcs_x(&disp_bot_top_x);
  EssentialBCs<double> bcs_y(Hermes::vector<EssentialBoundaryCondition<double> *>(&disp_bot_y, &disp_top_y));

  // Create x- and y- displacement space using the default H1 shapeset.
  SpaceSharedPtr<double> u_space(new H1Space<double>(mesh, &bcs_x, P_INIT_U));
  SpaceSharedPtr<double> v_space(new H1Space<double>(mesh, &bcs_y, P_INIT_V));
  SpaceSharedPtr<double> p_space(new H1Space<double>(mesh, P_INIT_P));
  Hermes::vector<SpaceSharedPtr<double>> spaces(u_space, v_space, p_space);
  // Set the spaces to adaptivity.
  adaptivity.set_spaces(spaces);

  // Initialise coarse and reference mesh solutions
  MeshFunctionSharedPtr<double> u_sln(new Solution<double>), v_sln(new Solution<double>), p_sln(new Solution<double>),
    u_ref_sln(new Solution<double>), v_ref_sln(new Solution<double>), p_ref_sln(new Solution<double>);

  Hermes::vector<MeshFunctionSharedPtr<double>> solutions(u_sln, v_sln, p_sln);
  Hermes::vector<MeshFunctionSharedPtr<double>> ref_solutions(u_ref_sln, v_ref_sln, p_ref_sln);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;
  SimpleGraph graph_dof_exact, graph_cpu_exact;

  // Initialize the weak formulation.
  CustomWeakFormLinearElasticity wf(E, nu, &Eo11d, &Eo12d, &Eo22d, &Eo33d, Eo11, Eo12, Eo22, Eo33, rho*g1, "Top", f0, f1);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, spaces);
  // Initialize Newton solver.
  NewtonSolver<double> newton(&dp);
  newton.set_verbose_output(true);
  newton.set_tolerance(1e-3, Hermes::Solvers::ResidualNormAbsolute);

  // Initialize views.
  Views::OrderView o_view_0("Mesh u", new Views::WinGeom(0, 0, 420, 350));
  Views::OrderView o_view_1("Mesh v", new Views::WinGeom(450, 0, 420, 350));
  Views::OrderView o_view_2("Mesh p", new Views::WinGeom(1330, 0, 420, 350));
  Views::ScalarView s_view_0("Solution u [m]", new Views::WinGeom(0, 380, 440, 350));
  Views::ScalarView s_view_1("Solution v [m]", new Views::WinGeom(470, 380, 440, 350));
  Views::ScalarView s_view_2("Solution p [Pa]", new Views::WinGeom(940, 380, 440, 350));
  s_view_0.show_mesh(false);
  s_view_1.show_mesh(false);
  s_view_2.show_mesh(false);

  Views::ScalarView viewEo22("Distortion Eo22 [1]", new Views::WinGeom(690, 710, 680, 400));
  Views::ScalarView viewTrEo("Distortion Tr(Eo) [1]", new Views::WinGeom(690, 710, 680, 400));
  viewEo22.show_mesh(false);
  viewTrEo.show_mesh(false);

  Views::ScalarView viewU("Displacement u [m]", new Views::WinGeom(0, 710, 680, 400));
  Views::ScalarView viewV("Displacement v [m]", new Views::WinGeom(690, 710, 680, 400));
  Views::ScalarView viewP("Pressure p [Pa]", new Views::WinGeom(390, 410, 680, 400));
  viewU.show_mesh(false);
  viewV.show_mesh(false);
  viewP.show_mesh(false);

  Views::ScalarView viewS11("Stress S11 [Pa]", new Views::WinGeom(20, 40, 680, 400));
  Views::ScalarView viewS12("Stress S12 [Pa]", new Views::WinGeom(120, 80, 680, 400));
  Views::ScalarView viewS22("Stress S22 [Pa]", new Views::WinGeom(220, 120, 680, 400));
  Views::ScalarView viewS33("Stress S33 [Pa]", new Views::WinGeom(320, 160, 680, 400));
  Views::ScalarView view_vM("Stress von Mises [Pa]", new Views::WinGeom(420, 200, 680, 400));
  viewS11.show_mesh(false);
  viewS12.show_mesh(false);
  viewS22.show_mesh(false);
  viewS33.show_mesh(false);
  view_vM.show_mesh(false);

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info(" ---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator u_ref_space_creator(u_space, ref_mesh);
    SpaceSharedPtr<double> u_ref_space = u_ref_space_creator.create_ref_space();
    Space<double>::ReferenceSpaceCreator v_ref_space_creator(v_space, ref_mesh);
    SpaceSharedPtr<double> v_ref_space = v_ref_space_creator.create_ref_space();
    Space<double>::ReferenceSpaceCreator p_ref_space_creator(p_space, ref_mesh);
    SpaceSharedPtr<double> p_ref_space = p_ref_space_creator.create_ref_space();
/*
    MeshView mvu("Mesh u", new WinGeom(0, 0, 580, 400));
    mvu.show(u_mesh);
    MeshView mvv("Mesh v", new WinGeom(500, 0, 580, 400));
    mvv.show(v_mesh);
    MeshView mvp("Mesh p", new WinGeom(800, 0, 580, 400));
    mvp.show(p_mesh);

    MeshView mvru("Mesh ref u", new WinGeom(0, 600, 580, 400));
    mvru.show(u_ref_mesh);
    MeshView mvrv("Mesh ref v", new WinGeom(500, 600, 580, 400));
    mvrv.show(v_ref_mesh);
    MeshView mvrp("Mesh ref p", new WinGeom(800, 600, 580, 400));
    mvrp.show(p_ref_mesh);
*/
    Hermes::vector<SpaceSharedPtr<double>> ref_spaces(u_ref_space, v_ref_space, p_ref_space);

    int ndof_ref = Space<double>::get_num_dofs(ref_spaces);

    // Initialize reference problem.
    Hermes::Mixins::Loggable::Static::info("Solving on reference mesh.");

    // Time measurement.
    cpu_time.tick();

    // Perform Newton's iteration.
    try
    {
      newton.set_spaces(ref_spaces);
      newton.solve();
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), ref_spaces, ref_solutions);

    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
    OGProjection<double>::project_global(spaces, ref_solutions, solutions);

    cpu_time.tick();
/*
    // View the coarse mesh solution and polynomial orders.
    s_view_0.show(u_sln);
    s_view_1.show(v_sln);
    s_view_2.show(p_sln);
    o_view_0.show(u_space);
    o_view_1.show(v_space);
    o_view_2.show(p_space);
*/
    // View the refined mesh solution and polynomial orders.
    s_view_0.show(u_ref_sln);
    s_view_1.show(v_ref_sln);
    s_view_2.show(p_ref_sln);
    o_view_0.show(u_ref_space);
    o_view_1.show(v_ref_space);
    o_view_2.show(p_ref_space);

    // Calculate error estimate.
    errorCalculator.calculate_errors(solutions, ref_solutions);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse[0]: %d, ndof_fine[0]: %d",
         u_space->get_num_dofs(), u_ref_space->get_num_dofs());
    Hermes::Mixins::Loggable::Static::info("ndof_coarse[1]: %d, ndof_fine[1]: %d",
         v_space->get_num_dofs(), v_ref_space->get_num_dofs());
    Hermes::Mixins::Loggable::Static::info("ndof_coarse[2]: %d, ndof_fine[2]: %d",
         p_space->get_num_dofs(), p_ref_space->get_num_dofs());
    Hermes::Mixins::Loggable::Static::info("ndof_coarse_total: %d, ndof_fine_total: %d",
         Space<double>::get_num_dofs(spaces), Space<double>::get_num_dofs(ref_spaces));
    Hermes::Mixins::Loggable::Static::info("err_est_rel_total: %g%%", err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs(spaces), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP)
      done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(Hermes::vector<RefinementSelectors::Selector<double>*>(&selector, &selector, &selector));
    }

    // Increase counter.
    as++;
  }
  while (as<2); // ! ! ! attention this is to control the number of iterations to be replaced by the line bottom
//  while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Distortion visualisation
  if (VIEW_DIST == true)
  {
    MeshFunctionSharedPtr<double> exact_Eo22(new ExactSolutionEo22(mesh, Eo22));
    viewEo22.show(exact_Eo22);

    MeshFunctionSharedPtr<double> exact_TrEo(new ExactSolutionTrEo(mesh, Eo11, Eo22, Eo33));
    viewTrEo.show(exact_TrEo);
  }

  // Stress visualisation
  if (VIEW_STRESS == true)
    {
    // First Lame constant.
    double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
    // Second Lame constant.
    double mu = E / (2 * (1 + nu));

    // Visualize the solution.
    MeshFunctionSharedPtr<double> S11(new CustomFilterS11(solutions, lambda, mu, &Eo11d, &Eo12d, &Eo22d, &Eo33d));
    MeshFunctionSharedPtr<double> S12(new CustomFilterS12(solutions, lambda, mu, &Eo11d, &Eo12d, &Eo22d, &Eo33d));
    MeshFunctionSharedPtr<double> S22(new CustomFilterS22(solutions, lambda, mu, &Eo11d, &Eo12d, &Eo22d, &Eo33d));
    MeshFunctionSharedPtr<double> S33(new CustomFilterS33(solutions, lambda, mu, &Eo11d, &Eo12d, &Eo22d, &Eo33d));
    MeshFunctionSharedPtr<double> vM(new CustomFilter_vM(solutions, lambda, mu, &Eo11d, &Eo12d, &Eo22d, &Eo33d));

    viewU.show(u_sln);
    viewV.show(v_sln);
    viewP.show(p_sln);
    viewS11.show(S11, HERMES_EPS_HIGH, H2D_FN_VAL_0, u_sln, v_sln, 1.0);
    viewS12.show(S12, HERMES_EPS_HIGH, H2D_FN_VAL_0, u_sln, v_sln, 1.0);
    viewS22.show(S22, HERMES_EPS_HIGH, H2D_FN_VAL_0, u_sln, v_sln, 1.0);
    viewS33.show(S33, HERMES_EPS_HIGH, H2D_FN_VAL_0, u_sln, v_sln, 1.0);
    view_vM.show(vM, HERMES_EPS_HIGH, H2D_FN_VAL_0, u_sln, v_sln, 1.0);
  }

  // Wait for the view to be closed.
  View::wait();

  return 0;
}
