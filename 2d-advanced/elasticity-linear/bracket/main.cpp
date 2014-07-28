#include "definitions.h"

// This example uses adaptive multimesh hp-FEM to solve a simple problem
// of linear elasticity. Note that since both displacement components
// have similar qualitative behavior, the advantage of the multimesh
// discretization is less striking than for example in the tutorial
// example P04-linear-adapt/02-system-adapt.
//
// PDE: Lame equations of linear elasticity, treated as a coupled system
//      of two PDEs.
//
// BC: u_1 = u_2 = 0 on Gamma_1 (right edge)
//     du_1/dn = f0 on Gamma_2 (top edge)
//     du_2/dn = f1 on Gamma_2 (top edge)
//     du_1/dn = du_2/dn = 0 elsewhere
//
// The following parameters can be changed:

// Initial polynomial degree of all mesh elements.
const int P_INIT = 2;
// MULTI = true  ... use multi-mesh,
// MULTI = false ... use single-mesh->
// Note: In the single mesh option, the meshes are
// forced to be geometrically the same but the
// polynomial degrees can still vary.
const bool MULTI = true;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD_MULTI = 0.35;
const double THRESHOLD_SINGLE = 0.7;
// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 2);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(MULTI ? THRESHOLD_MULTI : THRESHOLD_SINGLE);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;

// Problem parameters.
// Young modulus for steel: 200 GPa.
const double E = 200e9;
// Poisson ratio.
const double nu = 0.3;
// Density.
const double rho = 8000.0;
// Gravitational acceleration.
const double g1 = -9.81;
// Top surface force in x-direction.
const double f0 = 0;
// Top surface force in y-direction.
const double f1 = -1e3;

int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load the mesh->
  MeshSharedPtr u1_mesh(new Mesh), u2_mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", u1_mesh);

  // Initial mesh refinements.
  u1_mesh->refine_element_id(1);
  u1_mesh->refine_element_id(4);

  // Create initial mesh for the vertical displacement component.
  // This also initializes the multimesh hp-FEM.
  u2_mesh->copy(u1_mesh);

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> zero_disp("bdy_right", 0.0);
  EssentialBCs<double> bcs(&zero_disp);

  // Create x- and y- displacement space using the default H1 shapeset.
  SpaceSharedPtr<double> u1_space(new H1Space<double>(u1_mesh, &bcs, P_INIT));
  SpaceSharedPtr<double> u2_space(new H1Space<double>(u2_mesh, &bcs, P_INIT));
  Hermes::Mixins::Loggable::Static::info("ndof = %d.", Space<double>::get_num_dofs({u1_space, u2_space}));

  // Initialize the weak formulation.
  // NOTE; These weak forms are identical to those in example P01-linear/08-system.
  CustomWeakFormLinearElasticity wf(E, nu, rho*g1, "bdy_top", f0, f1);

  // Initialize the FE problem.
  std::vector<SpaceSharedPtr<double> > spaces(u1_space, u2_space);
  DiscreteProblem<double> dp(wf, spaces);

  // Initialize coarse and reference mesh solutions.
  MeshFunctionSharedPtr<double> u1_sln(new Solution<double>), u2_sln(new Solution<double>);
  MeshFunctionSharedPtr<double>u1_sln_ref(new Solution<double>), u2_sln_ref(new Solution<double>);

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  // Initialize views.
  ScalarView s_view_0("Solution (x-displacement)", new WinGeom(0, 0, 500, 350));
  s_view_0.show_mesh(false);
  ScalarView s_view_1("Solution (y-displacement)", new WinGeom(760, 0, 500, 350));
  s_view_1.show_mesh(false);
  OrderView  o_view_0("Mesh (x-displacement)", new WinGeom(410, 0, 440, 350));
  OrderView  o_view_1("Mesh (y-displacement)", new WinGeom(1170, 0, 440, 350));
  ScalarView mises_view("Von Mises stress [Pa]", new WinGeom(0, 405, 500, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator refMeshCreator1(u1_mesh);
    MeshSharedPtr ref_u1_mesh = refMeshCreator1.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreator1(u1_space, ref_u1_mesh);
    SpaceSharedPtr<double> ref_u1_space = refSpaceCreator1.create_ref_space();

    Mesh::ReferenceMeshCreator refMeshCreator2(u2_mesh);
    MeshSharedPtr ref_u2_mesh = refMeshCreator2.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreator2(u2_space, ref_u2_mesh);
    SpaceSharedPtr<double> ref_u2_space = refSpaceCreator2.create_ref_space();

    std::vector<SpaceSharedPtr<double> > ref_spaces(ref_u1_space, ref_u2_space);

    int ndof_ref = Space<double>::get_num_dofs(ref_spaces);

    // Initialize the FE problem.
    DiscreteProblem<double> dp(wf, ref_spaces);

    // Initialize Newton solver.
    NewtonSolver<double> newton(&dp);
    newton.set_verbose_output(true);

    // Time measurement.
    cpu_time.tick();

    // Perform Newton's iteration.
    Hermes::Mixins::Loggable::Static::info("Solving on reference mesh->");
    try
    {
      newton.solve();
    }
    catch (Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    }

    // Time measurement.
    cpu_time.tick();

    // Translate the resulting coefficient vector into the Solution sln->
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), ref_spaces,
      std::vector<MeshFunctionSharedPtr<double> >(u1_sln_ref, u2_sln_ref));

    // Project the fine mesh solution onto the coarse mesh->
    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh->");
    OGProjection<double> ogProjection; ogProjection.project_global({u1_space, u2_space},
      std::vector<MeshFunctionSharedPtr<double> >(u1_sln_ref, u2_sln_ref),
      std::vector<MeshFunctionSharedPtr<double> >(u1_sln, u2_sln));

    // View the coarse mesh solution and polynomial orders.
    s_view_0.show(u1_sln);
    o_view_0.show(u1_space);
    s_view_1.show(u2_sln);
    o_view_1.show(u2_space);
    // For von Mises stress Filter.
    double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));
    MeshFunctionSharedPtr<double> stress(new VonMisesFilter({u1_sln, u2_sln}, lambda, mu));
    MeshFunctionSharedPtr<double> limited_stress(new ValFilter(stress, 0.0, 2e5));
    mises_view.show(limited_stress, H2D_FN_VAL_0, u1_sln, u2_sln, 0);

    // Skip visualization time.
    cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

    // Initialize adaptivity.
    adaptivity.set_spaces({u1_space, u2_space});

    /*
    // Register custom forms for error calculation.
    errorCalculator.add_error_form(0, 0, bilinear_form_0_0<double, double>, bilinear_form_0_0<Ord, Ord>);
    errorCalculator.add_error_form(0, 1, bilinear_form_0_1<double, double>, bilinear_form_0_1<Ord, Ord>);
    errorCalculator.add_error_form(1, 0, bilinear_form_1_0<double, double>, bilinear_form_1_0<Ord, Ord>);
    errorCalculator.add_error_form(1, 1, bilinear_form_1_1<double, double>, bilinear_form_1_1<Ord, Ord>);
    */

    // Calculate error estimate for each solution component and the total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate and exact error.");
    errorCalculator.calculate_errors({u1_sln, u2_sln},
      std::vector<MeshFunctionSharedPtr<double> >(u1_sln_ref, u2_sln_ref));
    double err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse_total: %d, ndof_fine_total: %d, err_est_rel_total: %g%%",
      Space<double>::get_num_dofs({u1_space, u2_space}),
      Space<double>::get_num_dofs(ref_spaces), err_est_rel_total);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs({u1_space, u2_space}),
      err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh->
    if (err_est_rel_total < ERR_STOP)
      done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh->");
      done = adaptivity.adapt({&selector, &selector});
    }

    // Increase counter.
    as++;
  } while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  s_view_0.set_title("Fine mesh solution (x-displacement)");
  s_view_0.show(u1_sln_ref);
  s_view_1.set_title("Fine mesh solution (y-displacement)");
  s_view_1.show(u2_sln_ref);
  // For von Mises stress Filter.
  double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
  double mu = E / (2 * (1 + nu));
  MeshFunctionSharedPtr<double> stress(new VonMisesFilter({u1_sln_ref, u2_sln_ref}, lambda, mu));
  MeshFunctionSharedPtr<double> limited_stress(new ValFilter(stress, 0.0, 2e5));
  mises_view.show(limited_stress, H2D_FN_VAL_0, u1_sln_ref, u2_sln_ref, 0);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}