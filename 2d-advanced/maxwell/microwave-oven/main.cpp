#include "definitions.h"

// This example solves adaptively the electric field in a simplified microwave oven.
// The waves are generated using a harmonic surface current on the right-most edge.
// (Such small cavity is present in every microwave oven). There is a circular
// load located in the middle of the main cavity, defined through a different
// permittivity -- see function in_load(...). One can either use a mesh that is
// aligned to the load via curvilinear elements (ALIGN_MESH = true), or an unaligned
// mesh (ALIGN_MESH = false). Convergence graphs are saved both wrt. the dof number
// and cpu time.
//
// PDE: time-harmonic Maxwell's equations;
//      there is circular load in the middle of the large cavity, whose permittivity
//      is different from the rest of the domain.
//
// Domain: square cavity with another small square cavity attached from outside
//         on the right.
//
// Meshes: you can either use "oven_load_circle.mesh" containing curved elements
//         aligned with the circular load, or "oven_load_square.mesh" which is not
//         aligned.
//
// BC: perfect conductor on the boundary except for the right-most edge of the small
//     cavity, where a harmonic surface current is prescribed
//
// The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;
// Initial polynomial degree. NOTE: The meaning is different from
// standard continuous elements in the space H1. Here, P_INIT refers
// to the maximum poly order of the tangential component, and polynomials
// of degree P_INIT + 1 are present in element interiors. P_INIT = 0
// is for Whitney elements.
const int P_INIT = 1;
// if ALIGN_MESH == true, curvilinear elements aligned with the
// circular load are used, otherwise one uses a non-aligned mesh->
const bool ALIGN_MESH = false;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.6;
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
const double ERR_STOP = 5e-2;
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;

// Newton's method.
const double NEWTON_TOL = 1e-4;
const int NEWTON_MAX_ITER = 100;

// Problem parameters.
const double e_0 = 8.8541878176 * 1e-12;
const double mu_0 = 1.256 * 1e-6;
const double e_r = 1.0;
const double mu_r = 1.0;
const double rho = 3820.0;
const double Cp = 7.531000;
const double freq = 1.0*2450000000.0;
const double omega = 2 * M_PI * freq;
const double c = 1 / std::sqrt(e_0 * mu_0);
const double kappa = 2 * M_PI * freq * std::sqrt(e_0 * mu_0);
const double J = 0.0000033333;

//  Boundary markers.
const std::string BDY_PERFECT_CONDUCTOR = "b2";
const std::string BDY_CURRENT = "b1";

class CustomErrorCalculator : public ErrorCalculator < complex >
{
public:
  CustomErrorCalculator(CalculatedErrorType errorType) : ErrorCalculator<::complex>(errorType)
  {
    this->add_error_form(new CustomErrorForm(kappa));
  }
  virtual ~CustomErrorCalculator()
  {
  }
};

int main(int argc, char* argv[])
{
  // Load the mesh->
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  if (ALIGN_MESH)
    mloader.load("oven_load_circle.mesh", mesh);
  else
    mloader.load("oven_load_square.mesh", mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize boundary conditions
  DefaultEssentialBCConst<::complex> bc_essential(BDY_PERFECT_CONDUCTOR, std::complex<double>(0.0, 0.0));

  EssentialBCs<::complex> bcs(&bc_essential);

  // Create an Hcurl space with default shapeset.
  SpaceSharedPtr<::complex> space(new HcurlSpace<::complex>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();
  Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf(new CustomWeakForm(e_0, mu_0, mu_r, kappa, omega, J, ALIGN_MESH, mesh, BDY_CURRENT));

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<::complex>  sln(new Solution<::complex>), ref_sln(new Solution<::complex>);

  // Initialize refinements selector.
  HcurlProjBasedSelector<::complex> selector(CAND_LIST);

  // Error calculation.
  CustomErrorCalculator errorCalculator(RelativeErrorToGlobalNorm);

  // Stopping criterion for an adaptivity step.
  AdaptStoppingCriterionSingleElement<::complex> stoppingCriterion(THRESHOLD);

  // Adaptivity.
  Adapt<::complex> adaptivity(space, &errorCalculator, &stoppingCriterion);

  // Initialize views.
  ScalarView eview("Electric field", new WinGeom(0, 0, 580, 400));
  OrderView  oview("Polynomial orders", new WinGeom(590, 0, 550, 400));

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  Hermes::Hermes2D::NewtonSolver<::complex> newton(wf, space);
  newton.set_max_allowed_iterations(NEWTON_MAX_ITER);
  newton.set_tolerance(NEWTON_TOL, Hermes::Solvers::ResidualNormAbsolute);

  // Newton's iteration.
  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator refMeshCreator(mesh);
    MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

    Space<::complex>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh);
    SpaceSharedPtr<::complex> ref_space = refSpaceCreator.create_ref_space();
    int ndof_ref = Space<::complex>::get_num_dofs(ref_space);

    // Initialize reference problem.
    Hermes::Mixins::Loggable::Static::info("Solving on reference mesh->");
    // Time measurement.
    cpu_time.tick();
    newton.set_space(ref_space);

    try
    {
      newton.solve();
    }
    catch (Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    };
    // Translate the resulting coefficient vector into the Solution<::complex> sln->
    Hermes::Hermes2D::Solution<::complex>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh->
    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh->");
    OGProjection<::complex> ogProjection; ogProjection.project_global(space, ref_sln, sln);

    // View the coarse mesh solution and polynomial orders.
    MeshFunctionSharedPtr<double> real(new RealFilter(ref_sln));
    MeshFunctionSharedPtr<double> magn(new MagFilter<double>(real));
    MeshFunctionSharedPtr<double> limited_magn(new ValFilter(magn, 0.0, 4e3));
    char title[100];
    sprintf(title, "Electric field, adaptivity step %d", as);
    eview.set_title(title);
    //eview.set_min_max_range(0.0, 4e3);
    eview.show(limited_magn, HERMES_EPS_LOW);
    sprintf(title, "Polynomial orders, adaptivity step %d", as);
    oview.set_title(title);
    oview.show(space);

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");

    // Calculate error estimate.
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100.;

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      Space<::complex>::get_num_dofs(space),
      Space<::complex>::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space<::complex>::get_num_dofs(space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh->
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh->");
      done = adaptivity.adapt(&selector);
    }
    if (space->get_num_dofs() >= NDOF_STOP) done = true;

    // Increase counter.
    as++;
  } while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}