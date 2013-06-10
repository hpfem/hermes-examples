#include "definitions.h"

//  This problem describes the distribution of the vector potential in
//  a 2D domain comprising a wire carrying electrical current, air, and
//  an iron which is not under voltage.
//
//  PDE: -div(1/rho grad p) - omega**2 / (rho c**2) * p = 0.
//
//  Domain: Floor plan of an existing apartment.
//
//  BC: Prescribed pressure at the source.
//      Newton matched boundary on the walls.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Initial polynomial degree of mesh elements.
const int P_INIT = 2;                             
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.75;                     

// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;

// Problem parameters.
const double RHO = 1.25;
const double FREQ = 5e2;
const double OMEGA = 2 * M_PI * FREQ;
const double SOUND_SPEED = 353.0;
const std::complex<double> P_SOURCE(1.0, 0.0);

int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Error calculation & adaptivity.
  DefaultErrorCalculator<complex, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
  // Stopping criterion for an adaptivity step.
  AdaptStoppingCriterionSingleElement<complex> stoppingCriterion(THRESHOLD);
  // Adaptivity processor class.
  Adapt<complex> adaptivity(&errorCalculator, &stoppingCriterion);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<complex> bc_essential("Source", P_SOURCE);
  EssentialBCs<complex> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<complex> space(new H1Space<complex> (mesh, &bcs, P_INIT));
  int ndof = Space<complex>::get_num_dofs(space);
  Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakFormAcoustics wf("Wall", RHO, SOUND_SPEED, OMEGA);

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<complex>  sln(new Solution<complex>), ref_sln(new Solution<complex>);

  // Initialize refinement selector.
  H1ProjBasedSelector<complex> selector(CAND_LIST);

  // Initialize views.
  ScalarView sview("Acoustic pressure", new WinGeom(0, 0, 600, 350));
  ScalarView eview("Error", new WinGeom(0, 370, 600, 350));
  sview.show_mesh(false);
  sview.fix_scale_width(50);
  OrderView  oview("Polynomial orders", new WinGeom(610, 0, 600, 350));
  ScalarView ref_view("Refined elements", new WinGeom(610, 370, 600, 350));

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  // Adaptivity loop:
  int as = 1;
  adaptivity.set_space(space);
  bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator refMeshCreator(mesh);
    MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

    Space<complex>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh);
    SpaceSharedPtr<complex> ref_space = refSpaceCreator.create_ref_space();

    int ndof_ref = Space<complex>::get_num_dofs(ref_space);

    // Assemble the reference problem.
    Hermes::Mixins::Loggable::Static::info("Solving on reference mesh.");
    DiscreteProblem<complex> dp(&wf, ref_space);

    // Time measurement.
    cpu_time.tick();

    // Perform Newton's iteration.
    Hermes::Hermes2D::NewtonSolver<complex> newton(&dp);
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    };
    // Translate the resulting coefficient vector into the Solution<complex> sln->
    Hermes::Hermes2D::Solution<complex>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
    OGProjection<complex> ogProjection; ogProjection.project_global(space, ref_sln, sln);

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution and polynomial orders.
    MeshFunctionSharedPtr<double> acoustic_pressure(new RealFilter(sln));
    sview.show(acoustic_pressure);
    oview.show(space);

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    eview.show(errorCalculator.get_errorMeshFunction());

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      Space<complex>::get_num_dofs(space), Space<complex>::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space<complex>::get_num_dofs(space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(&selector);
      ref_view.show(adaptivity.get_refinementInfoMeshFunction());
      cpu_time.tick();
      char* cpu_time_string = cpu_time.accumulated_str().c_str();
      ref_view.titl
      ref_view.set_title();
    }
    
    // Increase counter.
    as++;
  }
  while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution magnitude");

  MeshFunctionSharedPtr<double>  ref_mag(new RealFilter(ref_sln));
  sview.show(ref_mag);

  // Output solution in VTK format.
  Linearizer lin;
  bool mode_3D = true;
  lin.save_solution_vtk(ref_mag, "sln.vtk", "Acoustic pressure", mode_3D);
  Hermes::Mixins::Loggable::Static::info("Solution in VTK format saved to file %s.", "sln.vtk");

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
