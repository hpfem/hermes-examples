#include "NIST-util.h"

bool adaptive_step_single_space(
  MeshSharedPtr& mesh,
  SpaceSharedPtr<double>& space, 
  MeshFunctionSharedPtr<double>& sln, 
  MySelector& selector,
  MeshFunctionSharedPtr<double>& ref_sln,
  Hermes::Mixins::TimeMeasurable& cpu_time,
  Solver<double>& solver,
  Views::ScalarView& sview,
  Views::OrderView & oview,
  SimpleGraph graph_dof_est,
  SimpleGraph graph_cpu_est,
  ErrorCalculator<double>& error_calculator,
  Adapt<double>& adaptivity,
  int& as,
  double ERR_STOP,
  MeshFunctionSharedPtr<double>& exact_sln,
  SimpleGraph graph_dof_exact,
  SimpleGraph graph_cpu_exact
  )
{
    cpu_time.tick();

    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator refMeshCreator(mesh);
    MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh);
    SpaceSharedPtr<double> ref_space = refSpaceCreator.create_ref_space();
    int ndof_ref = ref_space->get_num_dofs();

    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d (%d DOF):", as++, ndof_ref);
    cpu_time.tick();
    
    solver.set_space(ref_space);
    solver.solve();

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(solver.get_sln_vector(), ref_space, ref_sln);
    
    cpu_time.tick();
    Hermes::Mixins::Loggable::Static::info("Solution: %g s", cpu_time.last());
    
    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate and exact error.");
    OGProjection<double>::project_global(space, ref_sln, sln);

    // Calculate element errors and total error estimate.
    double err_exact_rel;
    if(exact_sln)
    {
      error_calculator.calculate_errors(sln, exact_sln);
      err_exact_rel = error_calculator.get_total_error_squared() * 100.0;
    }
    error_calculator.calculate_errors(sln, ref_sln);
    double err_est_rel = error_calculator.get_total_error_squared() * 100.0;

    cpu_time.tick();

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d", space->get_num_dofs(), ref_space->get_num_dofs());
    Hermes::Mixins::Loggable::Static::info("err_est_rel: %g%%.", err_est_rel);
    if(exact_sln)
      Hermes::Mixins::Loggable::Static::info("err_exact_rel: %g%%.", err_exact_rel);

    // Time measurement.
    cpu_time.tick();
    double accum_time = cpu_time.accumulated();
    
    // View the coarse mesh solution and polynomial orders.
    sview.show(ref_sln);
    oview.show(ref_space);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(space->get_num_dofs(), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(accum_time, err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");
    if(exact_sln)
    {
      graph_dof_exact.add_values(space->get_num_dofs(), err_exact_rel);
      graph_dof_exact.save("conv_dof_exact.dat");
      graph_cpu_exact.add_values(accum_time, err_exact_rel);
      graph_cpu_exact.save("conv_cpu_exact.dat");
    }
    cpu_time.tick();

    // If err_est too large, adapt the mesh. The NDOF test must be here, so that the solution may be visualized
    // after ending due to this criterion.
    if (err_est_rel < ERR_STOP) 
      return true;
    else
      return adaptivity.adapt(&selector);
   
  }