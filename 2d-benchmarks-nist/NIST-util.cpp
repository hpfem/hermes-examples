#include "NIST-util.h"

const char* thresholds[4] = { "Low", "Medium", "High" };
const double threshold_values[4] = { 0.1, 0.5, 0.95 };

bool adaptive_step_single_space(
  Hermes::Mixins::Loggable* logger,
  MeshSharedPtr& mesh,
  SpaceSharedPtr<double>& space, 
  MeshFunctionSharedPtr<double>& sln, 
  Selector<double>& selector,
  unsigned int order_increase,
  MeshFunctionSharedPtr<double>& ref_sln,
  Hermes::Mixins::TimeMeasurable& cpu_time,
  Solver<double>& solver,
  Views::ScalarView& sview,
  Views::OrderView & oview,
  ErrorCalculator<double>& error_calculator,
  Adapt<double>& adaptivity,
  int& as,
  double error_stop,
  double& error_reached,
  int& dof_reached,
  int& cache_searches,
  int& cache_record_found,
  int& cache_record_found_reinit,
  int& cache_record_not_found,
  double& exact_error_reached,
  MeshFunctionSharedPtr<double>& exact_sln
  )
{
  try
  {
    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator refMeshCreator(mesh);
    MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh, order_increase);
    SpaceSharedPtr<double> ref_space = refSpaceCreator.create_ref_space();
    int ndof_ref = ref_space->get_num_dofs();

    logger->info("---- Adaptivity step %d (%d DOF):", as++, ndof_ref);

    solver.set_report_cache_hits_and_misses();
    solver.set_space(ref_space);
    solver.solve();

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(solver.get_sln_vector(), ref_space, ref_sln);
    
    // Project the fine mesh solution onto the coarse mesh.
    logger->info("Calculating error estimate and exact error.");
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

    // Report results - skip time.
    cpu_time.tick();
    {
      logger->info("ndof_coarse: %d, ndof_fine: %d", space->get_num_dofs(), ref_space->get_num_dofs());
      logger->info("err_est_rel: %g%%.", err_est_rel);
      if(exact_sln)
        logger->info("err_exact_rel: %g%%.", err_exact_rel);
    
      // View the coarse mesh solution and polynomial orders.
      // sview.show(ref_sln);
      // oview.show(ref_space);

      error_reached = err_est_rel;
      if(exact_sln)
        exact_error_reached = err_exact_rel;
      dof_reached = ndof_ref;
      solver.get_cache_hits_and_misses(cache_searches, cache_record_found, cache_record_found_reinit, cache_record_not_found);
    }
    cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

    if (err_est_rel < error_stop)
    {
      return true;
    }
    else
    {
      adaptivity.adapt(&selector);
      return false;
    }
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
    return true;
  }
}