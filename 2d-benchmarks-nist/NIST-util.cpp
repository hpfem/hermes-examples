#include "NIST-util.h"

const char* thresholds[7] = { "Lowest", "Lower", "Low", "Medium", "High", "Higher", "Highest" };
const double threshold_values[7] = { 0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 0.95 };
extern const char* strategies[6] = {
  "noSelectionH",
  "noSelectionHP",
  "hXORpError",
  "hORpDOFs",
  "isoHPDOFs",
  "anisoHPDOFs"
};

char* process_arguments_main_comparison(int argc, char* argv[], Selector<double>*& selector, AdaptivityStoppingCriterion<double>*& stoppingCriterion)
{
  assert(argc > 2);

  char* toReturn = new char[1000];

  int selectorIndex = atoi(argv[1]);

  selector = new MySelector((hpAdaptivityStrategy)selectorIndex);

  int stoppingCriterionIndex = atoi(argv[2]);

  stoppingCriterion = new AdaptStoppingCriterionSingleElement<double>(threshold_values[stoppingCriterionIndex]);

  sprintf(toReturn, "%s-%s", strategies[selectorIndex], thresholds[stoppingCriterionIndex]);

  return toReturn;
}

bool adaptive_step_single_space(
  Hermes::Mixins::Loggable* logger,
  MeshSharedPtr& mesh,
  SpaceSharedPtr<double>& space, 
  MeshFunctionSharedPtr<double>& sln, 
  Selector<double>* selector,
  unsigned int order_increase,
  MeshFunctionSharedPtr<double>& ref_sln,
  Hermes::Mixins::TimeMeasurable& cpu_time,
  Solver<double>& solver,
  Views::ScalarView& sview,
  Views::OrderView & oview,
  ErrorCalculator<double>& errorCalculator,
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
  double& FactorizationSize,
  double& PeakMemoryUsage,
  double& Flops,
  MeshFunctionSharedPtr<double>& exact_sln
  )
{
  try
  {
    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator refMeshCreator(mesh);
    MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh, order_increase);
    SpaceSharedPtr<double> ref_space = refSpaceCreator.create_ref_space();
    int ndof_ref = ref_space->get_num_dofs();

    logger->info("---- Adaptivity step %d (%d DOF):", as++, ndof_ref);

    solver.set_report_cache_hits_and_misses();
    solver.set_space(ref_space);
    solver.solve();
    FactorizationSize = solver.get_UMFPACK_reporting_data(Solver<double>::FactorizationSize);
    PeakMemoryUsage = solver.get_UMFPACK_reporting_data(Solver<double>::PeakMemoryUsage);
    Flops = solver.get_UMFPACK_reporting_data(Solver<double>::Flops);

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(solver.get_sln_vector(), ref_space, ref_sln);
    
    // Project the fine mesh solution onto the coarse mesh.
    logger->info("Calculating error estimate and exact error.");
    OGProjection<double>::project_global(space, ref_sln, sln);

    // Calculate element errors and total error estimate.
    double err_exact_rel;
    if(exact_sln)
    {
      errorCalculator.calculate_errors(sln, exact_sln, false);
      err_exact_rel = errorCalculator.get_total_error_squared() * 100.0;
    }
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100.0;

    // Report results - skip time.
    cpu_time.tick();
    {
      logger->info("ndof_coarse: %d, ndof_fine: %d", space->get_num_dofs(), ref_space->get_num_dofs());
      logger->info("err_est_rel: %g%%.", err_est_rel);
      if(exact_sln)
        logger->info("err_exact_rel: %g%%.", err_exact_rel);
    
      // View the coarse mesh solution and polynomial orders.
       sview.show(ref_sln);
       oview.show(ref_space);

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
      adaptivity.adapt(selector);
      return false;
    }
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
    throw;
    return true;
  }
}