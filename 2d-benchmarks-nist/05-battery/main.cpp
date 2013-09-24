
#include "definitions.h"

using namespace RefinementSelectors;

//  This is the fifth in the series of NIST benchmarks with unknown exact solution.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: -\frac{\partial }{\partial x}\left(p(x, y)\frac{\partial u}{\partial x}\right)
//       -\frac{\partial }{\partial y}\left(q(x, y)\frac{\partial u}{\partial y}\right) - f = 0.
//
//  Exact solution: unknown.
//
//  Domain: square (0, 8.4) x (0, 24), see the file "battery.mesh".
//
//  BC: Zero Neumann on left edge, Newton on the rest of the boundary:
//
//  The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// This is a quantitative parameter of Adaptivity.
double THRESHOLD = 0.3;
// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;

// Newton tolerance
const double NEWTON_TOLERANCE = 1e-6;

int main(int argc, char* argv[])
{
  const char* THRESHOLD_STRING = "Custom";

  if(argc > 2)
  {
    CAND_LIST = CandList(atoi(argv[1]));
    THRESHOLD = threshold_values[atoi(argv[2])];
    THRESHOLD_STRING = thresholds[atoi(argv[2])];
  }

  sprintf(Hermes::Mixins::Loggable::staticLogFileName, "Logfile-%s-%s.log", get_cand_list_str(CAND_LIST), THRESHOLD_STRING);

  double ERR_STOP = 1.;


  // Load the mesh.
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("battery.mesh", mesh);
  MeshView m;
  m.show(mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  basemesh->copy(mesh);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, P_INIT));

  // Initialize weak formulation.
  CustomWeakFormPoisson wf("e1", "e2", "e3", "e4", "e5", 
                           "Bdy_left", "Bdy_top", "Bdy_right", "Bdy_bottom", mesh);

  // Initialize coarse and fine mesh solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>), ref_sln(new Solution<double>);
  
  // Initialize refinement selector.
  MySelector selector(hXORpSelectionBasedOnError);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 320, 600));
  OrderView  oview("Polynomial orders", new WinGeom(330, 0, 300, 600));
  Linearizer lin;
  Orderizer ord;
  char* filename = new char[1000];

  // Assemble the discrete problem.    
  NewtonSolver<double> newton;
  newton.set_weak_formulation(&wf);

  sprintf(filename, "Results-%s-%s.csv", get_cand_list_str(CAND_LIST), THRESHOLD_STRING);
  std::ofstream data(filename);
  data.precision(10);
  data.setf( std::ios::fixed, std::ios::floatfield );
  data << 
      "Iteration" << ';' <<
      "CPUTime" << ';' <<
      "AdaptivitySteps" << ';' <<
      "dof_reached" << ';' <<
      "dof_cumulative" << ';' <<
      "total_cache_searches" << ';' <<
      "total_cache_record_found" << ';' <<
      "total_cache_record_found_reinit" << ';' <<
      "total_cache_record_not_found" << ';' <<
      "error_stop" << ';' <<
      "error_reached" << ';' <<
      "exact_error_reached" <<
      std::endl;

  Hermes::Mixins::Loggable logger;
  logger.set_verbose_output(true);

  int iterations_count = 8;
  int error_levels_count = 5;
  double error_stop = ERR_STOP;

  for(int iteration = 0; iteration < iterations_count; iteration++)
  {
    for(int error_level = 0; error_level < error_levels_count; error_level++)
    {
      mesh->copy(basemesh);
      space->set_uniform_order(P_INIT);
      space->assign_dofs();

      error_stop = ERR_STOP / std::pow(4.0, (double)error_level);

      double factor = std::abs(std::sin( 0.5 * M_PI * std::pow((double)(iteration + 1) / (double)iterations_count, 4.0)));

      wf.p_1 = (25.0) * factor * factor;
      wf.p_2 = (7.0) * factor;
      wf.p_3 = (5.0) * factor;
      wf.p_4 = (0.2) * factor;
      wf.p_5 = (0.05) * factor;

      wf.q_1 = 25. - (25.0) * factor;
      wf.q_2 = 0.8 - (0.8) * factor;
      wf.q_3 = 0.0001 - (0.0001) * factor;
      wf.q_4 = 0.2 - (0.2) * factor;
      wf.q_5 = 0.05 - (0.05) * factor;

      /*
      wf.c_left = factor;
      wf.c_top = factor + 1.;
      wf.c_right = factor + 2.;
      wf.c_bottom = factor + 3.;

      wf.g_n_left(0.0),
      wf.g_n_top(3.0),
      wf.g_n_right(2.0),
      wf.g_n_bottom(1.0)
      */

      newton.set_weak_formulation(&wf);

      logger.info("Iteration: %i-%i, Error level: %g, Factor: %g.", iteration, error_level, error_stop, factor);

      // Cumulative.
      int dof_cumulative = 0;
      int total_cache_searches = 0;
      int total_cache_record_found = 0;
      int total_cache_record_found_reinit = 0;
      int total_cache_record_not_found = 0;

      // Max.
      int as = 1;
      int dof_reached;
      double error_reached;
      double exact_error_reached = 0;

      // One step.
      int cache_searches;
      int cache_record_found;
      int cache_record_found_reinit;
      int cache_record_not_found;
      double FactorizationSize;
      double PeakMemoryUsage;
      double Flops;

      // Time measurement.
      Hermes::Mixins::TimeMeasurable cpu_time;

      // Tick.
      cpu_time.tick();

      try
      {
        while (!adaptive_step_single_space(
          &logger,
          mesh, 
          space, 
          sln, 
          &selector,
          is_p(CAND_LIST) ? 1 : 0,
          ref_sln, 
          cpu_time,
          newton,
          sview,
          oview,
          errorCalculator,
          adaptivity,
          as,
          error_stop,
          error_reached,
          dof_reached,
          cache_searches,
          cache_record_found,
          cache_record_found_reinit,
          cache_record_not_found,
          exact_error_reached,
          FactorizationSize,
          PeakMemoryUsage,
          Flops))
        {
          dof_cumulative += dof_reached;

          total_cache_searches += cache_searches;
          total_cache_record_found += cache_record_found;
          total_cache_record_found_reinit += cache_record_found_reinit;
          total_cache_record_not_found += cache_record_not_found;
        }
      }
      catch(std::exception& e)
      {
        data.close();
        return -1;
      }

      dof_cumulative += dof_reached;

      total_cache_searches += cache_searches;
      total_cache_record_found += cache_record_found;
      total_cache_record_found_reinit += cache_record_found_reinit;
      total_cache_record_not_found += cache_record_not_found;

      cpu_time.tick();
      {
        sprintf(filename, "Solution-%s-%s-%i-%i.vtk", get_cand_list_str(CAND_LIST), THRESHOLD_STRING, error_level, iteration);
        lin.save_solution_vtk(ref_sln, filename, "sln", false, 1, HERMES_EPS_LOW);
        sprintf(filename, "Orders-%s-%s-%i-%i.vtk", get_cand_list_str(CAND_LIST), THRESHOLD_STRING, error_level, iteration);
        ord.save_orders_vtk(newton.get_space(0), filename);
        sprintf(filename, "Mesh-%s-%s-%i-%i.vtk", get_cand_list_str(CAND_LIST), THRESHOLD_STRING, error_level, iteration);
        ord.save_mesh_vtk(newton.get_space(0), filename);
      }
      cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

      data << 
        iteration << ';' <<
        cpu_time.accumulated() << ';' <<
        as - 1 << ';' <<
        dof_reached << ';' <<
        dof_cumulative << ';' <<
        total_cache_searches << ';' <<
        total_cache_record_found << ';' <<
        total_cache_record_found_reinit << ';' <<
        total_cache_record_not_found << ';' <<
        error_stop << ';' <<
        error_reached << ';' <<
        exact_error_reached <<
        std::endl;

      std::cout << std::endl << "Results:" << std::endl;
      std::cout << "CPU time: " << cpu_time.accumulated_str() << std::endl;
      std::cout << "Adaptivity steps: " << as - 1 << std::endl;
      std::cout << "dof_reached: " << dof_reached << std::endl;
      std::cout << "dof_cumulative: " << dof_cumulative << std::endl;

      std::cout << "total_cache_searches: " << total_cache_searches << std::endl;
      std::cout << "total_cache_record_found: " << total_cache_record_found << std::endl;
      std::cout << "total_cache_record_found_reinit: " << total_cache_record_found_reinit << std::endl;
      std::cout << "total_cache_record_not_found: " << total_cache_record_not_found << std::endl;

      std::cout << "error_stop: " << error_stop << std::endl;
      std::cout << "error_reached: " << error_reached << std::endl;
      std::cout << "exact_error_reached: " << exact_error_reached << std::endl;

    }
  }

  data.close();
  return 0;
}