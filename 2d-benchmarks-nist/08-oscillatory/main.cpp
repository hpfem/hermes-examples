#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace RefinementSelectors;

//  This is the eight in the series of NIST benchmarks with known exact solutions. It solves
//  the Helmholtz equation and has an oscillatory solution.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  The problem is made harder for adaptive algorithms by decreasing the parameter ALPHA.
//
//  PDE: -Laplace u - u/(ALPHA + r(x, y)) + f = 0 where r(x, y) = sqrt(x*x + y*y)
//
//  Known exact solution, see functions fn() and fndd().
//
//  Domain: unit square (0, 1) x (0, 1), see the file square.mesh->
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

double alpha = 2. / M_PI;

// Initial polynomial degree of mesh elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;

// Error stop type.
const CalculatedErrorType errorType = RelativeErrorToGlobalNorm;
// Maximum allowed level of hanging node  s.
const int MESH_REGULARITY = -1;

// Error stop value (in percent).
double ERR_STOP = 1.0;

int main(int argc, char* argv[])
{
  Selector<double>* refinement_selector;
  AdaptivityStoppingCriterion<double>* stoppingCriterion;
  char* resultStringIdentification;
  if(argc > 2)
    resultStringIdentification = process_arguments_main_comparison(argc, argv, refinement_selector, stoppingCriterion);
  else
  {
    refinement_selector = new MySelector(hXORpSelectionBasedOnError);
    stoppingCriterion = new AdaptStoppingCriterionSingleElement<double>(0.5);
    resultStringIdentification = "Custom";
  }

  sprintf(Hermes::Mixins::Loggable::logFileName, "Logfile-%s.log", resultStringIdentification);
  
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square_quad.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();
  basemesh->copy(mesh);

  // Set exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh, alpha));

  // Define right-hand side.
  CustomRightHandSide f(alpha);

  // Initialize weak formulation.
  CustomWeakForm wf(&f);

  // Initialize boundary conditions
  DefaultEssentialBCNonConst<double> bc_essential("Bdy", exact_sln);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Initialize approximate solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>), ref_sln(new Solution<double>);
  
  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 440, 350));
  Views::OrderView  oview("Polynomial orders", new Views::WinGeom(450, 0, 420, 350));
  Linearizer lin;
  Orderizer ord;
  char* filename = new char[1000];

  // Adaptivity loop.
  DefaultErrorCalculator<double, HERMES_H1_NORM> error_calculator(errorType, 1);
  Adapt<double> adaptivity(space, &error_calculator);
  adaptivity.set_strategy(stoppingCriterion);

  sprintf(filename, "%s.csv", resultStringIdentification);
  std::ofstream data(filename);
  data.precision(10);
  data.setf( std::ios::fixed, std::ios::floatfield );
  data << 
      "Iteration" << ';' <<
      "ErrorLevel" << ';' <<
      "CPUTime" << ';' <<
      "AdaptivitySteps" << ';' <<
      "dof_reached" << ';' <<
      "dof_cumulative" << ';' <<
      "total_cache_searches" << ';' <<
      "total_cache_record_found" << ';' <<
      "total_cache_record_found_reinit" << ';' <<
      "total_cache_record_not_found" << ';' <<
      "max_FactorizationSize" << ';' <<
      "total_PeakMemoryUsage" << ';' <<
      "total_Flops" << ';' <<
      "error_stop" << ';' <<
      "error_reached" << ';' <<
      "exact_error_reached" <<
      std::endl;

  Hermes::Mixins::Loggable logger;
  logger.set_verbose_output(true);

  int iterations_count = 10;
  int error_levels_count = 5;
  double error_stop = ERR_STOP;

  for(int iteration = 0; iteration < iterations_count; iteration++)
  {
    for(int error_level = 0; error_level < error_levels_count; error_level++)
    {
      // Assemble the discrete problem.    
      NewtonSolver<double> newton;
      newton.set_weak_formulation(&wf);
      newton.set_UMFPACK_output(true, false);

      mesh->copy(basemesh);
      space->set_uniform_order(P_INIT);
      space->assign_dofs();

      alpha = (1. / M_PI) * std::pow(0.75, iteration);
      f.alpha = alpha;
      ((CustomExactSolution*)exact_sln.get())->alpha = alpha;

      error_stop = ERR_STOP / std::pow(2.0, (double)error_level);
      
      logger.info("Iteration: %i-%i, Error level: %g, Alpha: %g%.", iteration, error_level, error_stop, alpha);

      // Cumulative.
      int dof_cumulative = 0;
      int total_cache_searches = 0;
      int total_cache_record_found = 0;
      int total_cache_record_found_reinit = 0;
      int total_cache_record_not_found = 0;
      double total_PeakMemoryUsage = 0;
      double total_Flops = 0;

      // Max.
      int as = 1;
      int dof_reached;
      double error_reached;
      double exact_error_reached;
      double max_FactorizationSize = 0;

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
          refinement_selector,
          (argc > 2 && atoi(argv[1]) == 0) ? 0 : 1,
          ref_sln, 
          cpu_time,
          newton,
          sview,
          oview,
          error_calculator,
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
          Flops,
          exact_sln))
        {
          dof_cumulative += dof_reached;

          total_cache_searches += cache_searches;
          total_cache_record_found += cache_record_found;
          total_cache_record_found_reinit += cache_record_found_reinit;
          total_cache_record_not_found += cache_record_not_found;

          max_FactorizationSize = std::max(max_FactorizationSize, FactorizationSize);
          total_PeakMemoryUsage += PeakMemoryUsage;
          total_Flops += Flops;
        }
      }
      catch(std::exception& e)
      {
        logger.info(e.what());
        data << 
        iteration << ';' <<
        error_level << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. << ';' <<
        -1. <<
        std::endl;
        continue;
      }

      dof_cumulative += dof_reached;

      total_cache_searches += cache_searches;
      total_cache_record_found += cache_record_found;
      total_cache_record_found_reinit += cache_record_found_reinit;
      total_cache_record_not_found += cache_record_not_found;

      max_FactorizationSize = std::max(max_FactorizationSize, FactorizationSize);
      total_PeakMemoryUsage += PeakMemoryUsage;
      total_Flops += Flops;

      cpu_time.tick();
      {
        sprintf(filename, "Solution-%s-%i-%i.vtk", resultStringIdentification, error_level, iteration);
        lin.save_solution_vtk(ref_sln, filename, "sln", false, 1, HERMES_EPS_LOW);
        sprintf(filename, "Orders-%s-%i-%i.vtk", resultStringIdentification, error_level, iteration);
        ord.save_orders_vtk(newton.get_space(0), filename);
        sprintf(filename, "Mesh-%s-%i-%i.vtk", resultStringIdentification, error_level, iteration);
        ord.save_mesh_vtk(newton.get_space(0), filename);
      }
      cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

      data << 
        iteration << ';' <<
        error_level << ';' <<
        cpu_time.accumulated() << ';' <<
        as - 1 << ';' <<
        dof_reached << ';' <<
        dof_cumulative << ';' <<
        total_cache_searches << ';' <<
        total_cache_record_found << ';' <<
        total_cache_record_found_reinit << ';' <<
        total_cache_record_not_found << ';' <<
        max_FactorizationSize << ';' <<
        total_PeakMemoryUsage << ';' <<
        total_Flops << ';' <<
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

      std::cout << "max_FactorizationSize: " << max_FactorizationSize << std::endl;
      std::cout << "total_PeakMemoryUsage: " << total_PeakMemoryUsage << std::endl;
      std::cout << "total_Flops: " << total_Flops << std::endl;


      std::cout << "error_stop: " << error_stop << std::endl;
      std::cout << "error_reached: " << error_reached << std::endl;
      std::cout << "exact_error_reached: " << exact_error_reached << std::endl;

    }
  }

  data.close();
  return 0;
}

