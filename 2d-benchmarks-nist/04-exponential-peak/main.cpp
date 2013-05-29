
#include "definitions.h"

using namespace RefinementSelectors;

//  This is the fourth in the series of NIST benchmarks with known exact solutions.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: -Laplace u + f = 0.
//
//  Known exact solution: exp(-alpha * (pow(x - x_loc, 2) + pow(y - y_loc, 2))).
//  See functions CustomExactSolution::value and CustomExactSolution::derivatives in "exact_solution.cpp".
//
//  Domain: unit square (0, 1)x(0, 1), see the file "square_tri" or "square_quad.mesh".
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

// This problem has and exponential peak in the interior of the domain.
// (x_loc, y_loc) is the location of the peak, and alpha determines the strenghth of the peak. 
double alpha = 10;        
double x_loc = 0.5;         
double y_loc = 0.5;

// Initial polynomial degree of mesh elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;

// Error stop type.
const CalculatedErrorType errorType = RelativeErrorToGlobalNorm;
// Maximum allowed level of hanging node  s.
const int MESH_REGULARITY = -1;

// Error stop value (in percent).
double ERR_STOP = 0.1;

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
  // Quadrilaterals.
  mloader.load("square_quad.mesh", mesh);
  // Triangles.
  // mloader.load("square_tri.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i<INIT_REF_NUM; i++)
    mesh->refine_all_elements();
  basemesh->copy(mesh);

  // Set exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh, alpha, x_loc, y_loc));

  // Define right-hand side.
  CustomRightHandSide f(alpha, x_loc, y_loc);

  // Initialize weak formulation.
  WeakFormsH1::DefaultWeakFormPoissonLinear<double> wf(HERMES_ANY, &f);

  // Initialize boundary conditions
  DefaultEssentialBCNonConst<double> bc_essential("Bdy", exact_sln);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Initialize approximate solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>());
  MeshFunctionSharedPtr<double> ref_sln(new Solution<double>());

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 440, 350));
  Views::OrderView oview("Polynomial orders", new Views::WinGeom(450, 0, 420, 350));
  Linearizer lin;
  Orderizer ord;
  char* filename = new char[1000];

  // Assemble the discrete problem.    
  LinearSolver<double> linear_solver;
  linear_solver.set_weak_formulation(&wf);
  linear_solver.set_UMFPACK_output(true, false);

  // Adaptivity loop.
  DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(errorType, 1);
  Adapt<double> adaptivity(space, &errorCalculator);
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
      mesh->copy(basemesh);
      space->set_uniform_order(P_INIT);
      space->assign_dofs();

      double factor = std::abs(std::sin( 0.5 * M_PI * std::pow((double)(iteration + 1) / (double)iterations_count, 4.0)));
      alpha = 10. + factor * 5000.;
      f.alpha = alpha;
      ((CustomExactSolution*)exact_sln.get())->alpha = alpha;
      
      error_stop = ERR_STOP / std::pow(4.0, (double)error_level);
      
      logger.info("Iteration: %i-%i, Error level: %g, Factor: %g, Alpha: %g%.", iteration, error_level, error_stop, factor, alpha);

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
          linear_solver,
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
        /*
        sprintf(filename, "Solution-%s-%i-%i.vtk", resultStringIdentification, error_level, iteration);
        lin.save_solution_vtk(ref_sln, filename, "sln", false, 1, HERMES_EPS_LOW);
        sprintf(filename, "Orders-%s-%i-%i.vtk", resultStringIdentification, error_level, iteration);
        ord.save_orders_vtk(linear_solver.get_space(0), filename);
        sprintf(filename, "Mesh-%s-%i-%i.vtk", resultStringIdentification, error_level, iteration);
        ord.save_mesh_vtk(linear_solver.get_space(0), filename);
        */
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