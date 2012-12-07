#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;
using namespace RefinementSelectors;

//  This is the IAEA EIR-2 benchmark problem. Note the way of handling different material
//  parameters. This is an alternative to how this is done in tutorial examples 07 and 12
//  and in example "iron-water".
//
//  PDE: -div(D(x,y)grad\Phi) + \Sigma_a(x,y)\Phi = Q_{ext}(x,y)
//  where D(x, y) is the diffusion coefficient, \Sigma_a(x,y) the absorption cross-section,
//  and Q_{ext}(x,y) external sources.
//
//  Domain: square (0, L)x(0, L) where L = 30c (see mesh file domain.mesh).
//
//  BC:  Zero Dirichlet for the right and top edges ("vacuum boundary").
//       Zero Neumann for the left and bottom edges ("reflection boundary").
//
//  The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;                       
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.6;                     
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
const int STRATEGY = 0;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO_H;        
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
// Stopping criterion for adaptivity (rel. error tolerance between the
// reference mesh and coarse mesh solution in percent).
const double ERR_STOP = 0.1;                      
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 100000;                     
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Problem parameters
// Total horizontal length.
double LH = 96;                                   
// First horizontal length.
double LH0 = 18;                                  
// Second horizontal length.
double LH1 = 48;                                  
// Third horizontal length.
double LH2 = 78;                                  
// Total vertical length.
double LV = 96;                                   
// First vertical length.
double LV0 = 18;                                  
// Second vertical length.
double LV1 = 48;                                  
// Third vertical length.
double LV2 = 78;                                  
// Total cross-sections.
double SIGMA_T_1 = 0.60;                          
double SIGMA_T_2 = 0.48;
double SIGMA_T_3 = 0.70;
double SIGMA_T_4 = 0.85;
double SIGMA_T_5 = 0.90;
// Scattering cross sections.
double SIGMA_S_1 = 0.53;                          
double SIGMA_S_2 = 0.20;
double SIGMA_S_3 = 0.66;
double SIGMA_S_4 = 0.50;
double SIGMA_S_5 = 0.89;
// Nonzero sources in regions 1 and 3 only.
// Sources in other regions are zero.
double Q_EXT_1 = 1;                               
double Q_EXT_3 = 1;                               

// Additional constants
// Diffusion coefficients.
double D_1 = 1/(3.*SIGMA_T_1);                    
double D_2 = 1/(3.*SIGMA_T_2);
double D_3 = 1/(3.*SIGMA_T_3);
double D_4 = 1/(3.*SIGMA_T_4);
double D_5 = 1/(3.*SIGMA_T_5);
// Absorption coefficients.
double SIGMA_A_1 = SIGMA_T_1 - SIGMA_S_1;         
double SIGMA_A_2 = SIGMA_T_2 - SIGMA_S_2;
double SIGMA_A_3 = SIGMA_T_3 - SIGMA_S_3;
double SIGMA_A_4 = SIGMA_T_4 - SIGMA_S_4;
double SIGMA_A_5 = SIGMA_T_5 - SIGMA_S_5;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);
  
  // Perform initial uniform mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Set essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("right", "top"), 0.0);
  EssentialBCs<double> bcs(&bc_essential);
  
  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);

  // Associate element markers (corresponding to physical regions) 
  // with material properties (diffusion coefficient, absorption 
  // cross-section, external sources).
  Hermes::vector<std::string> regions("e1", "e2", "e3", "e4", "e5");
  Hermes::vector<double> D_map(D_1, D_2, D_3, D_4, D_5);
  Hermes::vector<double> Sigma_a_map(SIGMA_A_1, SIGMA_A_2, SIGMA_A_3, SIGMA_A_4, SIGMA_A_5);
  Hermes::vector<double> Sources_map(Q_EXT_1, 0.0, Q_EXT_3, 0.0, 0.0);
  
  // Initialize the weak formulation.
  WeakFormsNeutronics::Monoenergetic::Diffusion::DefaultWeakFormFixedSource<double>
    wf(regions, D_map, Sigma_a_map, Sources_map);

  // Initialize coarse and reference mesh solution.
  Solution<double> sln, ref_sln;
  
  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  OrderView  oview("Polynomial orders", new WinGeom(450, 0, 400, 350));
  
  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;
  
  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);
    
    // Time measurement.
    cpu_time.tick();

    // Construct globally refined mesh and setup fine mesh space.
    Mesh::ReferenceMeshCreator refMeshCreator(&mesh);
    Mesh* ref_mesh = refMeshCreator.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreator(&space, ref_mesh);
    Space<double>* ref_space = refSpaceCreator.create_ref_space();
    int ndof_ref = ref_space->get_num_dofs();

    // Initialize fine mesh problem.
    Hermes::Mixins::Loggable::Static::info("Solving on fine mesh.");
    DiscreteProblem<double> dp(&wf, ref_space);
    
    NewtonSolver<double> newton(&dp);
    //newton.set_verbose_output(false);

    // Perform Newton's iteration.
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    }

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);
    
    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting fine mesh solution on coarse mesh.");
    OGProjection<double> ogProjection; ogProjection.project_global(&space, &ref_sln, &sln);

    // Time measurement.
    cpu_time.tick();

    // Visualize the solution and mesh.
    sview.show(&sln);
    oview.show(&space);

    // Skip visualization time.
    cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
    Adapt<double> adaptivity(&space);
    bool solutions_for_adapt = true;
    double err_est_rel = adaptivity.calc_err_est(&sln, &ref_sln, solutions_for_adapt,
                         HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      space.get_num_dofs(), ref_space->get_num_dofs(), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    cpu_time.tick();    
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    graph_dof.add_values(space.get_num_dofs(), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    
    // Skip the time spent to save the convergence graphs.
    cpu_time.tick(Hermes::Mixins::TimeMeasurable::HERMES_SKIP);

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) 
      done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of performed adaptivity steps.
      if (done == false)  
        as++;
    }
    if (space.get_num_dofs() >= NDOF_STOP) 
      done = true;

    // Keep the mesh from final step to allow further work with the final fine mesh solution.
    if(done == false)
    {
      delete ref_space->get_mesh(); 
      delete ref_space;
    }
  }
  while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Show the fine mesh solution - final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(&ref_sln);

  // Wait for all views to be closed.
  Views::View::wait();

  return 0;
}

