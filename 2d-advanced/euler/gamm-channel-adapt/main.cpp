#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// This example solves the compressible Euler equations using 
// Discontinuous Galerkin method of higher order with adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: GAMM channel, see mesh file GAMM-channel.mesh
//
// BC: Solid walls, inlet, no outlet.
//
// IC: Constant state identical to inlet.
//
// The following parameters can be changed:

// Visualization.
// Set to "true" to enable Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;
// Set visual output for every nth step.
const unsigned int EVERY_NTH_STEP = 1;
// Shock capturing.
enum shockCapturingType
{
  KUZMIN,
  KRIVODONOVA
};
bool SHOCK_CAPTURING = false;
shockCapturingType SHOCK_CAPTURING_TYPE = KUZMIN;
// Quantitative parameter of the discontinuity detector in case of Krivodonova.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;
// Quantitative parameter of the shock capturing in case of Feistauer.
const double NU_1 = 0.1;
const double NU_2 = 0.1;

// Initial polynomial degree. 
const int P_INIT = 0;                                                  
// Number of initial uniform mesh refinements.  
const int INIT_REF_NUM = 2;
// CFL value.
double CFL_NUMBER = 0.5;                          
// Initial time step.
double time_step_n = 1E-6;                        

// Adaptivity.
// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 10;

// Number of mesh refinements between two unrefinements.
// The mesh is not unrefined unless there has been a refinement since
// last unrefinement.
int REFINEMENT_COUNT = 0;                         

// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies (see below).
const double THRESHOLD = 0.3;                     

// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
const int STRATEGY = 1;

// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
CandList CAND_LIST = H2D_HP_ANISO;                

// Maximum polynomial degree used. -1 for unlimited.
const int MAX_P_ORDER = -1;                       

// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   

// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0. 
const double CONV_EXP = 1;                        

// Stopping criterion for adaptivity.
const double ERR_STOP = 0.5;                     

// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 100000;                   

// Matrix solver for orthogonal projections: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 2.5;
// Inlet density (dimensionless).   
const double RHO_EXT = 1.0;                       
// Inlet x-velocity (dimensionless).
const double V1_EXT = 1.25;                       
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;                        
// Kappa.
const double KAPPA = 1.4;                         

// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh->
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("GAMM-channel.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) 
    mesh->refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT));
  SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT));
  SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT));
  SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT));
  SpaceSharedPtr<double> space_stabilization(new L2Space<double>(mesh, 0));
  Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);
  int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  MeshFunctionSharedPtr<double> sln_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> sln_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> sln_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> sln_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));

  MeshFunctionSharedPtr<double> prev_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> prev_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));

  MeshFunctionSharedPtr<double> rsln_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> rsln_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> rsln_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> rsln_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));

  // Filters for visualization of Mach number, pressure and entropy.
  MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), KAPPA));
  MeshFunctionSharedPtr<double>  pressure(new PressureFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), KAPPA));
  MeshFunctionSharedPtr<double>  entropy(new EntropyFilter(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), KAPPA, RHO_EXT, P_EXT));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  ScalarView s1("prev_rho", new WinGeom(0, 0, 600, 300));
  ScalarView s2("prev_rho_v_x", new WinGeom(700, 0, 600, 300));
  ScalarView s3("prev_rho_v_y", new WinGeom(0, 400, 600, 300));
  ScalarView s4("prev_e", new WinGeom(700, 400, 600, 300));

  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, MAX_P_ORDER);
  selector.set_error_weights(1.0, 1.0, 1.0);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Time stepping loop.
  int iteration = 1;
  double t = 0.0;
  for(; t < 3.5; t += time_step_n)
  {
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0) 
    {
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      REFINEMENT_COUNT = 0;

      space_rho->unrefine_all_mesh_elements(true);

      space_rho->adjust_element_order(-1, P_INIT);
      space_rho_v_x->adjust_element_order(-1, P_INIT);
      space_rho_v_y->adjust_element_order(-1, P_INIT);
      space_e->adjust_element_order(-1, P_INIT);
      Space<double>::assign_dofs(spaces);
    }

    // Adaptivity loop:
    int as = 1; 
    int ndofs_prev = 0;
    bool done = false;
    do
    {
      Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      int order_increase = CAND_LIST == H2D_HP_ANISO ? 1 : 0;

      Mesh::ReferenceMeshCreator refMeshCreatorFlow(mesh);
      MeshSharedPtr ref_mesh = refMeshCreatorFlow.create_ref_mesh();

      Space<double>::ReferenceSpaceCreator refSpaceCreatorRho(space_rho, ref_mesh, order_increase);
      SpaceSharedPtr<double> ref_space_rho = refSpaceCreatorRho.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVx(space_rho_v_x, ref_mesh, order_increase);
      SpaceSharedPtr<double> ref_space_rho_v_x = refSpaceCreatorRhoVx.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVy(space_rho_v_y, ref_mesh, order_increase);
      SpaceSharedPtr<double> ref_space_rho_v_y = refSpaceCreatorRhoVy.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorE(space_e, ref_mesh, order_increase);
      SpaceSharedPtr<double> ref_space_e = refSpaceCreatorE.create_ref_space();

      Hermes::vector<SpaceSharedPtr<double>  > ref_spaces(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e);

      SpaceSharedPtr<double> refspace_stabilization(new L2Space<double>(ref_space_rho->get_mesh(), 0));

      if(ndofs_prev != 0)
        if(Space<double>::get_num_dofs(ref_spaces) == ndofs_prev)
          selector.set_error_weights(2.0 * selector.get_error_weight_h(), 1.0, 1.0);
        else
          selector.set_error_weights(1.0, 1.0, 1.0);

      ndofs_prev = Space<double>::get_num_dofs(ref_spaces);

      // Project the previous time level solution onto the new fine mesh->
      Hermes::Mixins::Loggable::Static::info("Projecting the previous time level solution onto the new fine mesh.");
      OGProjection<double> ogProjection; ogProjection.project_global(ref_spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), 
          Hermes::vector<MeshFunctionSharedPtr<double> >(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e), Hermes::vector<Hermes::Hermes2D::ProjNormType>(), iteration > 1);
      
      // Initialize weak formulation.
      Hermes::vector<std::string> solid_wall_markers(BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP);
      Hermes::vector<std::string> inlet_markers;
      inlet_markers.push_back(BDY_INLET);
      Hermes::vector<std::string> outlet_markers;
      outlet_markers.push_back(BDY_OUTLET);

      EulerEquationsWeakFormSemiImplicit wf(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT,solid_wall_markers, inlet_markers, outlet_markers, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);

      EulerEquationsWeakFormStabilization wf_stabilization(prev_rho);

      
      // Assemble the reference problem.
      Hermes::Mixins::Loggable::Static::info("Solving on reference mesh.");
      Space<double>::assign_dofs(ref_spaces);
      DiscreteProblem<double> dp(&wf, ref_spaces);
      dp.set_linear();
      DiscreteProblem<double> dp_stabilization(&wf_stabilization, refspace_stabilization);
      bool* discreteIndicator = NULL;

      SparseMatrix<double>* matrix = create_matrix<double>();
      Vector<double>* rhs = create_vector<double>();
      Vector<double>* rhs_stabilization = create_vector<double>();
      LinearMatrixSolver<double>* solver = create_linear_solver<double>( matrix, rhs);

      // Set the current time step.
      wf.set_current_time_step(time_step_n);

      FluxLimiter* flux_limiter;
      dp.assemble(matrix, rhs);

      // Solve the matrix problem.
      Hermes::Mixins::Loggable::Static::info("Solving the matrix problem.");
      if(solver->solve())
      {
        if(!SHOCK_CAPTURING)
        {
          Solution<double>::vector_to_solutions(solver->get_sln_vector(), ref_spaces, 
            Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e));
        }
        else
        {
          if(SHOCK_CAPTURING_TYPE == KUZMIN)
            flux_limiter = new FluxLimiter(FluxLimiter::Kuzmin, solver->get_sln_vector(), ref_spaces);
          else
            flux_limiter = new FluxLimiter(FluxLimiter::Krivodonova, solver->get_sln_vector(), ref_spaces);
          if(SHOCK_CAPTURING_TYPE == KUZMIN)
            flux_limiter->limit_second_orders_according_to_detector(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
            space_rho_v_y, space_e));

          flux_limiter->limit_according_to_detector(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
            space_rho_v_y, space_e));

          flux_limiter->get_limited_solutions(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e));
        }
      }
      else
        throw Hermes::Exceptions::Exception("Matrix solver failed.\n");

      // Project the fine mesh solution onto the coarse mesh.
      Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
      ogProjection.project_global(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
        space_rho_v_y, space_e), Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), 
        Hermes::vector<MeshFunctionSharedPtr<double> >(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e), 
        Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
      Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
        space_rho_v_y, space_e), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<MeshFunctionSharedPtr<double> >(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e),
        Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e)) * 100;

      CFL.calculate(Hermes::vector<MeshFunctionSharedPtr<double> >(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e), (ref_spaces)[0]->get_mesh(), time_step_n);

      // Report results.
      Hermes::Mixins::Loggable::Static::info("err_est_rel: %g%%", err_est_rel_total);

      // If err_est too large, adapt the mesh->
      if (err_est_rel_total < ERR_STOP)
        done = true;
      else
      {
        Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector, &selector, &selector), 
          THRESHOLD, STRATEGY, MESH_REGULARITY);

        REFINEMENT_COUNT++;
        if (Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(space_rho, space_rho_v_x, 
          space_rho_v_y, space_e)) >= NDOF_STOP) 
          done = true;
        else
          as++;
      }

      // Visualization and saving on disk.
      if(done && (iteration - 1) % EVERY_NTH_STEP == 0 && iteration > 1)
      {
        // Hermes visualization.
        if(HERMES_VISUALIZATION)
        {        
          Mach_number->reinit();
          pressure->reinit();
          entropy->reinit();
          pressure_view.show(pressure, 1);
          entropy_production_view.show(entropy, 1);
          Mach_number_view.show(Mach_number, 1);
        }
        // Output solution in VTK format.
        if(VTK_VISUALIZATION)
        {
          pressure->reinit();
          Mach_number->reinit();
          entropy->reinit();
          Linearizer lin;
          char filename[40];
          sprintf(filename, "Pressure-%i.vtk", iteration - 1);
          lin.save_solution_vtk(pressure, filename, "Pressure", false);
          sprintf(filename, "Mach number-%i.vtk", iteration - 1);
          lin.save_solution_vtk(Mach_number, filename, "MachNumber", false);
          sprintf(filename, "Entropy-%i.vtk", iteration - 1);
          lin.save_solution_vtk(entropy, filename, "Entropy", false);
        }
      }

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete rhs_stabilization;
      delete adaptivity;
    }
    while (done == false);

    // Copy the solutions into the previous time level ones.
    prev_rho->copy(rsln_rho);
    prev_rho_v_x->copy(rsln_rho_v_x);
    prev_rho_v_y->copy(rsln_rho_v_y);
    prev_e->copy(rsln_e);
  }

  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();

  return 0;
}
