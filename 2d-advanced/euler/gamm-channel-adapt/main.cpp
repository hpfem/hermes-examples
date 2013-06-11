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
double CFL_NUMBER = 0.3;                          
// Initial time step.
double time_step_n = 1E-6;                        

// Adaptivity.
// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 10;

// Number of mesh refinements between two unrefinements.
// The mesh is not unrefined unless there has been a refinement since
// last unrefinement.
int REFINEMENT_COUNT = 1;                         

// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies (see below).
const double THRESHOLD = 0.6;                     

// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
CandList CAND_LIST = H2D_HP_ANISO;                

// Stopping criterion for adaptivity.
const double ERR_STOP = 0.007;                     

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
  #pragma region 1. Load mesh and initialize spaces.
    // Load the mesh.
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
  #pragma endregion

  #pragma region 2. Initialize solutions, set initial conditions.
    MeshFunctionSharedPtr<double> sln_rho(new ConstantSolution<double>(mesh, RHO_EXT));
    MeshFunctionSharedPtr<double> sln_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
    MeshFunctionSharedPtr<double> sln_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
    MeshFunctionSharedPtr<double> sln_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
    Hermes::vector<MeshFunctionSharedPtr<double> > slns(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e);

    MeshFunctionSharedPtr<double> prev_rho(new ConstantSolution<double>(mesh, RHO_EXT));
    MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
    MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
    MeshFunctionSharedPtr<double> prev_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
    Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);

    MeshFunctionSharedPtr<double> rsln_rho(new ConstantSolution<double>(mesh, RHO_EXT));
    MeshFunctionSharedPtr<double> rsln_rho_v_x(new ConstantSolution<double> (mesh, RHO_EXT * V1_EXT));
    MeshFunctionSharedPtr<double> rsln_rho_v_y(new ConstantSolution<double> (mesh, RHO_EXT * V2_EXT));
    MeshFunctionSharedPtr<double> rsln_e(new ConstantSolution<double> (mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));
    Hermes::vector<MeshFunctionSharedPtr<double> > rslns(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e);
  #pragma endregion

  #pragma region 3. Filters for visualization of Mach number, pressure + visualization setup.
    MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(rslns, KAPPA));
    MeshFunctionSharedPtr<double>  pressure(new PressureFilter(rslns, KAPPA));

    ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
    pressure_view.show_contours(.1);
    pressure_view.show_mesh(false);
    ScalarView Mach_number_view("Mach number", new WinGeom(650, 0, 600, 300));
    Mach_number_view.show_contours(.02);
    Mach_number_view.show_mesh(false);
    ScalarView eview("Error - density", new WinGeom(0, 330, 600, 300));
    ScalarView eview1("Error - momentum", new WinGeom(0, 660, 600, 300));
    OrderView order_view("Orders", new WinGeom(650, 330, 600, 300));
  #pragma endregion

  #pragma region 4. Adaptivity setup.
    // Initialize refinement selector.
    L2ProjBasedSelector<double> selector(CAND_LIST);
    selector.set_dof_score_exponent(2.0);

    //selector.set_error_weights(1.0, 1.0, 1.0);

    // Error calculation.
    DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 4);
    // Stopping criterion for an adaptivity step.
    AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
    Adapt<double> adaptivity(spaces, &errorCalculator, &stoppingCriterion);
  #pragma endregion

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  #pragma region 5. Initialize weak formulation and solver.
    Hermes::vector<std::string> solid_wall_markers(BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP);
    Hermes::vector<std::string> inlet_markers;
    inlet_markers.push_back(BDY_INLET);
    Hermes::vector<std::string> outlet_markers;
    outlet_markers.push_back(BDY_OUTLET);
    EulerEquationsWeakFormSemiImplicit wf(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, solid_wall_markers, inlet_markers, outlet_markers, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
  
    // Solver.
    LinearSolver<double> solver(&wf, spaces);
  #pragma endregion
  
  #pragma region 6. Time stepping loop.
    int iteration = 1;
    for(double t = 0.0; t < 50.; t += time_step_n)
    {
      Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

      #pragma region 6.1. Periodic global derefinements.
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
      #pragma endregion

      #pragma region 7. Adaptivity loop.
        int as = 1; int ndofs_prev = 0; bool done = false;
        do
        {
          // Info.
          Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);
          // Set the current time step.
          wf.set_current_time_step(time_step_n);

          #pragma region 7.1. Construct globally refined reference mesh and setup reference space.
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
            solver.set_spaces(ref_spaces);

            if(ndofs_prev != 0)
              if(Space<double>::get_num_dofs(ref_spaces) == ndofs_prev)
                selector.set_error_weights(2.0 * selector.get_error_weight_h(), 1.0, 1.0);
              else
                selector.set_error_weights(1.0, 1.0, 1.0);

            ndofs_prev = Space<double>::get_num_dofs(ref_spaces);

            // Project the previous time level solution onto the new fine mesh
            Hermes::Mixins::Loggable::Static::info("Projecting the previous time level solution onto the new fine mesh.");
            OGProjection<double>::project_global(ref_spaces, prev_slns, prev_slns);
          #pragma endregion

          // Solve the problem.
          solver.solve();

          #pragma region *. Get the solution with optional shock capturing.
            if(!SHOCK_CAPTURING)
              Solution<double>::vector_to_solutions(solver.get_sln_vector(), ref_spaces, rslns);
            else
            {
              FluxLimiter* flux_limiter;
              if(SHOCK_CAPTURING_TYPE == KUZMIN)
                flux_limiter = new FluxLimiter(FluxLimiter::Kuzmin, solver.get_sln_vector(), ref_spaces);
              else
                flux_limiter = new FluxLimiter(FluxLimiter::Krivodonova, solver.get_sln_vector(), ref_spaces);
          
              if(SHOCK_CAPTURING_TYPE == KUZMIN)
                flux_limiter->limit_second_orders_according_to_detector(spaces);

              flux_limiter->limit_according_to_detector(spaces);

              flux_limiter->get_limited_solutions(rslns);

              delete flux_limiter;
            }
          #pragma endregion

          // Calculate time step according to CFL condition.
          CFL.calculate(rslns, (ref_spaces)[0]->get_mesh(), time_step_n);

          #pragma region 7.2. Project to coarse mesh -> error estimation -> space adaptivity
            // Project the fine mesh solution onto the coarse mesh.
            Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
            OGProjection<double>::project_global(spaces, rslns, slns, Hermes::vector<NormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

            // Calculate element errors and total error estimate.
            Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
            errorCalculator.calculate_errors(slns, rslns);
            double err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

            // Report results.
            Hermes::Mixins::Loggable::Static::info("err_est_rel: %g%%", err_est_rel_total);

            // If err_est too large, adapt the mesh.
            if (err_est_rel_total < ERR_STOP)
              done = true;
            else
            {
              Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
              done = adaptivity.adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector, &selector, &selector));

              REFINEMENT_COUNT++;
              if (Space<double>::get_num_dofs(spaces) >= NDOF_STOP) 
                done = true;
              else
                as++;
            }
          #pragma endregion
              
          #pragma region 7.3. Visualization and saving on disk.
            if(done && (iteration - 1) % EVERY_NTH_STEP == 0 && iteration > 2)
            {
              // Hermes visualization.
              if(HERMES_VISUALIZATION)
              {        
                Mach_number->reinit();
                pressure->reinit();
                pressure_view.show(pressure, 1);
                Mach_number_view.show(Mach_number, 1);
                order_view.show((ref_spaces)[0]);
              }
            }
          #pragma endregion
        }
        while (done == false);
      #pragma endregion

      // Copy the solutions into the previous time level ones.
      prev_rho->copy(rsln_rho);
      prev_rho_v_x->copy(rsln_rho_v_x);
      prev_rho_v_y->copy(rsln_rho_v_y);
      prev_e->copy(rsln_e);
    }
  #pragma endregion

  pressure_view.close();
  Mach_number_view.close();

  return 0;
}
