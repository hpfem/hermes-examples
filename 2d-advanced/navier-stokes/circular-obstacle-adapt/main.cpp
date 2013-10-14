

#include "definitions.h"

// The time-dependent laminar incompressible Navier-Stokes equations are
// discretized in time via the implicit Euler method. The Newton's method 
// is used to solve the nonlinear problem at each time step. We show how
// to use discontinuous ($L^2$) elements for pressure and thus make the
// velocity discreetely divergence free. Comparison to approximating the
// pressure with the standard (continuous) Taylor-Hood elements is enabled.
// The Reynolds number Re = 200 which is embarrassingly low. You
// can increase it but then you will need to make the mesh finer, and the
// computation will take more time.
//
// PDE: incompressible Navier-Stokes equations in the form
// \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
// div v = 0
//
// BC: u_1 is a time-dependent constant and u_2 = 0 on Gamma_4 (inlet)
//     u_1 = u_2 = 0 on Gamma_1 (bottom), Gamma_3 (top) and Gamma_5 (obstacle)
//     "do nothing" on Gamma_2 (outlet)
//
// Geometry: Rectangular channel containing an off-axis circular obstacle. The
//           radius and position of the circle, as well as other geometry
//           parameters can be changed in the mesh file "domain.mesh".
//
// The following parameters can be changed:

// For application of Stokes flow (creeping flow).
const bool STOKES = false;                        
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Number of initial mesh refinements towards boundary.
const int INIT_REF_NUM_BDY = 3;                   
// If this is defined, the pressure is approximated using
// discontinuous L2 elements (making the velocity discreetely
// divergence-free, more accurate than using a continuous
// pressure approximation). Otherwise the standard continuous
// elements are used. The results are striking - check the
// tutorial for comparisons.
//#define PRESSURE_IN_L2                            
// Initial polynomial degree for velocity components.
// Note: P_INIT_VEL should always be greater than
// P_INIT_PRESSURE because of the inf-sup condition.
const int P_INIT_VEL = 2;                         
// Initial polynomial degree for pressure.
const int P_INIT_PRESSURE = 1;                    

// Adaptivity
// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 1;                         
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;                     
// Error calculation & adaptivity.
class CustomErrorCalculator : public DefaultErrorCalculator<double, HERMES_H1_NORM>
{
public:
  // Two first components are in the H1 space - we can use the classic class for that, for the last component, we will manually add the L2 norm for pressure.
  CustomErrorCalculator(CalculatedErrorType errorType) : DefaultErrorCalculator<double, HERMES_H1_NORM>(errorType, 2)
  {
    this->add_error_form(new DefaultNormFormVol<double>(2, 2, HERMES_L2_NORM));
  }
} errorCalculator(RelativeErrorToGlobalNorm);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);      
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_H_ISO;           
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-3;

// Problem parameters
// Reynolds number.
const double RE = 200.0;                          
// Inlet velocity (reached after STARTUP_TIME).
const double VEL_INLET = 1.0;                     
// During this time, inlet velocity increases gradually
// from 0 to VEL_INLET, then it stays constant.
const double STARTUP_TIME = 1.0;                  
// Time step.
const double TAU = 0.001;                          
// Time interval length.
const double T_FINAL = 30000.0;                   
// Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL = 0.05;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 20;                   
// Domain height (necessary to define the parabolic
// velocity profile at inlet).
const double H = 5;                                    

// Current time (defined as global since needed in weak forms).
double TIME = 0;

// Boundary markers.
const std::string BDY_BOTTOM = "b1";
const std::string BDY_RIGHT = "b2";
const std::string BDY_TOP = "b3";
const std::string BDY_LEFT = "b4";
const std::string BDY_OBSTACLE = "b5";

// Current time (used in weak forms).
double current_time = 0;

/*// Boundary condition values for x-velocity
double essential_bc_values_xvel(double x, double y, double time) {
// time-dependent inlet velocity (parabolic profile)
double val_y = VEL_INLET * y*(H-y) / (H/2.)/(H/2.); //parabolic profile with peak VEL_INLET at y = H/2
if (time <= STARTUP_TIME) return val_y * time/STARTUP_TIME;
else return val_y;
}
*/

/*
void mag(int n, double* a, double* dadx, double* dady,
double* b, double* dbdx, double* dbdy,
double* out, double* outdx, double* outdy)
{
for (int i = 0; i < n; i++)
{
out[i] = sqrt(sqr(a[i]) + sqr(b[i]));
outdx[i] = (0.5 / out[i]) * (2.0 * a[i] * dadx[i] + 2.0 * b[i] * dbdx[i]);
outdy[i] = (0.5 / out[i]) * (2.0 * a[i] * dady[i] + 2.0 * b[i] * dbdy[i]);
}
}
*/

int main(int argc, char* argv[])
{
  HermesCommonApi.set_integral_param_value(numThreads, 1);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", basemesh);

  // Initial mesh refinements.
  //mesh->refine_all_elements();
  basemesh->refine_towards_boundary(BDY_OBSTACLE, 1, false);
  // '4' is the number of levels.
  basemesh->refine_towards_boundary(BDY_TOP, 1, true);     
  // 'true' stands for anisotropic refinements.
  basemesh->refine_towards_boundary(BDY_BOTTOM, 1, true);  

  mesh->copy(basemesh);

  // Initialize boundary conditions.
  EssentialBCNonConst bc_left_vel_x(BDY_LEFT, VEL_INLET, H, STARTUP_TIME);
  DefaultEssentialBCConst<double> bc_other_vel_x(Hermes::vector<std::string>(BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  EssentialBCs<double> bcs_vel_x(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_left_vel_x, &bc_other_vel_x));
  DefaultEssentialBCConst<double> bc_vel_y(Hermes::vector<std::string>(BDY_LEFT, BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  EssentialBCs<double> bcs_vel_y(&bc_vel_y);

  SpaceSharedPtr<double> xvel_space(new H1Space<double>(mesh, &bcs_vel_x, P_INIT_VEL));
  SpaceSharedPtr<double> yvel_space(new H1Space<double>(mesh, &bcs_vel_y, P_INIT_VEL));
#ifdef PRESSURE_IN_L2
  SpaceSharedPtr<double> p_space(new L2Space<double>(mesh, P_INIT_PRESSURE));
#else
  SpaceSharedPtr<double> p_space(new H1Space<double>(mesh, P_INIT_PRESSURE));
#endif
  Hermes::vector<SpaceSharedPtr<double> > spaces(xvel_space, yvel_space, p_space);

  // Calculate and report the number of degrees of freedom.
  int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof = %d.", ndof);

  // Define projection norms.
  NormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  NormType p_proj_norm = HERMES_L2_NORM;
#else
  NormType p_proj_norm = HERMES_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  Hermes::Mixins::Loggable::Static::info("Setting initial conditions.");

  // Define initial conditions on the
  MeshFunctionSharedPtr<double>  xvel_ref_sln(new Solution<double>());
  MeshFunctionSharedPtr<double>  yvel_ref_sln(new Solution<double>());
  MeshFunctionSharedPtr<double>  p_ref_sln(new Solution<double>());

  MeshFunctionSharedPtr<double>  xvel_prev_time(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double>  yvel_prev_time(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double>  p_prev_time(new ZeroSolution<double>(mesh));

  MeshFunctionSharedPtr<double>  xvel_sln(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double>  yvel_sln(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double>  p_sln(new ZeroSolution<double>(mesh));

  // Initialize weak formulation.
  WeakFormNSNewton wf(STOKES, RE, TAU, xvel_prev_time, yvel_prev_time);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, spaces);

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  // Initialize views.
  VectorView vview("velocity [m/s]", new WinGeom(0, 0, 750, 240));
  ScalarView pview("pressure [Pa]", new WinGeom(0, 290, 750, 240));
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  pview.fix_scale_width(80);
  pview.show_mesh(true);

  // Time-stepping loop:
  char title[100];
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    current_time += TAU;
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time = %g:", ts, current_time);

    // Update time-dependent essential BCs.
    Hermes::Mixins::Loggable::Static::info("Updating time-dependent essential BC.");
    Space<double>::update_essential_bc_values(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space), current_time);

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
      mesh->copy(basemesh);
      xvel_space->set_uniform_order(P_INIT_VEL);
      yvel_space->set_uniform_order(P_INIT_VEL);
      p_space->set_uniform_order(P_INIT_PRESSURE);

      xvel_space->assign_dofs();
      yvel_space->assign_dofs();
      p_space->assign_dofs();
    }

    DiscreteProblem<double> dp(&wf, spaces);
    Hermes::Hermes2D::NewtonSolver<double> newton(&dp);

    // Spatial adaptivity loop. Note: xvel_prev_time, yvel_prev_time and pvel_prev_time
    // must not be changed during spatial adaptivity. 
    bool done = false; int as = 1;
    double err_est;
    do {
      Hermes::Mixins::Loggable::Static::info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Mesh::ReferenceMeshCreator refMeshCreator(mesh);
      MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

      Space<double>::ReferenceSpaceCreator refSpaceCreatorX(xvel_space, ref_mesh, 0);
      SpaceSharedPtr<double> ref_xvel_space = refSpaceCreatorX.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorY(yvel_space, ref_mesh, 0);
      SpaceSharedPtr<double> ref_yvel_space = refSpaceCreatorY.create_ref_space();
      Space<double>::ReferenceSpaceCreator refSpaceCreatorP(p_space, ref_mesh, 0);
      SpaceSharedPtr<double> ref_p_space = refSpaceCreatorP.create_ref_space();

      Hermes::vector<SpaceSharedPtr<double> > ref_spaces(ref_xvel_space, ref_yvel_space, ref_p_space);

      // Calculate initial coefficient vector for Newton on the fine mesh.
      double* coeff_vec = new double[Space<double>::get_num_dofs(ref_spaces)];

      if (ts == 1) {
        Hermes::Mixins::Loggable::Static::info("Projecting coarse mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection<double> ogProj; ogProj.project_global(ref_spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(xvel_sln, yvel_sln, p_sln), 
          coeff_vec);
      }
      else {
        Hermes::Mixins::Loggable::Static::info("Projecting previous fine mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection<double> ogProj; ogProj.project_global(ref_spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(xvel_prev_time, yvel_prev_time, p_prev_time), 
          coeff_vec);
      }

      // Perform Newton's iteration.
      Hermes::Mixins::Loggable::Static::info("Solving nonlinear problem:");
      try
      {
        newton.set_spaces(ref_spaces);
        newton.set_max_allowed_iterations(NEWTON_MAX_ITER);
        newton.set_tolerance(NEWTON_TOL, Hermes::Solvers::ResidualNormAbsolute);
        if(as == 2)
          newton.output_matrix();
        newton.solve(coeff_vec);
      }
      catch(Hermes::Exceptions::Exception e)
      {
        e.print_msg();
      };

      // Update previous time level solutions.
      Solution<double>::vector_to_solutions(newton.get_sln_vector(), ref_spaces, Hermes::vector<MeshFunctionSharedPtr<double> >(xvel_ref_sln, yvel_ref_sln, p_ref_sln));

      // Project the fine mesh solution onto the coarse mesh.
      Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
      OGProjection<double> ogProj; ogProj.project_global(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space), 
        Hermes::vector<MeshFunctionSharedPtr<double> >(xvel_ref_sln, yvel_ref_sln, p_ref_sln), 
        Hermes::vector<MeshFunctionSharedPtr<double> >(xvel_sln, yvel_sln, p_sln), 
        Hermes::vector<NormType>(vel_proj_norm, vel_proj_norm, p_proj_norm) );

      // Calculate element errors and total error estimate.
      Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
      adaptivity.set_spaces(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space));

      errorCalculator.calculate_errors(Hermes::vector<MeshFunctionSharedPtr<double> >(xvel_sln, yvel_sln, p_sln), 
        Hermes::vector<MeshFunctionSharedPtr<double> >(xvel_ref_sln, yvel_ref_sln, p_ref_sln));
      double err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

      // Report results.
      Hermes::Mixins::Loggable::Static::info("ndof: %d, ref_ndof: %d, err_est_rel: %g%%", 
        Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space)), 
        Space<double>::get_num_dofs(ref_spaces), err_est_rel_total);

      // Show the solution at the end of time step.
      sprintf(title, "Velocity, time %g", TIME);
      vview.set_title(title);
      vview.show(xvel_ref_sln, yvel_ref_sln);
      sprintf(title, "Pressure, time %g", TIME);
      pview.set_title(title);
      pview.show(p_ref_sln);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
      else 
      {
        Hermes::Mixins::Loggable::Static::info("Adapting the coarse mesh.");
        done = adaptivity.adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector, &selector));

        // Increase the counter of performed adaptivity steps.
        as++;
      }

      // Clean up.
      delete [] coeff_vec;
    }
    while (done == false);

    // Copy new time level reference solution into prev_time.
    xvel_prev_time->copy(xvel_ref_sln);
    yvel_prev_time->copy(yvel_ref_sln);
    p_prev_time->copy(p_ref_sln);
  }

  ndof = Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space));
  Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
