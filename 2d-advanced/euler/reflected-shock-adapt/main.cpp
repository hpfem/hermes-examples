#define HERMES_REPORT_INFO
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// This example solves the compressible Euler equations using Discontinuous Galerkin method of higher order with adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: A rectangular channel, see channel.mesh->
//
// BC: Solid walls, inlet / outlet. In this case there are two inlets (the left, and the upper wall).
//
// IC: Constant state.
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
  FEISTAUER,
  KUZMIN,
  KRIVODONOVA
};
bool SHOCK_CAPTURING = true;
shockCapturingType SHOCK_CAPTURING_TYPE = KUZMIN;
// Quantitative parameter of the discontinuity detector in case of Krivodonova.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;
// Quantitative parameter of the shock capturing in case of Feistauer.
const double NU_1 = 0.1;
const double NU_2 = 0.1;

// Initial polynomial degree.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;
// CFL value.
double CFL_NUMBER = 0.3;
// Initial time step.
double time_step_n = 1E-6;
double TIME_INTERVAL_LENGTH = 20.;

// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 10;

// Number of mesh refinements between two unrefinements.
// The mesh is not unrefined unless there has been a refinement since
// last unrefinement.
int REFINEMENT_COUNT = 0;

// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies (see below).
const double THRESHOLD = 0.3;

// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
CandList CAND_LIST = H2D_H_ANISO;

// Maximum polynomial degree used. -1 for unlimited.
// See User Documentation for details.
const int MAX_P_ORDER = 1;

// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;

// Stopping criterion for adaptivity.
double adaptivityErrorStop(int iteration)
{
  return 0.02;
}

// Equation parameters.
const double KAPPA = 1.4;
const double RHO_LEFT = 1.0;
const double RHO_TOP = 1.7;

const double V1_LEFT = 2.9;
const double V1_TOP = 2.619334;

const double V2_LEFT = 0.0;
const double V2_TOP = -0.5063;

const double PRESSURE_LEFT = 0.714286;
const double PRESSURE_TOP = 1.52819;

// Initial values
const double RHO_INIT = 1.0;
const double V1_INIT = 2.9;
const double V2_INIT = 0.0;
const double PRESSURE_INIT = 0.714286;

// Mesh filename.
const std::string MESH_FILENAME = "channel.mesh";
// Boundary markers.
const std::string BDY_SOLID_WALL = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_INLET_TOP = "3";
const std::string BDY_INLET_LEFT = "4";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
#pragma region 1. Load mesh and initialize spaces.
  // Load the mesh->
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load(MESH_FILENAME, mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT, new L2ShapesetTaylor));
  std::vector<SpaceSharedPtr<double> > spaces({ space_rho, space_rho_v_x, space_rho_v_y, space_e });
  int ndof = Space<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);
#pragma endregion

#pragma region 2. Initialize solutions.
  MeshFunctionSharedPtr<double> sln_rho(new Solution<double>(mesh));
  MeshFunctionSharedPtr<double> sln_rho_v_x(new Solution<double>(mesh));
  MeshFunctionSharedPtr<double> sln_rho_v_y(new Solution<double>(mesh));
  MeshFunctionSharedPtr<double> sln_e(new Solution<double>(mesh));
  std::vector<MeshFunctionSharedPtr<double> > slns({ sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e });

  MeshFunctionSharedPtr<double> rsln_rho(new Solution<double>(mesh));
  MeshFunctionSharedPtr<double> rsln_rho_v_x(new Solution<double>(mesh));
  MeshFunctionSharedPtr<double> rsln_rho_v_y(new Solution<double>(mesh));
  MeshFunctionSharedPtr<double> rsln_e(new Solution<double>(mesh));
  std::vector<MeshFunctionSharedPtr<double> > rslns({ rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e });
#pragma endregion

#pragma region 3. Filters for visualization of Mach number, pressure + visualization setup.
  MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(rslns, KAPPA));
  MeshFunctionSharedPtr<double>  pressure(new PressureFilter(rslns, KAPPA));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(650, 0, 600, 300));
  ScalarView eview("Error - density", new WinGeom(0, 330, 600, 300));
  ScalarView eview1("Error - momentum", new WinGeom(0, 660, 600, 300));
  OrderView order_view("Orders", new WinGeom(650, 330, 600, 300));
#pragma endregion

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  Vector<double>* rhs_stabilization = create_vector<double>(HermesCommonApi.get_integral_param_value(matrixSolverType));

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

  // Set initial conditions.
  MeshFunctionSharedPtr<double> prev_rho(new ConstantSolution<double>(mesh, RHO_INIT));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double>(mesh, RHO_INIT * V1_INIT));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double>(mesh, RHO_INIT * V2_INIT));
  MeshFunctionSharedPtr<double> prev_e(new ConstantSolution<double>(mesh, QuantityCalculator::calc_energy(RHO_INIT, RHO_INIT * V1_INIT, RHO_INIT * V2_INIT, PRESSURE_INIT, KAPPA)));

  // Initialize weak formulation.
  std::vector<std::string> solid_wall_markers;
  solid_wall_markers.push_back(BDY_SOLID_WALL);

  std::vector<std::string> inlet_markers({ BDY_INLET_LEFT, BDY_INLET_TOP });
  std::vector<double> rho_ext({ RHO_LEFT, RHO_TOP });
  std::vector<double> v1_ext({ V1_LEFT, V1_TOP });
  std::vector<double> v2_ext({ V2_LEFT, V2_TOP });
  std::vector<double> pressure_ext({ PRESSURE_LEFT, PRESSURE_TOP });

  std::vector<std::string> outlet_markers({ BDY_OUTLET });

  WeakFormSharedPtr<double> wf(new EulerEquationsWeakFormSemiImplicit(KAPPA, rho_ext, v1_ext, v2_ext, pressure_ext, solid_wall_markers, inlet_markers, outlet_markers, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e));
#include "../euler-time-loop-space-adapt.cpp"
}