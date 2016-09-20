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
// Domain: A square, see file square.mesh->
//
// BC: Solid walls, inlet, no outlet.
//
// IC: Constant state identical to inlet, only with higher pressure.
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
bool SHOCK_CAPTURING = false;
shockCapturingType SHOCK_CAPTURING_TYPE = KUZMIN;
// Quantitative parameter of the discontinuity detector in case of Krivodonova.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;
// Quantitative parameter of the shock capturing in case of Feistauer.
const double NU_1 = 0.1;
const double NU_2 = 0.1;

// Initial polynomial degree.
// Do not change.
const int P_INIT = 0;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;
// CFL value.
double CFL_NUMBER = 0.1;
// Initial time step.
double time_step_n = 1E-6;
double TIME_INTERVAL_LENGTH = 20.;

// Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_FREQ = 5;

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
CandList CAND_LIST = H2D_HP_ANISO;

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
  if (iteration > 49)
    return 0.01;
  else
    return 0.5 - iteration * 0.01;
}

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 2.0;
// Initial pressure (dimensionless).
const double P_INITIAL_HIGH = 1.5;
// Initial pressure (dimensionless).
const double P_INITIAL_LOW = 1.0;
// Inlet density (dimensionless).
const double RHO_EXT = 1.0;
// Initial density (dimensionless).
const double RHO_INITIAL_HIGH = 0.5;
// Initial density (dimensionless).
const double RHO_INITIAL_LOW = 0.3;
// Inlet x-velocity (dimensionless).
const double V1_EXT = 0.0;
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;
// Kappa.
const double KAPPA = 1.4;

const double MESH_SIZE = 3.0;

// Mesh filename.
const std::string MESH_FILENAME = "square.mesh";
// Boundary markers.
const std::string BDY_INLET = "Inlet";
const std::string BDY_SOLID_WALL = "Solid";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
#include "../euler-init-main-adapt.cpp"

  // Set initial conditions.
  MeshFunctionSharedPtr<double> prev_rho(new InitialSolutionLinearProgress(mesh, RHO_INITIAL_HIGH, RHO_INITIAL_LOW, MESH_SIZE));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double>(mesh, 0.0));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double>(mesh, 0.0));
  MeshFunctionSharedPtr<double> prev_e(new InitialSolutionLinearProgress(mesh, QuantityCalculator::calc_energy(RHO_INITIAL_HIGH, RHO_INITIAL_HIGH * V1_EXT, RHO_INITIAL_HIGH * V2_EXT, P_INITIAL_HIGH, KAPPA), QuantityCalculator::calc_energy(RHO_INITIAL_LOW, RHO_INITIAL_LOW * V1_EXT, RHO_INITIAL_LOW * V2_EXT, P_INITIAL_LOW, KAPPA), MESH_SIZE));

  // Initialize weak formulation.
  std::vector<std::string> solid_wall_markers({ BDY_SOLID_WALL });
  std::vector<std::string> inlet_markers({ BDY_INLET });
  std::vector<std::string> outlet_markers;

  WeakFormSharedPtr<double> wf(new EulerEquationsWeakFormSemiImplicit(KAPPA, {RHO_EXT}, {V1_EXT}, {V2_EXT}, {P_EXT}, solid_wall_markers,
    inlet_markers, outlet_markers, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e));
#include "../euler-time-loop-space-adapt.cpp"
}