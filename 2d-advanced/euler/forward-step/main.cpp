#define HERMES_REPORT_INFO

#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// This example solves the compressible Euler equations using a basic
// piecewise-constant finite volume method, or Discontinuous Galerkin method of higher order with no adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: forward facing step, see mesh file ffs.mesh
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
// Do not change this.
const int P_INIT = 1;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;
// Number of initial localized mesh refinements.
const int INIT_REF_NUM_STEP = 2;
// CFL value.
double CFL_NUMBER = 0.25;
// Initial time step.
double time_step_n = 1E-6;

// Equation parameters.
// Exterior pressure (dimensionless).
const double P_EXT = 1.0;
// Inlet density (dimensionless).
const double RHO_EXT = 1.4;
// Inlet x-velocity (dimensionless).
const double V1_EXT = 3.0;
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;
// Kappa.
const double KAPPA = 1.4;

double TIME_INTERVAL_LENGTH = 20.;

// Mesh filename.
const std::string MESH_FILENAME = "ffs.mesh";
// Boundary markers.
const std::string BDY_SOLID_WALL_BOTTOM = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_TOP = "3";
const std::string BDY_INLET = "4";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

// Criterion for mesh refinement.
int refinement_criterion(Element* e)
{
  if (e->vn[2]->y <= 0.4 && e->vn[1]->x <= 0.6)
    return 0;
  else
    return -1;
}

int main(int argc, char* argv[])
{
#include "../euler-init-main.cpp"

  // Set initial conditions.
  MeshFunctionSharedPtr<double> prev_rho(new ConstantSolution<double>(mesh, RHO_EXT));
  MeshFunctionSharedPtr<double> prev_rho_v_x(new ConstantSolution<double>(mesh, RHO_EXT * V1_EXT));
  MeshFunctionSharedPtr<double> prev_rho_v_y(new ConstantSolution<double>(mesh, RHO_EXT * V2_EXT));
  MeshFunctionSharedPtr<double> prev_e(new ConstantSolution<double>(mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA)));

  // Initialize weak formulation.
  std::vector<std::string> solid_wall_markers(BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP);
  std::vector<std::string> inlet_markers;
  inlet_markers.push_back(BDY_INLET);
  std::vector<std::string> outlet_markers;
  outlet_markers.push_back(BDY_OUTLET);

  EulerEquationsWeakFormSemiImplicit wf(KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, solid_wall_markers,
    inlet_markers, outlet_markers, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, (P_INIT == 0));

#include "../euler-time-loop.cpp"
}