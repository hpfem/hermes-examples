#include "definitions.h"

//  This example shows how to set an arbitrary initial guess for the
//  Newton's method, and nonzero Dirichlet boundary conditions.
//
//  PDE: magnetostatics with nonlinear magnetic permeability
//  curl[1/mu curl u] = current_density.
//
//  Domain: unit square (-10,10)^2.
//
//  BC: Dirichlet, see function dir_lift() below.
//
//  The following parameters can be changed:

// Initial polynomial degree.
const int P_INIT = 2;
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-5;
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 1000;
// Number between 0 and 1 to damp Newton's iterations.
const double NEWTON_DAMPING = 1.0;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;

// Problem parameters.
double MU_VACUUM = 4. * M_PI * 1e-7;
// Constant initial condition for the magnetic potential.
double INIT_COND = 0.0;
// Volume source term.
double CURRENT_DENSITY = 1e6;

// Material and boundary markers.
const std::string MAT_AIR = "2";
const std::string MAT_IRON_1 = "0";
const std::string MAT_IRON_2 = "3";
const std::string MAT_COPPER = "1";
const std::vector<std::string> BDY_DIRICHLET = { 
"16",
"17",
"18",
"19",
"21",
"22"
};

int main(int argc, char* argv[])
{
  // Define nonlinear magnetic permeability via a cubic spline.
  std::vector<double> mu_inv_pts({ 0.,
	  0.0622676,
	  0.0699801,
	  0.0781449,
	  0.108779,
	  0.140252,
	  0.589805,
	  0.777016,
	  0.796086,
	  0.895488,
	  1.01564,
	  1.11306,
	  1.15599,
	  1.19261,
	  1.27821,
	  1.37143,
	  1.47334,
	  1.50096,
	  1.51305,
	  1.53748,
	  1.53964,
	  1.56167,
	  1.57542,
	  1.57922,
	  1.60284,
	  1.61321,
	  1.62708,
	  1.70197,
	  1.73799,
	  1.81436,
	  1.91181,
	  1.9277,
	  1.94428 });
  std::vector<double> mu_inv_val({ 14617.8,
	  18943.1,
	  20424.2,
	  21880.9,
	  26343.9,
	  30625.4,
	  63874.2,
	  59235.4,
	  58231.,
	  51103.9,
	  38703.5,
	  26056.4,
	  20664.7,
	  16448.8,
	  7526.81,
	  2687.02,
	  688.82,
	  435.214,
	  334.89,
	  180.578,
	  166.409,
	  61.6026,
	  40.5858,
	  37.0509,
	  24.3051,
	  21.1503,
	  18.0643,
	  10.3419,
	  8.66857,
	  6.55724,
	  5.11165,
	  4.94426,
	  4.78435 });

  // Create the cubic spline (and plot it for visual control).
  double bc_left = 0.0;
  double bc_right = 0.0;
  bool first_der_left = false;
  bool first_der_right = false;
  bool extrapolate_der_left = false;
  bool extrapolate_der_right = false;
  CubicSpline mu_inv_iron(mu_inv_pts, mu_inv_val, bc_left, bc_right, first_der_left, first_der_right,
    extrapolate_der_left, extrapolate_der_right);
  mu_inv_iron.calculate_coeffs();

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  std::vector<MeshSharedPtr> meshes({ mesh });
  mloader.load("sit.msh", meshes);

  // View the mesh.
  MeshView m_view;
  m_view.show(mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_essential("16", 0.0);
  DefaultEssentialBCConst<double> bc_essential1("17", 0.0);
  DefaultEssentialBCConst<double> bc_essential2("18", 0.0);
  DefaultEssentialBCConst<double> bc_essential3("19", 0.0);
  DefaultEssentialBCConst<double> bc_essential4("21", 0.0);
  DefaultEssentialBCConst<double> bc_essential5("22", 0.0);
  EssentialBCs<double> bcs({ &bc_essential, &bc_essential1, &bc_essential2 });

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();
  Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

  Hermes::Hermes2D::Views::BaseView<double> b;
  b.show(space);

  // Initialize the weak formulation
  // This determines the increase of integration order
  // for the axisymmetric term containing 1/r. Default is 3.
  int order_inc = 3;
  WeakFormSharedPtr<double> wf(new CustomWeakFormMagnetostatics(MAT_IRON_1, MAT_IRON_2, &mu_inv_iron, MAT_AIR, MAT_COPPER, MU_VACUUM, CURRENT_DENSITY, order_inc));

  // Initialize the FE problem.
  DiscreteProblem<double> dp(wf, space);

  // Initialize the solution.
  MeshFunctionSharedPtr<double> sln(new ConstantSolution<double>(mesh, INIT_COND));

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  Hermes::Mixins::Loggable::Static::info("Projecting to obtain initial vector for the Newton's method.");

  // Perform Newton's iteration.
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp);
  newton.set_manual_damping_coeff(true, 0.1);
  newton.set_sufficient_improvement_factor_jacobian(0.5);
  newton.set_max_steps_with_reused_jacobian(5);
  newton.set_max_allowed_iterations(NEWTON_MAX_ITER);
  newton.set_tolerance(NEWTON_TOL, Hermes::Solvers::ResidualNormAbsolute);
  
  try
  {
    newton.solve(sln);
  }
  catch (Hermes::Exceptions::Exception e)
  {
    e.print_msg();
    throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    return 1;
  };

  // Translate the resulting coefficient vector into the Solution sln.
  Solution<double>::vector_to_solution(newton.get_sln_vector(), space, sln);

  // Visualise the solution and mesh.
  ScalarView s_view1("Vector potential", new WinGeom(0, 0, 350, 450));
  MeshFunctionSharedPtr<double> vector_potential(new FilterVectorPotential(std::vector<MeshFunctionSharedPtr<double> >({ sln, sln }), { H2D_FN_VAL, H2D_FN_VAL }));
  s_view1.show(sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}