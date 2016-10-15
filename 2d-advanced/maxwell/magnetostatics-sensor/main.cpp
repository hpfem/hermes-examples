#include "definitions.h"

//  This example shows how to handle stiff nonlinear problems.
//
//  PDE: magnetostatics with nonlinear magnetic permeability
//  curl[1/mu curl u] = current_density.

//  The following parameters can be changed:

// Initial polynomial degree.
const int P_INIT = 2;
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-9;
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 1000;

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
const std::vector<std::string> BDY_DIRICHLET = { "0","16","17","18","19","21","22","23" };

int main(int argc, char* argv[])
{
	// Define nonlinear magnetic permeability via a cubic spline.
	std::vector<double> mu_inv_pts({ 0.,
		0.0622676,		0.0699801,		0.0781449,		0.108779,		0.140252,		0.589805,		0.777016,		0.796086,
		0.895488,		1.01564,		1.11306,		1.15599,		1.19261,		1.27821,		1.37143,		1.47334,
		1.50096,		1.51305,		1.53748,		1.53964,		1.56167,		1.57542,		1.57922,		1.60284,
		1.61321,		1.62708,		1.70197,		1.73799,		1.81436,		1.91181,		1.9277,		1.94428 });
	std::vector<double> mu_inv_val({ 0.0000684097470207555,
		0.0000527896701173514,		0.0000489615260328434,		0.0000457019592429927,		0.0000379594517136794,
		0.0000326526347410973,		0.0000156557733795492,		0.0000168817970335306,		0.0000171729834624169,
		0.0000195679781777907,		0.0000258374565607762,		0.0000383782871002901,		0.0000483917017909769,
		0.0000607947084285784,		0.0001328584088079810,		0.0003721594926721800,		0.0014517580790337100,
		0.0022977202020155600,		0.0029860551225775600,		0.0055377731506606600,		0.0060092903629010500,
		0.0162330810712535000,		0.0246391595089908000,		0.0269898976812979000,		0.0411436282920046000,
		0.0472806532295050000,		0.0553578051737405000,		0.0966940310774616000,		0.1153592807118130000,
		0.1525031873166150000,		0.1956315475433570000,		0.2022547357946390000,		0.2090148086991960000
	});

	// Create the cubic spline (and plot it for visual control).
	double bc_left = 0.;
	double bc_right = 0.;
	bool first_der_left = false;
	bool first_der_right = false;
	bool extrapolate_der_left = true;
	bool extrapolate_der_right = true;
	CubicSpline mu_inv_iron(mu_inv_pts, mu_inv_val, bc_left, bc_right, first_der_left, first_der_right,
		extrapolate_der_left, extrapolate_der_right);
	mu_inv_iron.calculate_coeffs();

	// Load the mesh.
	MeshSharedPtr mesh(new Mesh);
	MeshReaderH2DXML mloader;
	mloader.load("sensor.msh", mesh);

	// View the mesh.
	MeshView m_view("Mesh", new WinGeom(50, 0, 400, 400));
	m_view.show(mesh);

	// Initialize boundary conditions.
	DefaultEssentialBCConst<double> bc_essential(BDY_DIRICHLET, 0.0);
	EssentialBCs<double> bcs({ &bc_essential });

	// Create an H1 space with default shapeset.
	SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
	int ndof = space->get_num_dofs();
	Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);

	// Initialize the weak formulation
	int order_inc = 3;
	WeakFormSharedPtr<double> wf(new CustomWeakFormMagnetostatics(MAT_IRON_1, MAT_IRON_2, &mu_inv_iron, MAT_AIR, MAT_COPPER, MU_VACUUM, CURRENT_DENSITY, order_inc));

	// Initialize the FE problem.
	DiscreteProblem<double> dp(wf, space);

	// Initialize the solution.
	MeshFunctionSharedPtr<double> sln(new Solution<double>());

	// Perform Newton's iteration.
	Hermes::Hermes2D::NewtonSolver<double> newton(&dp);
	newton.set_sufficient_improvement_factor(1.5);
	newton.set_necessary_successful_steps_to_increase(1);
	newton.set_max_steps_with_reused_jacobian(0);
	newton.set_max_allowed_iterations(NEWTON_MAX_ITER);
	newton.set_tolerance(NEWTON_TOL, Hermes::Solvers::ResidualNormAbsolute);

	newton.solve();

	// Translate the resulting coefficient vector into the Solution sln.
	Solution<double>::vector_to_solution(newton.get_sln_vector(), space, sln);

	// Visualise the solution and mesh.
	ScalarView Sview1("Vector Potential", new WinGeom(450, 0, 400, 400));
	Sview1.show(sln);
	ScalarView Sview2("Flux Density", new WinGeom(850, 0, 400, 400));
	MeshFunctionSharedPtr<double> flux_density(new FilterFluxDensity(std::vector<MeshFunctionSharedPtr<double> >({ sln, sln })));
	Sview2.show(flux_density);

	// Wait for all views to be closed.
	View::wait();
	return 0;
}