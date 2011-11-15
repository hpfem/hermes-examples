#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace Hermes::Hermes2D;

// This example is analogous to P01-linear/03-poisson.
//
// PDE: Poisson equation -(LAMBDA u')' - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: Interval (0, 2*pi).
//
// The following parameters can be changed:

// Set to "false" to suppress Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;           
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;              
// Uniform polynomial degree of mesh elements.
const int P_INIT = 5;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 3;                       

// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  

// Problem parameters.
// Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_AL = 236.0;            
// Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            
// Volume heat sources generated (for example) by electric current.  
const double VOLUME_HEAT_SRC = 5e3;              
// Fixed temperature on the boundary.
const double FIXED_BDY_TEMP = 20.0;        

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH1DXML mloader;
  mloader.load("domain.xml", &mesh);

  // Perform initial mesh refinements (optional).
  // Split elements vertically.
  int refinement_type = 2;            
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(refinement_type);

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Al", new Hermes::Hermes1DFunction<double>(LAMBDA_AL), "Cu", 
      new Hermes::Hermes1DFunction<double>(LAMBDA_CU), 
      new Hermes::Hermes2DFunction<double>(-VOLUME_HEAT_SRC));

  // Initialize essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Left", "Right"), 
      FIXED_BDY_TEMP);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);
  
  // Show the mesh and poly degrees.
  Views::OrderView oview("Mesh", new Views::WinGeom(0, 0, 900, 250));
  if (HERMES_VISUALIZATION) oview.show(&space);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Solution<double> sln;
  NewtonSolver<double> newton(&dp, matrix_solver);
  try
  {
    newton.solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.printMsg();
    error("Newton's iteration failed.");
  };
  Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);

  // Get info about time spent during assembling in its respective parts.
  //dp.get_all_profiling_output(std::cout);

  // VTK output.
  if (VTK_VISUALIZATION)
  {
    // Output solution in VTK format.
    Views::Linearizer lin;
    bool mode_3D = true;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D);
    info("Solution in VTK format saved to file %s.", "sln.vtk");

    // Output mesh and element orders in VTK format.
    Views::Orderizer ord;
    ord.save_orders_vtk(&space, "ord.vtk");
    info("Element orders in VTK format saved to file %s.", "ord.vtk");
  }

  // Visualize the solution.
  if (HERMES_VISUALIZATION)
  {
    Views::ScalarView view("Solution", new Views::WinGeom(0, 300, 900, 350));
    // Hermes uses adaptive FEM to approximate higher-order FE solutions with linear
    // triangles for OpenGL. The second parameter of View::show() sets the error 
    // tolerance for that. Options are HERMES_EPS_LOW, HERMES_EPS_NORMAL (default), 
    // HERMES_EPS_HIGH and HERMES_EPS_VERYHIGH. The size of the graphics file grows 
    // considerably with more accurate representation, so use it wisely.
    view.show(&sln, Views::HERMES_EPS_HIGH);
    Views::View::wait();
  }

  return 0;
}

