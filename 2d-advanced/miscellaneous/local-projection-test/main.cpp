#define HERMES_REPORT_ALL
#include "definitions.h"

// This example shows how to solve a simple PDE that describes stationary 
// heat transfer in an object consisting of two materials (aluminum and 
// copper). The object is heated by constant volumetric heat sources
// (generated, for example, by a DC electric current). The temperature 
// on the boundary is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - select a matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format 
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh). 
//
// The following parameters can be changed:

// Read the original or XML mesh file.
const bool USE_XML_FORMAT = false;                 
// Set to "false" to suppress Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;           
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;              
// Uniform polynomial degree of mesh elements.
const int P_INIT = 5;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;                       
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

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
  if (USE_XML_FORMAT == true)
  {
    MeshReaderH2DXML mloader;  
    info("Reading mesh in XML format.");
    mloader.load("domain.xml", &mesh);
  }
  else 
  {
    MeshReaderH2D mloader;
    info("Reading mesh in original format.");
    mloader.load("domain.mesh", &mesh);
  }

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) 
      mesh.refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Aluminum", new Hermes1DFunction<double>(LAMBDA_AL), 
      "Copper", new Hermes1DFunction<double>(LAMBDA_CU), 
      new Hermes2DFunction<double>(-VOLUME_HEAT_SRC));
  
  // Initialize essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential(
      Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"), FIXED_BDY_TEMP);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Initialize Newton solver.
  NewtonSolver<double> newton(&dp, matrix_solver);

  // Perform Newton's iteration.
  try
  {
    newton.solve();
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.printMsg();
    error("Newton's iteration failed.");
  }

  // Translate the resulting coefficient vector into a Solution.
  Solution<double> sln;
  Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);

  Solution<double> sln_proj;
  LocalProjection<double>::project_local(&space, &sln, &sln_proj, matrix_solver, HERMES_H1_NORM);
  ScalarView view1("Projection", new WinGeom(0, 0, 440, 350));
  view1.show(&sln_proj, HERMES_EPS_HIGH);
  View::wait();

  // VTK output.
  if (VTK_VISUALIZATION) 
  {
    // Output solution in VTK format.
    Linearizer lin;
    bool mode_3D = true;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D);
    info("Solution in VTK format saved to file %s.", "sln.vtk");

    // Output mesh and element orders in VTK format.
    Orderizer ord;
    ord.save_orders_vtk(&space, "ord.vtk");
    info("Element orders in VTK format saved to file %s.", "ord.vtk");
  }

  // Visualize the solution.
  if (HERMES_VISUALIZATION) 
  {
    ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
    // Hermes uses adaptive FEM to approximate higher-order FE solutions with linear
    // triangles for OpenGL. The second parameter of View::show() sets the error 
    // tolerance for that. Options are HERMES_EPS_LOW, HERMES_EPS_NORMAL (default), 
    // HERMES_EPS_HIGH and HERMES_EPS_VERYHIGH. The size of the graphics file grows 
    // considerably with more accurate representation, so use it wisely.
    view.show(&sln, HERMES_EPS_HIGH);
    View::wait();
  }

  return 0;
}

