#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

const int P_INIT = 2;

const double time_step = 4e-5;
const double end_time = 1.;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  std::vector<MeshSharedPtr> meshes({ mesh });
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mloader.load("domain.mesh", meshes);

  MeshView m;
  m.show(mesh);

  std::vector<std::string> soft_boundaries({ "2" });

  CustomBCValue custom_bc_value({ "3", "4", "5", "6" }, 1., 1000.);
  DefaultEssentialBCConst<double> default_bc_value(soft_boundaries, 0.0);
  Hermes::Hermes2D::EssentialBCs<double> bcs_value({ &custom_bc_value, &default_bc_value });

  CustomBCDerivative custom_bc_derivative({ "3", "4", "5", "6" }, 1., 1000.);
  DefaultEssentialBCConst<double> default_bc_derivative(soft_boundaries, 0.0);
  Hermes::Hermes2D::EssentialBCs<double> bcs_derivative({ &custom_bc_derivative, &default_bc_derivative });

  SpaceSharedPtr<double> space_value(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs_value, P_INIT));
  SpaceSharedPtr<double> space_derivative(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs_derivative, P_INIT));
  std::vector<SpaceSharedPtr<double> > spaces({ space_value, space_derivative });

  BaseView<double> b;
  b.show(space_value);

  // Initialize the solution.
  MeshFunctionSharedPtr<double> sln_value(new ZeroSolution<double>(mesh));
  MeshFunctionSharedPtr<double> sln_derivative(new ZeroSolution<double>(mesh));
  std::vector<MeshFunctionSharedPtr<double> > slns({ sln_value, sln_derivative });

  Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));
  //viewS.set_min_max_range(-1.0, 1.0);

  WeakFormSharedPtr<double> wf(new MyWeakForm({ "0", "1", "7" }, { 345 * 1.25, 345 * 1.25, 345 * 1.25 }, slns));
  wf->set_current_time_step(time_step);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(wf, spaces);
  linear_solver.set_jacobian_constant();

  // Solve the linear problem.
  unsigned int iteration = 0;
  for (double time = time_step; time < end_time; time += time_step)
  {
    try
    {
      Hermes::Mixins::Loggable::Static::info("Iteration: %u, Time: %g.", iteration, time);
      Space<double>::update_essential_bc_values(spaces, time);
      linear_solver.solve();

      // Get the solution vector.
      double* sln_vector = linear_solver.get_sln_vector();

      // Translate the solution vector into the previously initialized Solution.
      Hermes::Hermes2D::Solution<double>::vector_to_solutions(sln_vector, spaces, slns);

      // Visualize the solution.
      viewS.show(slns[0]);

      //viewS.wait_for_keypress();
    }

    catch (std::exception& e)
    {
      std::cout << e.what();
      return -1;
    }

    iteration++;
  }

  return 0;
}