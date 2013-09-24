#include "hermes2d.h"
#include "../NIST-util.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(MeshSharedPtr mesh, double alpha)
      : ExactSolutionScalar<double>(mesh), alpha(alpha) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord (double x, double y) const;

  double get_angle(double y, double x) const;

  MeshFunction<double>* clone() const { return new CustomExactSolution(mesh, alpha); }

  double alpha;
};



