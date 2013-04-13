#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

/* Right-hand side */

class CustomRightHandSide : public Hermes::Hermes2DFunction<double>
{
public:
  CustomRightHandSide(double alpha)
      : Hermes::Hermes2DFunction<double>(), alpha(alpha) {};

  virtual double value(double x, double y) const;

  virtual Ord value (Ord x, Ord y) const;
  
  double alpha;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(MeshSharedPtr mesh, double alpha)
      : ExactSolutionScalar<double>(mesh), alpha(alpha) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord (Ord x, Ord y) const; 

  double alpha;
};
