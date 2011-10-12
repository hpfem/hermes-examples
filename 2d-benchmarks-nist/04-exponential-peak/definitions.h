#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

/* Right-hand side */

class CustomRightHandSide : public Hermes::Hermes2DFunction<double>
{
public:
  CustomRightHandSide(double alpha, double x_loc, double y_loc)
      : Hermes::Hermes2DFunction<double>(), alpha(alpha), x_loc(x_loc), y_loc(y_loc) {};

  virtual double value(double x, double y) const;

  virtual Ord value (Ord x, Ord y) const;
  
  double alpha;
  double x_loc;
  double y_loc;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(Mesh* mesh, double alpha, double x_loc, double y_loc)
      : ExactSolutionScalar<double>(mesh), alpha(alpha), x_loc(x_loc), y_loc(y_loc) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord (Ord x, Ord y) const; 

  double alpha;
  double x_loc;
  double y_loc;
};
