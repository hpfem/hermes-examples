#include "hermes2d.h"

using namespace Hermes::Hermes2D;
using namespace WeakFormsH1;
using Hermes::Ord;

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(Mesh* mesh, double slope)
            : ExactSolutionScalar<double>(mesh), slope(slope) {};

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double slope;
};

/* Custom function f */

class CustomFunction: public HermesFunction<double>
{
public:
  CustomFunction(double slope)
    : HermesFunction<double>(), slope(slope) {};

  virtual double value(double x, double y) const;

  virtual Ord value_ord(Ord x, Ord y) const;

  double slope;
};

