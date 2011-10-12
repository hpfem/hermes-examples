#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

class CustomExactFunction 
{
public:
  CustomExactFunction(double k, double alpha): k(k), alpha(alpha) {};

  double fn(double x, double y);

  double k;
  double alpha;
};

/* Right-hand side */

class CustomRightHandSide : public Hermes::Hermes2DFunction<double>
{
public:
  CustomRightHandSide(double k, double alpha);

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  ~CustomRightHandSide();

  CustomExactFunction* cef;

  double k;
  double alpha;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(Mesh* mesh, double k, double alpha);

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  ~CustomExactSolution();

  CustomExactFunction* cef;

  double k;
  double alpha;
};
