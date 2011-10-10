#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

/* Right-hand side */

class CustomRightHandSide : public Hermes::Hermes2DFunction<double>
{
public:
  CustomRightHandSide(double poly_deg)
      : Hermes::Hermes2DFunction<double>(), poly_deg(poly_deg) {};

  virtual double value(double x, double y) const;
  virtual Ord value (Ord x, Ord y) const 
  { 
    return Ord(8); 
  }
  
  double poly_deg;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(Mesh* mesh, double poly_deg)
      : ExactSolutionScalar<double>(mesh), poly_deg(poly_deg) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord (Ord x, Ord y) const 
  { 
    return Ord(Ord::get_max_order()); 
  }

  double poly_deg;
};



