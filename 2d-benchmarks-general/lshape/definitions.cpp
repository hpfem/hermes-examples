#include "definitions.h"

double CustomExactSolution::value(double x, double y) const 
{
  double r = Hermes::sqrt(x*x + y*y);
  double a = Hermes::atan2(x, y);
  return Hermes::pow(r, 2.0/3.0) * Hermes::sin(2.0*a/3.0 + M_PI/3);
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const 
{
  double t1 = 2.0/3.0*Hermes::atan2(x, y) + M_PI/3;
  double t2 = Hermes::pow(x*x + y*y, 1.0/3.0);
  double t3 = x*x * ((y*y)/(x*x) + 1);
  dx = 2.0/3.0*x*Hermes::sin(t1)/(t2*t2) + 2.0/3.0*y*t2*Hermes::cos(t1)/t3;
  dy = 2.0/3.0*y*Hermes::sin(t1)/(t2*t2) - 2.0/3.0*x*t2*Hermes::cos(t1)/t3;
}

Ord CustomExactSolution::ord(double x, double y) const 
{
  return Ord(20);
}
