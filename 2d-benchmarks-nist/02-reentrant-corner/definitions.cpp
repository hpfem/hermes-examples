#include "definitions.h"

double CustomExactSolution::value(double x, double y) const 
{
  return (Hermes::pow(Hermes::sqrt(x*x + y*y), alpha) * Hermes::sin(alpha * get_angle(y, x)));
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  double a = Hermes::sqrt(x*x + y*y);
  double b = Hermes::pow(a, (alpha - 1.0));
  double c = Hermes::pow(a, alpha);
  double d = ((y*y)/(x*x) + 1.0 );

  dx = (((alpha * x * Hermes::sin(alpha * get_angle(y,x)) * b)/a)
       - ((alpha * y * Hermes::cos(alpha * get_angle(y, x)) * c)/(Hermes::pow(x, 2.0) *d)));
  dy = (((alpha * Hermes::cos(alpha * get_angle(y, x)) * c)/(x * d))
       + ((alpha * y * Hermes::sin(alpha* get_angle(y, x)) * b)/a));
}

Ord CustomExactSolution::ord (double x, double y) const
{
  return Ord(10);
}

double CustomExactSolution::get_angle(double y, double x) const 
{
  double theta = Hermes::atan2(y, x);
  if (theta < 0)
    theta += 2 * M_PI;
  return theta;
}
