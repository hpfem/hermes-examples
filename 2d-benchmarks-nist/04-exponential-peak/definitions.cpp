#include "definitions.h"

double CustomRightHandSide::value(double x, double y) const
{
  double a_P = (-alpha * Hermes::pow((x - x_loc), 2) - alpha * Hermes::pow((y - y_loc), 2));

  return -(4 * exp(a_P) * alpha * (alpha * (x - x_loc) * (x - x_loc)
         + alpha * (y - y_loc) * (y - y_loc) - 1));
}

Ord CustomRightHandSide::value (Ord x, Ord y) const
{
  return Ord(8);
}

double CustomExactSolution::value(double x, double y) const 
{
  return exp(-alpha * (Hermes::pow((x - x_loc), 2) + Hermes::pow((y - y_loc), 2)));
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  double a = -alpha * ( (x - x_loc) * (x - x_loc) + (y - y_loc) * (y - y_loc));
  dx = -exp(a) * (2 * alpha * (x - x_loc));
  dy = -exp(a) * (2 * alpha * (y - y_loc));
}

Ord CustomExactSolution::ord (Ord x, Ord y) const
{
  return Ord(8);
}
