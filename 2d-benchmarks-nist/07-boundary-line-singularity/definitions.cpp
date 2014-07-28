#include "definitions.h"

double CustomRightHandSide::value(double x, double y) const
{
  return alpha * (alpha - 1.) * Hermes::pow(x, alpha - 2.);
}

Ord CustomRightHandSide::value(Ord x, Ord y) const
{
  return Ord((int)(alpha + 3.1));
}

double CustomExactSolution::value(double x, double y) const
{
  return Hermes::pow(x, alpha);
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = alpha * Hermes::pow(x, alpha - 1.);
  dy = 0;
}

Ord CustomExactSolution::ord(double x, double y) const
{
  return Ord((int)(alpha + 1));
}