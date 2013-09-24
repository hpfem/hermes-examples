#include "definitions.h"

double CustomExactFunction::fn(double x, double y) 
{
  if (x <= 0) return Hermes::cos(k * y);
  else return Hermes::cos(k * y) + Hermes::pow(x, alpha);
}

CustomRightHandSide::CustomRightHandSide(double k, double alpha) : Hermes::Hermes2DFunction<double>(), k(k), alpha(alpha)
{
  cef = new CustomExactFunction(k, alpha);
}

CustomRightHandSide::~CustomRightHandSide()
{
  delete cef;
}

double CustomRightHandSide::value(double x, double y) const
{
  if (x < 0) return -cef->fn(x, y) * k * k;
  else return -(cef->fn(x, y) * k * k
              - alpha *(alpha - 1) * Hermes::pow(x, alpha - 2.)
              - k * k * Hermes::pow(x, alpha));
}

Ord CustomRightHandSide::value(Ord x, Ord y) const
{
  return Ord(16);
}

CustomExactSolution::CustomExactSolution(MeshSharedPtr mesh, double k, double alpha) : ExactSolutionScalar<double>(mesh), k(k), alpha(alpha)
{
  cef = new CustomExactFunction(k, alpha);
}

CustomExactSolution::~CustomExactSolution() 
{ 
  delete cef; 
}

double CustomExactSolution::value(double x, double y) const 
{
  return cef->fn(x, y);
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  if (x <= 0) dx = 0;
  else dx = alpha * Hermes::pow(x, alpha - 1);
  dy = -Hermes::sin(k * y) * k;
}

Ord CustomExactSolution::ord(double x, double y) const
{
  return Ord(16);
}
