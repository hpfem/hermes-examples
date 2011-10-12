#include "definitions.h"

double CustomRightHandSide::value(double x, double y) const
{
  double a = Hermes::pow(2.0, 4.0*poly_deg);
  double b = Hermes::pow(x-1.0, 8.0);
  double c = (38.0*Hermes::pow(x, 2.0) - 38.0*x + 9.0);
  double d = Hermes::pow(y-1.0, poly_deg);
  double e = Hermes::pow(y-1.0, 8.0);
  double f = (38.0*Hermes::pow(y, 2.0) - 38.0*y + 9.0);
  double g = Hermes::pow(x-1.0, poly_deg);

  return (poly_deg*a*Hermes::pow(x, 8.0)*b*c*Hermes::pow(y, poly_deg)*d
	       + poly_deg*a*Hermes::pow(y, 8.0)*e*f*Hermes::pow(x, poly_deg)*g);
}

Ord CustomRightHandSide::value (Ord x, Ord y) const
{
  return Ord(8);
}

double CustomExactSolution::value(double x, double y) const 
{
  return Hermes::pow(2, 4 * poly_deg) * Hermes::pow(x, poly_deg) * Hermes::pow(1 - x, poly_deg)
         * Hermes::pow(y, poly_deg) * Hermes::pow(1 - y, poly_deg);
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  double A = Hermes::pow((1.0-y), poly_deg);
  double B = Hermes::pow((1.0-x), poly_deg);
  double C = Hermes::pow(y, poly_deg);
  double D = Hermes::pow(x, poly_deg);

  dx = ((poly_deg * Hermes::pow(16.0, poly_deg)*A*C) / (x-1.0)
       + (poly_deg*Hermes::pow(16.0, poly_deg)*A*C)/x)*B*D;
  dy = ((poly_deg*Hermes::pow(16.0, poly_deg)*B*D)/(y-1.0)
       + (poly_deg*Hermes::pow(16.0, poly_deg)*B*D)/y)*A*C;
}

Ord CustomExactSolution::ord (Ord x, Ord y) const
{
  return Ord(8);
}
