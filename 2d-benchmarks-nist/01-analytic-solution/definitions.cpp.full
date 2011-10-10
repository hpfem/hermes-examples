#include "definitions.h"

double CustomWeakForm::Jacobian::value(int n, double* wt, 
                                       Func< double >* u_ext[], Func< double >* u, Func< double >* v, 
                                       Geom< double >* e, ExtData< double >* ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  
  return result;
}

Ord CustomWeakForm::Jacobian::ord(int n, double* wt, 
                                  Func< Ord >* u_ext[], Func< Ord >* u, Func< Ord >* v, 
                                  Geom< Ord >* e, ExtData< Ord >* ext) const
{ 
  return u->dx[0] * v->dx[0] + u->dy[0] * v->dy[0];
}

double CustomWeakForm::Residual::value(int n, double* wt, Func< double >* u_ext[], Func< double >* v,
                                       Geom< double >* e, ExtData< double >* ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ( u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i] + rhs->value(e->x[i], e->y[i]) * v->val[i] );
  
  return result;
}

Ord CustomWeakForm::Residual::ord(int n, double* wt, Func< Ord >* u_ext[], Func< Ord >* v, 
                                  Geom< Ord >* e, ExtData< Ord >* ext) const
{ 
  return u_ext[0]->dx[0] * v->dx[0] + u_ext[0]->dy[0] * v->dy[0] + rhs->value(e->x[0], e->y[0]) * v->val[0];
}

double CustomRightHandSide::value(double x, double y) const
{
  double a = pow(2.0, 4.0*poly_deg);
  double b = pow(x-1.0, 8.0);
  double c = (38.0*pow(x, 2.0) - 38.0*x + 9.0);
  double d = pow(y-1.0, poly_deg);
  double e = pow(y-1.0, 8.0);
  double f = (38.0*pow(y, 2.0) - 38.0*y + 9.0);
  double g = pow(x-1.0, poly_deg);

  return (poly_deg*a*pow(x, 8.0)*b*c*pow(y, poly_deg)*d
	       + poly_deg*a*pow(y, 8.0)*e*f*pow(x, poly_deg)*g);
}

double CustomExactSolution::value(double x, double y) const 
{
  return pow(2, 4 * poly_deg) * pow(x, poly_deg) * pow(1 - x, poly_deg)
         * pow(y, poly_deg) * pow(1 - y, poly_deg);
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  double A = pow((1.0-y), poly_deg);
  double B = pow((1.0-x), poly_deg);
  double C = pow(y, poly_deg);
  double D = pow(x, poly_deg);

  dx = ((poly_deg * pow(16.0, poly_deg)*A*C) / (x-1.0)
       + (poly_deg*pow(16.0, poly_deg)*A*C)/x)*B*D;
  dy = ((poly_deg*pow(16.0, poly_deg)*B*D)/(y-1.0)
       + (poly_deg*pow(16.0, poly_deg)*B*D)/y)*A*C;
}
