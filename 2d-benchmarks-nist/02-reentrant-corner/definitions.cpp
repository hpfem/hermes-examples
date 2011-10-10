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
    result += wt[i] * ( u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
  
  return result;
}

Ord CustomWeakForm::Residual::ord(int n, double* wt, Func< Ord >* u_ext[], Func< Ord >* v, 
                                  Geom< Ord >* e, ExtData< Ord >* ext) const
{ 
  return u_ext[0]->dx[0] * v->dx[0] + u_ext[0]->dy[0] * v->dy[0];
}

double CustomExactSolution::value(double x, double y) const 
{
  return (pow(sqrt(x*x + y*y), alpha) * sin(alpha * get_angle(y, x)));
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  double a = sqrt(x*x + y*y);
  double b = pow(a, (alpha - 1.0));
  double c = pow(a, alpha);
  double d = ((y*y)/(x*x) + 1.0 );

  dx = (((alpha* x* sin(alpha * get_angle(y,x)) *b)/a)
       - ((alpha *y *cos(alpha * get_angle(y, x)) * c)/(pow(x, 2.0) *d)));
  dy = (((alpha* cos(alpha* get_angle(y, x)) *c)/(x * d))
       + ((alpha* y* sin(alpha* get_angle(y, x)) *b)/a));
}


