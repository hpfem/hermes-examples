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
  double a_P = (-alpha * pow((x - x_loc), 2) - alpha * pow((y - y_loc), 2));

  return (4 * exp(a_P) * alpha * (alpha * (x - x_loc) * (x - x_loc)
         + alpha * (y - y_loc) * (y - y_loc) - 1));
}

double CustomExactSolution::value(double x, double y) const 
{
  return exp(-alpha * (pow((x - x_loc), 2) + pow((y - y_loc), 2)));
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  double a = -alpha * ( (x - x_loc) * (x - x_loc) + (y - y_loc) * (y - y_loc));
  dx = -exp(a) * (2 * alpha * (x - x_loc));
  dy = -exp(a) * (2 * alpha * (y - y_loc));
}
