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

MatrixFormVol<double>* CustomWeakForm::Jacobian::clone()
{
  return new CustomWeakForm::Jacobian(*this);
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

VectorFormVol<double>* CustomWeakForm::Residual::clone()
{
  return new CustomWeakForm::Residual(*this);
}

double CustomRightHandSide::value(double x, double y) const
{
  double a = pow(x - x_loc, 2);
  double b = pow(y - y_loc, 2);
  double c = sqrt(a + b);
  double d = ((alpha*x - (alpha * x_loc)) * (2*x - (2 * x_loc)));
  double e = ((alpha*y - (alpha * y_loc)) * (2*y - (2 * y_loc)));
  double f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);
  double g = (alpha * c - (alpha * r_zero));
  
  return +( ((alpha/(c * f)) - (d/(2 * pow(a + b, 1.5) * f))
          - ((alpha * d * g)/((a + b) * pow(f, 2))) + (alpha/(c * f)) 
          - (e/(2 * pow(a + b, 1.5) * f))
          - ((alpha * e * g)/((a + b) * pow(f, 2)))));
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  double a = pow(x - x_loc, 2);
  double b = pow(y - y_loc, 2);
  double c = sqrt(a + b);
  double d = (alpha*x - (alpha * x_loc));
  double e = (alpha*y - (alpha * y_loc));
  double f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);
  
  dx = (d/(c * f));
  dy = (e/(c * f));
}
