#include "definitions.h"

double CustomRightHandSide::value(double x, double y) const
{
  return -epsilon*(-2*Hermes::pow(M_PI,2)*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*Hermes::cos(M_PI*(x + y))
         + 2*M_PI*(1 - exp(-(1 - x)/epsilon))*exp(-(1 - y)/epsilon)*Hermes::sin(M_PI*(x + y))/epsilon
         + 2*M_PI*(1 - exp(-(1 - y)/epsilon))*exp(-(1 - x)/epsilon)*Hermes::sin(M_PI*(x + y))/epsilon
         - (1 - exp(-(1 - y)/epsilon))*Hermes::cos(M_PI*(x + y))*exp(-(1 - x)/epsilon)/Hermes::pow(epsilon,2)
         - (1 - exp(-(1 - x)/epsilon))*Hermes::cos(M_PI*(x + y))*exp(-(1 - y)/epsilon)/Hermes::pow(epsilon,2))
         - 3*M_PI*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*Hermes::sin(M_PI*(x + y))
         - 2*(1 - exp(-(1 - y)/epsilon))*Hermes::cos(M_PI*(x + y))*exp(-(1 - x)/epsilon)/epsilon
         - (1 - exp(-(1 - x)/epsilon))*Hermes::cos(M_PI*(x + y))*exp(-(1 - y)/epsilon)/epsilon;
}

Ord CustomRightHandSide::value(Ord x, Ord y) const
{
  return Ord(8);
}

double CustomExactSolution::value(double x, double y) const 
{
  return (1 - exp(-(1-x)/epsilon)) * (1 - exp(-(1-y)/epsilon)) * Hermes::cos(M_PI * (x + y));
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = -M_PI*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*Hermes::sin(M_PI*(x + y))
       - (1 - exp(-(1 - y)/epsilon))*Hermes::cos(M_PI*(x + y))*exp(-(1 - x)/epsilon)/epsilon;
  dy = -M_PI*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*Hermes::sin(M_PI*(x + y))
       - (1 - exp(-(1 - x)/epsilon))*Hermes::cos(M_PI*(x + y))*exp(-(1 - y)/epsilon)/epsilon;
}

Ord CustomExactSolution::ord (Ord x, Ord y) const
{
  return Ord(8);
}


CustomWeakForm::CustomWeakForm(CustomRightHandSide* f) : WeakForm<double>(1) 
{
  // Jacobian.
  add_matrix_form(new CustomMatrixFormVol(0, 0, f->epsilon));
  // Residual.
  add_vector_form(new CustomVectorFormVol(0, f));
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormVol::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],     
    Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar val = Scalar(0);
  for (int i = 0; i < n; i++) 
  {
    val += wt[i] * epsilon * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    val += wt[i] * (2*u->dx[i] + u->dy[i]) * v->val[i];
  }

  return val;
}

double CustomWeakForm::CustomMatrixFormVol::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* CustomWeakForm::CustomMatrixFormVol::clone()
{
  return new CustomWeakForm::CustomMatrixFormVol(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormVol::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar val = Scalar(0);
  for (int i = 0; i < n; i++) 
  {
    val += wt[i] * f->epsilon * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);
    val += wt[i] * (2*u_ext[0]->dx[i] + u_ext[0]->dy[i]) * v->val[i];
    val -= wt[i] * f->value(e->x[i], e->y[i]) * v->val[i]; 
  }

  return val;
}

double CustomWeakForm::CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakForm::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

VectorFormVol<double>* CustomWeakForm::CustomVectorFormVol::clone()
{
  return new CustomWeakForm::CustomVectorFormVol(*this);
}
