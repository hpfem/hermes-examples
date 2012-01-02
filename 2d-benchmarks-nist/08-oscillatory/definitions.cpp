#include "definitions.h"

double CustomRightHandSide::value(double x, double y) const
{
  return (-Hermes::sin(1.0/(alpha + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))))/Hermes::pow((alpha
           + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))),4)
           + 2*Hermes::cos(1.0/(alpha + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))))/(Hermes::pow((alpha
           + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))),2)*Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0)))
           + Hermes::pow(x,2)*Hermes::sin(1.0/(alpha + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))))/(Hermes::pow((alpha
           + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))),4)*(Hermes::pow(x,2) + Hermes::pow(y,2)))
           + Hermes::pow(y,2)*Hermes::sin(1.0/(alpha + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))))/(Hermes::pow((alpha
           + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))),4)*(Hermes::pow(x,2) + Hermes::pow(y,2)))
           - Hermes::pow(x,2)*Hermes::cos(1.0/(alpha + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))))/(Hermes::pow((alpha
           + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))),2)*Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(3.0/2.0)))
           - Hermes::pow(y,2)*Hermes::cos(1.0/(alpha + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))))/(Hermes::pow((alpha
           + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))),2)*Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(3.0/2.0)))
           - 2*Hermes::pow(x,2)*Hermes::cos(1.0/(alpha + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))))/(Hermes::pow((alpha
           + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))),3)*(Hermes::pow(x,2) + Hermes::pow(y,2)))
           - 2*Hermes::pow(y,2)*Hermes::cos(1.0/(alpha + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))))/(Hermes::pow((alpha
           + Hermes::pow((Hermes::pow(x,2) + Hermes::pow(y,2)),(1.0/2.0))),3)*(Hermes::pow(x,2) + Hermes::pow(y,2))));
}

Ord CustomRightHandSide::value(Ord x, Ord y) const
{
  return Ord(10);
}

double CustomExactSolution::value(double x, double y) const 
{
  double r = Hermes::sqrt(x*x + y*y);
  return Hermes::sin(1/(alpha + r));
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  double r = Hermes::sqrt(x*x + y*y);
  double h = 1/(alpha + r);
  dx = -Hermes::cos(h) * h * h * x / r;
  dy = -Hermes::cos(h) * h * h * y / r;
}

Ord CustomExactSolution::ord (Ord x, Ord y) const
{
  return Ord(10);
}


CustomWeakForm::CustomWeakForm(CustomRightHandSide* f) : WeakForm<double>(1) 
{
  // Jacobian.
  add_matrix_form(new CustomMatrixFormVol(0, 0, f->alpha));
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
    Scalar x = e->x[i];
    Scalar y = e->y[i];
    Scalar r = Hermes::sqrt(x*x + y*y);
    Scalar h = 1/(alpha + r);
    Scalar grad_u_grad_v = u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i];
    val += wt[i] * (grad_u_grad_v - Hermes::pow(h, 4) * u->val[i] * v->val[i]);
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
    Scalar x = e->x[i];
    Scalar y = e->y[i];
    Scalar r = Hermes::sqrt(x*x + y*y);
    Scalar h = 1/(f->alpha + r);
    Scalar grad_u_grad_v = u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i];
    val += wt[i] * (grad_u_grad_v - Hermes::pow(h, 4) * u_ext[0]->val[i] * v->val[i]);
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
