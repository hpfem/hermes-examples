#include "definitions.h"

template<typename Real>
Real CustomRightHandSide::value(Real x, Real y) 
{
  Real u = Hermes::atan(this->K * x);
  Real dudx = 1. / (1 + (this->K * x) * (this->K * x)) * this->K;
  Real ddudxx = - this->K / (1 + (this->K * x) * (this->K * x)) / (1 + (this->K * x) * (this->K * x)) * 2. * this->K * this->K * x;
  return - ddudxx + u;
}


double CustomExactSolution::value(double x, double y) const 
{
  return std::atan(this->K * x);
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = 1./(1 + (this->K * x)*(this->K * x)) * this->K;
  dy = 0;
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(20);
}


CustomWeakForm::CustomWeakForm(CustomRightHandSide* rhs, std::string bdy_left_right, double K) : WeakForm<double>(1) 
{
  add_matrix_form(new CustomMatrixFormVol(0, 0));
  add_vector_form(new CustomVectorFormVol(0, rhs));
  add_vector_form_surf(new CustomVectorFormSurfRight(0, K, bdy_left_right));
  add_vector_form_surf(new CustomVectorFormSurfLeft(0, K, bdy_left_right));
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomMatrixFormVol::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                                                        Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + 
         int_u_v<Real, Scalar>(n, wt, u, v);
}

double CustomWeakForm::CustomMatrixFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                                                  Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakForm::CustomMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                             Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormVol::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                                        Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * rhs->value(e->x[i], e->y[i]) * v->val[i];
  return result;
}

double CustomWeakForm::CustomVectorFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                                  Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakForm::CustomVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                             Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormSurfRight::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                                              Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar dfdx_at_1 = Scalar(1. / (1 + (this->K * 1.) * (this->K * 1.)) * this->K);
  return - dfdx_at_1 * int_v<Real>(n, wt, v);
}

double CustomWeakForm::CustomVectorFormSurfRight::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                                        Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakForm::CustomVectorFormSurfRight::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar CustomWeakForm::CustomVectorFormSurfLeft::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                                             Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar dfdx_at_minus_1 = Scalar(-1. / (1 + (-this->K * 1.) * (-this->K * 1.)) * this->K);
  return - dfdx_at_minus_1 * int_v<Real>(n, wt, v);
}

double CustomWeakForm::CustomVectorFormSurfLeft::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                                       Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakForm::CustomVectorFormSurfLeft::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                                  Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}


