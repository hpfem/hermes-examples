#include "definitions.h"

using namespace WeakFormsHcurl;

Scalar2<double> CustomInitialConditionWave::value (double x, double y) const 
{
  return Scalar2<double>(std::sin(x) * std::cos(y), -std::cos(x) * std::sin(y));
}

void CustomInitialConditionWave::derivatives (double x, double y, Scalar2<double>& dx, Scalar2<double>& dy) const 
{
    dx[0] = std::cos(x) * std::cos(y);
    dx[1] = std::sin(x) * std::sin(y);
    dy[0] = -std::sin(x) * std::sin(y);
    dy[1] = -std::cos(x) * std::cos(y);
}

Ord CustomInitialConditionWave::ord(Ord x, Ord y) const
{
  return Ord(10);
}

MeshFunction<double>* CustomInitialConditionWave::clone()
{
  return new CustomInitialConditionWave(this->mesh);
}


CustomWeakFormWave::CustomWeakFormWave(double c_squared) : WeakForm<double>(2) 
{
  // Jacobian.
  add_matrix_form(new MatrixFormVolWave_0_1(c_squared));
  add_matrix_form(new MatrixFormVolWave_1_0);

  // Residual.
  add_vector_form(new VectorFormVolWave_0(c_squared));
  add_vector_form(new VectorFormVolWave_1());
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::MatrixFormVolWave_0_1::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = Scalar(0);
  for (int i=0; i < n; i++) 
  {
    result += wt[i] * (u->dy[i] * v->val0[i] - u->dx[i] * v->val1[i]);
  }
  return c_squared * result;
}

double CustomWeakFormWave::MatrixFormVolWave_0_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormWave::MatrixFormVolWave_0_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* CustomWeakFormWave::MatrixFormVolWave_0_1::clone() 
{
  return new MatrixFormVolWave_0_1(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::MatrixFormVolWave_1_0::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = Scalar(0);
  for (int i=0; i < n; i++) 
  {
    result -= wt[i] * u->curl[i] * v->val[i];
  }
  return result;
}

double CustomWeakFormWave::MatrixFormVolWave_1_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormWave::MatrixFormVolWave_1_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* CustomWeakFormWave::MatrixFormVolWave_1_0::clone() 
{
  return new MatrixFormVolWave_1_0(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::VectorFormVolWave_0::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                                            Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = Scalar(0);
  Func<Scalar>* sln_prev = u_ext[1];

  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * (sln_prev->dy[i] * v->val0[i] - sln_prev->dx[i] * v->val1[i]);
  }
  return c_squared * result;
}

double CustomWeakFormWave::VectorFormVolWave_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                                      Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormWave::VectorFormVolWave_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                                                 ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

VectorFormVol<double>* CustomWeakFormWave::VectorFormVolWave_0::clone() 
{
  return new VectorFormVolWave_0(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::VectorFormVolWave_1::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                                            Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = Scalar(0);
  Func<Scalar>* sln_prev = u_ext[0];
      
  for (int i = 0; i < n; i++) 
  {
    result -= wt[i] * sln_prev->curl[i] * v->val[i];
  }
  return result;
}

double CustomWeakFormWave::VectorFormVolWave_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                                      Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormWave::VectorFormVolWave_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

VectorFormVol<double>* CustomWeakFormWave::VectorFormVolWave_1::clone() 
{
  return new VectorFormVolWave_1(*this);
}

