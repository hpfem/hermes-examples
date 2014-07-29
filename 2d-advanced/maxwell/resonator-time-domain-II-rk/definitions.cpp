#include "definitions.h"

template<typename Real, typename Scalar>
static Scalar int_e_f(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
  return result;
}

template<typename Real, typename Scalar>
static Scalar int_curl_e_curl_f(int n, double *wt, Func<Real> *u, Func<Real> *v)
{
  Scalar result = Scalar(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->curl[i] * conj(v->curl[i]));
  return result;
}

Scalar2<double> CustomInitialConditionWave::value(double x, double y) const
{
  return Scalar2<double>(std::sin(x) * std::cos(y), -std::cos(x) * std::sin(y));
}

void CustomInitialConditionWave::derivatives(double x, double y, Scalar2<double>& dx, Scalar2<double>& dy) const
{
  dx[0] = std::cos(x) * std::cos(y);
  dx[1] = std::sin(x) * std::sin(y);
  dy[0] = -std::sin(x) * std::sin(y);
  dy[1] = -std::cos(x) * std::cos(y);
}

Ord CustomInitialConditionWave::ord(double x, double y) const
{
  return Ord(10);
}

MeshFunction<double>* CustomInitialConditionWave::clone() const
{
  return new CustomInitialConditionWave(this->mesh);
}

CustomWeakFormWaveRK::CustomWeakFormWaveRK(double c_squared) : WeakForm<double>(2)
{
  // Stationary Jacobian.
  add_matrix_form(new MatrixFormVolWave_0_1);
  add_matrix_form(new MatrixFormVolWave_1_0(c_squared));

  // Stationary Residual.
  add_vector_form(new VectorFormVolWave_0());
  add_vector_form(new VectorFormVolWave_1(c_squared));
}

double CustomWeakFormWaveRK::MatrixFormVolWave_0_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
  Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const
{
  return int_e_f<double, double>(n, wt, u, v);
}

Ord CustomWeakFormWaveRK::MatrixFormVolWave_0_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
  GeomVol<Ord> *e, Func<Ord>* *ext) const
{
  return int_e_f<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* CustomWeakFormWaveRK::MatrixFormVolWave_0_1::clone() const
{
  return new MatrixFormVolWave_0_1(*this);
}

double CustomWeakFormWaveRK::MatrixFormVolWave_1_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
  Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const
{
  return -c_squared * int_curl_e_curl_f<double, double>(n, wt, u, v);
}

Ord CustomWeakFormWaveRK::MatrixFormVolWave_1_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
  GeomVol<Ord> *e, Func<Ord>* *ext) const
{
  return -c_squared * int_curl_e_curl_f<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* CustomWeakFormWaveRK::MatrixFormVolWave_1_0::clone() const
{
  return new MatrixFormVolWave_1_0(*this);
}

double CustomWeakFormWaveRK::VectorFormVolWave_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
  GeomVol<double> *e, Func<double>* *ext) const
{
  Func<double>* F_prev_newton = u_ext[1];
  return int_e_f<double, double>(n, wt, F_prev_newton, v);
}

Ord CustomWeakFormWaveRK::VectorFormVolWave_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e,
  Func<Ord>* *ext) const
{
  Func<Ord>* F_prev_newton = u_ext[1];
  return int_e_f<Ord, Ord>(n, wt, F_prev_newton, v);
}

VectorFormVol<double>* CustomWeakFormWaveRK::VectorFormVolWave_0::clone() const
{
  return new VectorFormVolWave_0(*this);
}

double CustomWeakFormWaveRK::VectorFormVolWave_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
  GeomVol<double> *e, Func<double>* *ext) const
{
  Func<double>* E_prev_newton = u_ext[0];
  return -c_squared * int_curl_e_curl_f<double, double>(n, wt, E_prev_newton, v);
}

Ord CustomWeakFormWaveRK::VectorFormVolWave_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const
{
  Func<Ord>* E_prev_newton = u_ext[0];
  return -c_squared * int_curl_e_curl_f<Ord, Ord>(n, wt, E_prev_newton, v);
}

VectorFormVol<double>* CustomWeakFormWaveRK::VectorFormVolWave_1::clone() const
{
  return new VectorFormVolWave_1(*this);
}