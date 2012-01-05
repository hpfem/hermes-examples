#include "definitions.h"

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
                                                        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return int_e_f<double, double>(n, wt, u, v);
}

Ord CustomWeakFormWaveRK::MatrixFormVolWave_0_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return int_e_f<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* CustomWeakFormWaveRK::MatrixFormVolWave_0_1::clone() 
{
  return new MatrixFormVolWave_0_1(*this);
}

double CustomWeakFormWaveRK::MatrixFormVolWave_1_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return -c_squared * int_curl_e_curl_f<double, double>(n, wt, u, v);
}

Ord CustomWeakFormWaveRK::MatrixFormVolWave_1_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return -c_squared * int_curl_e_curl_f<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* CustomWeakFormWaveRK::MatrixFormVolWave_1_0::clone() 
{
  return new MatrixFormVolWave_1_0(*this);
}

double CustomWeakFormWaveRK::VectorFormVolWave_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                                      Geom<double> *e, ExtData<double> *ext) const 
{
  Func<double>* F_prev_newton = u_ext[1];
  return int_e_f<double, double>(n, wt, F_prev_newton, v);
}

Ord CustomWeakFormWaveRK::VectorFormVolWave_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                                                 ExtData<Ord> *ext) const 
{
  Func<Ord>* F_prev_newton = u_ext[1];
  return int_e_f<Ord, Ord>(n, wt, F_prev_newton, v);
}

VectorFormVol<double>* CustomWeakFormWaveRK::VectorFormVolWave_0::clone() 
{
  return new VectorFormVolWave_0(*this);
}

double CustomWeakFormWaveRK::VectorFormVolWave_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                                      Geom<double> *e, ExtData<double> *ext) const 
{
  Func<double>* E_prev_newton = u_ext[0];
  return -c_squared * int_curl_e_curl_f<double, double>(n, wt, E_prev_newton, v);
}

Ord CustomWeakFormWaveRK::VectorFormVolWave_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Func<Ord>* E_prev_newton = u_ext[0];
  return -c_squared * int_curl_e_curl_f<Ord, Ord>(n, wt, E_prev_newton, v);
}

VectorFormVol<double>* CustomWeakFormWaveRK::VectorFormVolWave_1::clone() 
{
  return new VectorFormVolWave_1(*this);
}

