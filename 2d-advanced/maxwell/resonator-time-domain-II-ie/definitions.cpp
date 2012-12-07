#include "definitions.h"

CustomWeakFormWaveIE::CustomWeakFormWaveIE(double tau, double c_squared, Solution<double>* E_prev_sln, Solution<double>* F_prev_sln) : WeakForm<double>(2) 
{
  // Jacobian.
  add_matrix_form(new MatrixFormVolWave_0_0(tau));
  add_matrix_form(new MatrixFormVolWave_0_1);
  add_matrix_form(new MatrixFormVolWave_1_0(c_squared));
  add_matrix_form(new MatrixFormVolWave_1_1(tau));

  // Residual.
  VectorFormVolWave_0* vector_form_0 = new VectorFormVolWave_0(tau);
  vector_form_0->set_ext(E_prev_sln);
  add_vector_form(vector_form_0);
  VectorFormVolWave_1* vector_form_1 = new VectorFormVolWave_1(tau, c_squared);
  vector_form_1->set_ext(F_prev_sln);
  add_vector_form(vector_form_1);
}

double CustomWeakFormWaveIE::MatrixFormVolWave_0_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return int_e_f<double, double>(n, wt, u, v) / tau;
}

Ord CustomWeakFormWaveIE::MatrixFormVolWave_0_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return int_e_f<Ord, Ord>(n, wt, u, v) / tau;
}

MatrixFormVol<double>* CustomWeakFormWaveIE::MatrixFormVolWave_0_0::clone() const 
{
  return new CustomWeakFormWaveIE::MatrixFormVolWave_0_0(*this);
}

double CustomWeakFormWaveIE::MatrixFormVolWave_0_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return -int_e_f<double, double>(n, wt, u, v);
}

Ord CustomWeakFormWaveIE::MatrixFormVolWave_0_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return -int_e_f<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* CustomWeakFormWaveIE::MatrixFormVolWave_0_1::clone() const 
{
  return new CustomWeakFormWaveIE::MatrixFormVolWave_0_1(*this);
}

double CustomWeakFormWaveIE::MatrixFormVolWave_1_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return c_squared * int_curl_e_curl_f<double, double>(n, wt, u, v);
}

Ord CustomWeakFormWaveIE::MatrixFormVolWave_1_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return c_squared * int_curl_e_curl_f<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* CustomWeakFormWaveIE::MatrixFormVolWave_1_0::clone() const 
{
  return new CustomWeakFormWaveIE::MatrixFormVolWave_1_0(*this);
}

double CustomWeakFormWaveIE::MatrixFormVolWave_1_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  return int_e_f<double, double>(n, wt, u, v) / tau;
}

Ord CustomWeakFormWaveIE::MatrixFormVolWave_1_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return int_e_f<Ord, Ord>(n, wt, u, v) / tau;
}

MatrixFormVol<double>* CustomWeakFormWaveIE::MatrixFormVolWave_1_1::clone() const 
{
  return new CustomWeakFormWaveIE::MatrixFormVolWave_1_1(*this);
}

double CustomWeakFormWaveIE::VectorFormVolWave_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                                      Geom<double> *e, Func<double>* *ext) const 
{
  Func<double>* E_prev_newton = u_ext[0];
  Func<double>* F_prev_newton = u_ext[1];
  Func<double>* E_prev_time = ext[0];
  return   int_e_f<double, double>(n, wt, E_prev_newton, v) / tau 
         - int_e_f<double, double>(n, wt, E_prev_time, v) / tau
         - int_e_f<double, double>(n, wt, F_prev_newton, v);
}

Ord CustomWeakFormWaveIE::VectorFormVolWave_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                                 Geom<Ord> *e, Func<Ord>* *ext) const 
{
  Func<Ord>* E_prev_newton = u_ext[0];
  Func<Ord>* F_prev_newton = u_ext[1];
  Func<Ord>* E_prev_time = ext[0];
  return   int_e_f<Ord, Ord>(n, wt, E_prev_newton, v) / tau 
         - int_e_f<Ord, Ord>(n, wt, E_prev_time, v) / tau
         - int_e_f<Ord, Ord>(n, wt, F_prev_newton, v);
}

VectorFormVol<double>* CustomWeakFormWaveIE::VectorFormVolWave_0::clone() const 
{
  return new CustomWeakFormWaveIE::VectorFormVolWave_0(*this);
}

double CustomWeakFormWaveIE::VectorFormVolWave_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                                      Geom<double> *e, Func<double>* *ext) const 
{
  Func<double>* E_prev_newton = u_ext[0];
  Func<double>* F_prev_newton = u_ext[1];
  Func<double>* F_prev_time = ext[0];
  return   int_e_f<double, double>(n, wt, F_prev_newton, v) / tau 
         - int_e_f<double, double>(n, wt, F_prev_time, v) / tau
         + c_squared * int_curl_e_curl_f<double, double>(n, wt, E_prev_newton, v);
}

Ord CustomWeakFormWaveIE::VectorFormVolWave_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  Func<Ord>* E_prev_newton = u_ext[0];
  Func<Ord>* F_prev_newton = u_ext[1];
  Func<Ord>* F_prev_time = ext[0];
  return   int_e_f<Ord, Ord>(n, wt, F_prev_newton, v) / tau 
    - int_e_f<Ord, Ord>(n, wt, F_prev_time, v) / tau
    + c_squared * int_curl_e_curl_f<Ord, Ord>(n, wt, E_prev_newton, v);
}

VectorFormVol<double>* CustomWeakFormWaveIE::VectorFormVolWave_1::clone() const 
{
  return new CustomWeakFormWaveIE::VectorFormVolWave_1(*this);
}

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
  return Ord(20);
}

MeshFunction<double>* CustomInitialConditionWave::clone() const 
{
  return new CustomInitialConditionWave(this->mesh);
}


