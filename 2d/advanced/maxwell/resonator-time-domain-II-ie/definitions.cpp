#include "definitions.h"

CustomWeakFormWave::CustomWeakFormWave(double tau, double c_squared, Solution<double>* E_prev_sln, Solution<double>* F_prev_sln) : WeakForm(2) 
{
  add_matrix_form(new MatrixFormVolWave_0_0(tau));
  add_matrix_form(new MatrixFormVolWave_0_1);
  add_matrix_form(new MatrixFormVolWave_1_0(c_squared));
  add_matrix_form(new MatrixFormVolWave_1_1(tau));

  VectorFormVolWave_0* vector_form_0 = new VectorFormVolWave_0(tau);
  vector_form_0->ext.push_back(E_prev_sln);
  add_vector_form(vector_form_0);

  VectorFormVolWave_1* vector_form_1 = new VectorFormVolWave_1(tau);
  vector_form_1->ext.push_back(F_prev_sln);
  add_vector_form(vector_form_1);
}


template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::MatrixFormVolWave_0_0::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return int_e_f<Real, Scalar>(n, wt, u, v) / tau;
}

double CustomWeakFormWave::MatrixFormVolWave_0_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormWave::MatrixFormVolWave_0_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}


template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::MatrixFormVolWave_0_1::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return -int_e_f<Real, Scalar>(n, wt, u, v);
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

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::MatrixFormVolWave_1_0::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return c_squared * int_curl_e_curl_f<Real, Scalar>(n, wt, u, v);
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

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::MatrixFormVolWave_1_1::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                              Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return int_e_f<Real, Scalar>(n, wt, u, v) / tau;
}

double CustomWeakFormWave::MatrixFormVolWave_1_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormWave::MatrixFormVolWave_1_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}


template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::VectorFormVolWave_0::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                                            Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = Scalar(0);
  Func<Scalar>* sln_prev_time = ext->fn[0];

  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * (sln_prev_time->val0[i] * v->val0[i] + sln_prev_time->val1[i] * v->val1[i]);
  }
  return result / tau;
}

double CustomWeakFormWave::VectorFormVolWave_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                                                      Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormWave::VectorFormVolWave_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                                                 Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}


template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::VectorFormVolWave_1::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                                                            Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = Scalar(0);
  Func<Scalar>* sln_prev_time = ext->fn[0];
      
  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * (sln_prev_time->val0[i] * v->val0[i] + sln_prev_time->val1[i] * v->val1[i]);
  }
  return result / tau;
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


