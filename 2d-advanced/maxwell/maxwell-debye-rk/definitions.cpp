#include "definitions.h"

/* Global function alpha */

double alpha(double omega, double k)
{
  return Hermes::sqr(omega) - omega + Hermes::sqr(k) * Hermes::sqr(M_PI);
}

// **************
Scalar2<double> CustomInitialConditionE::value (double x, double y) const 
{
  double val_0 = -omega*k_y*exp(-omega*time) * std::cos(k_x * M_PI * x) * std::sin(k_y * M_PI * y);
  double val_1 = omega*k_x*exp(-omega*time) * std::sin(k_x * M_PI * x) * std::cos(k_y * M_PI * y);
  return Scalar2<double>(val_0, val_1);
}

void CustomInitialConditionE::derivatives (double x, double y, Scalar2<double>& dx, Scalar2<double>& dy) const 
{
  dx[0] = -omega*k_y*exp(-omega*time) * (-1.0*std::sin(k_x * M_PI * x)*k_x*M_PI) * std::sin(k_y * M_PI * y);
  dx[1] = omega*k_x*exp(-omega*time) * std::cos(k_x * M_PI * x)*k_x*M_PI * std::cos(k_y * M_PI * y);
  dy[0] = -omega*k_y*exp(-omega*time) * std::cos(k_x * M_PI * x) * std::cos(k_y * M_PI * y)*k_y*M_PI;
  dy[1] = omega*k_x*exp(-omega*time) * std::sin(k_x * M_PI * x) * (-1.0*std::sin(k_y * M_PI * y)*k_y*M_PI);
}

Ord CustomInitialConditionE::ord(Ord x, Ord y) const 
{
  return Ord(20);
}

// **************
double CustomInitialConditionH::value (double x, double y) const 
{
  double k_squared = Hermes::sqr(k_x) + Hermes::sqr(k_y);
  return k_squared * M_PI * exp(-omega*time) * std::cos(k_x * M_PI * x) * std::cos(k_y * M_PI * y);
}

void CustomInitialConditionH::derivatives (double x, double y, double& dx, double& dy) const 
{
  double k_squared = Hermes::sqr(k_x) + Hermes::sqr(k_y);
  dx = k_squared * M_PI * exp(-omega*time) * (-1.0*std::sin(k_x * M_PI * x)*k_x*M_PI) * std::cos(k_y * M_PI * y);
  dy = k_squared * M_PI * exp(-omega*time) * std::cos(k_x * M_PI * x) * (-1.0*std::sin(k_y * M_PI * y)*k_y*M_PI);
}

Ord CustomInitialConditionH::ord(Ord x, Ord y) const 
{
  return Ord(20);
}

// **************
Scalar2<double> CustomInitialConditionP::value (double x, double y) const 
{
  double k_squared = Hermes::sqr(k_x) + Hermes::sqr(k_y);
  double k = std::sqrt(k_squared);
  double val_0 = alpha(omega, k)*k_y*exp(-omega*time) * std::cos(k_x * M_PI * x) * std::sin(k_y * M_PI * y);
  double val_1 = -k_x*alpha(omega, k)*exp(-omega*time) * std::sin(k_x * M_PI * x) * std::cos(k_y * M_PI * y);
  return Scalar2<double>(val_0, val_1);
}

void CustomInitialConditionP::derivatives (double x, double y, Scalar2<double>& dx, Scalar2<double>& dy) const 
{
  double k_squared = std::abs(Hermes::sqr(k_x) + Hermes::sqr(k_y));
  double k = std::sqrt(k_squared);
  dx[0] = alpha(omega, k)*k_y*exp(-omega*time) * (-1.0*std::sin(k_x * M_PI * x)*k_x*M_PI) * std::sin(k_y * M_PI * y);
  dx[1] = -k_x*alpha(omega, k)*exp(-omega*time) * std::cos(k_x * M_PI * x)*k_x*M_PI * std::cos(k_y * M_PI * y);
  dy[0] = alpha(omega, k)*k_y*exp(-omega*time) * std::cos(k_x * M_PI * x) * std::cos(k_y * M_PI * y)*k_y*M_PI;
  dy[1] = -k_x*alpha(omega, k)*exp(-omega*time) * std::sin(k_x * M_PI * x) * (-1.0*std::sin(k_y * M_PI * y)*k_y*M_PI);
}

Ord CustomInitialConditionP::ord(Ord x, Ord y) const 
{
  return Ord(20);
}

// **************
CustomWeakFormMD::CustomWeakFormMD(double omega, double k_x, double k_y, double mu_0, 
    double eps_0, double eps_inf, double eps_q, double tau) : WeakForm<double>(3) 
{
  // Stationary Jacobian.
  add_matrix_form(new MatrixFormVolMD_0_0(eps_q, tau));
  add_matrix_form(new MatrixFormVolMD_0_1(eps_0, eps_inf));
  add_matrix_form(new MatrixFormVolMD_0_2(eps_0, eps_inf, tau));
  add_matrix_form(new MatrixFormVolMD_1_0(mu_0));
  add_matrix_form(new MatrixFormVolMD_2_0(eps_0, eps_inf, eps_q, tau));
  add_matrix_form(new MatrixFormVolMD_2_2(tau));

  // Stationary Residual.
  add_vector_form(new VectorFormVolMD_0(eps_0, eps_inf, eps_q, tau));
  add_vector_form(new VectorFormVolMD_1(mu_0));
  add_vector_form(new VectorFormVolMD_2(eps_0, eps_inf, eps_q, tau));
}

// **************
double CustomWeakFormMD::MatrixFormVolMD_0_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return -(eps_q - 1)/tau * int_e_f<double, double>(n, wt, u, v);
}

Ord CustomWeakFormMD::MatrixFormVolMD_0_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
     Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return -(eps_q - 1)/tau * int_e_f<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* CustomWeakFormMD::MatrixFormVolMD_0_0::clone() 
{
  return new MatrixFormVolMD_0_0(*this);
}

// **************
double CustomWeakFormMD::MatrixFormVolMD_0_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  double result = 0;
  // int \curl u \dot v
  for (int i=0; i < n; i++)
  {
    result += wt[i] * (-u->dy[i] * v->val0[i] + u->dx[i] * v->val1[i]);   
  }
  return result * (1.0 / eps_0 / eps_inf);
}

Ord CustomWeakFormMD::MatrixFormVolMD_0_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                   Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Ord result = Ord(0);
  for (int i=0; i < n; i++)
  {
    result += wt[i] * (-u->dy[i] * v->val0[i] + u->dx[i] * v->val1[i]);   
  }
  return result * (1.0 / eps_0 / eps_inf);
}

MatrixFormVol<double>* CustomWeakFormMD::MatrixFormVolMD_0_1::clone() 
{
  return new MatrixFormVolMD_0_1(*this);
}

// **************
double CustomWeakFormMD::MatrixFormVolMD_0_2::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return (1.0 / eps_0 / eps_inf / tau) * int_e_f<double, double>(n, wt, u, v);
}

Ord CustomWeakFormMD::MatrixFormVolMD_0_2::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
    Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return (1.0 / eps_0 / eps_inf / tau) * int_e_f<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* CustomWeakFormMD::MatrixFormVolMD_0_2::clone() 
{
  return new MatrixFormVolMD_0_2(*this);
}

// **************
double CustomWeakFormMD::MatrixFormVolMD_1_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  double result = 0;
  // int \curl u \dot v
  for (int i=0; i < n; i++)
  {
    result += wt[i] * (u->curl[i] * v->val[i]);   
  }
  return result * (-1.0 / mu_0);
}

Ord CustomWeakFormMD::MatrixFormVolMD_1_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
    Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  Ord result = Ord(0);
  // int \curl u \dot v
  for (int i=0; i < n; i++)
  {
    result += wt[i] * (u->curl[i] * v->val[i]);   
  }
  return result * (-1.0 / mu_0); 
}

MatrixFormVol<double>* CustomWeakFormMD::MatrixFormVolMD_1_0::clone() 
{
  return new MatrixFormVolMD_1_0(*this);
}

// **************
double CustomWeakFormMD::MatrixFormVolMD_2_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return (eps_q - 1)*eps_0*eps_inf/tau * int_e_f<double, double>(n, wt, u, v);
}

Ord CustomWeakFormMD::MatrixFormVolMD_2_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
     Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return (eps_q - 1)*eps_0*eps_inf/tau * int_e_f<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* CustomWeakFormMD::MatrixFormVolMD_2_0::clone() 
{
  return new MatrixFormVolMD_2_0(*this);
}

// **************
double CustomWeakFormMD::MatrixFormVolMD_2_2::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return (-1.0 / tau) * int_e_f<double, double>(n, wt, u, v);
}

Ord CustomWeakFormMD::MatrixFormVolMD_2_2::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
    Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return (-1.0 / tau) * int_e_f<Ord, Ord>(n, wt, u, v);
}

MatrixFormVol<double>* CustomWeakFormMD::MatrixFormVolMD_2_2::clone() 
{
  return new MatrixFormVolMD_2_2(*this);
}

// **************
double CustomWeakFormMD::VectorFormVolMD_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
    Geom<double> *e, ExtData<double> *ext) const 
{
  Func<double>* E_prev_newton = u_ext[0];
  Func<double>* H_prev_newton = u_ext[1];
  Func<double>* P_prev_newton = u_ext[2];
  double int_curl_H_v = 0;
  for (int i=0; i < n; i++)
  {
    int_curl_H_v += wt[i] * (-H_prev_newton->dy[i] * v->val0[i] + H_prev_newton->dx[i] * v->val1[i]);   
  }

  return -(eps_q - 1)/tau * int_e_f<double, double>(n, wt, E_prev_newton, v)  
         + (1.0 / eps_0 / eps_inf) * int_curl_H_v
         + (1.0 / eps_0 / eps_inf / tau) * int_e_f<double, double>(n, wt, P_prev_newton, v);
}

Ord CustomWeakFormMD::VectorFormVolMD_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                                                 ExtData<Ord> *ext) const 
{
  Func<Ord>* E_prev_newton = u_ext[0];
  Func<Ord>* H_prev_newton = u_ext[1];
  Func<Ord>* P_prev_newton = u_ext[2];
  Ord int_curl_H_v = Ord(0);
  for (int i=0; i < n; i++)
  {
    int_curl_H_v += wt[i] * (-H_prev_newton->dy[i] * v->val0[i] + H_prev_newton->dx[i] * v->val1[i]);   
  }

  return -(eps_q - 1)/tau * int_e_f<Ord, Ord>(n, wt, E_prev_newton, v)  
         + (1.0 / eps_0 / eps_inf) * int_curl_H_v
         + (1.0 / eps_0 / eps_inf / tau) * int_e_f<Ord, Ord>(n, wt, P_prev_newton, v);
}

VectorFormVol<double>* CustomWeakFormMD::VectorFormVolMD_0::clone() 
{
  return new VectorFormVolMD_0(*this);
}

// **************
double CustomWeakFormMD::VectorFormVolMD_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
    Geom<double> *e, ExtData<double> *ext) const 
{
  Func<double>* E_prev_newton = u_ext[0];
  double int_curl_E_v = 0;
  for (int i=0; i < n; i++)
  {
    int_curl_E_v += wt[i] * (E_prev_newton->curl[i] * v->val[i]);
  }
  return (-1.0 / mu_0) * int_curl_E_v;
}

Ord CustomWeakFormMD::VectorFormVolMD_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                                                 ExtData<Ord> *ext) const 
{
  Func<Ord>* E_prev_newton = u_ext[0];
  Ord int_curl_E_v = Ord(0);
  for (int i=0; i < n; i++)
  {
    int_curl_E_v += wt[i] * (E_prev_newton->curl[i] * v->val[i]);
  }
  return (-1.0 / mu_0) * int_curl_E_v;
}

VectorFormVol<double>* CustomWeakFormMD::VectorFormVolMD_1::clone() 
{
  return new VectorFormVolMD_1(*this);
}

// **************
double CustomWeakFormMD::VectorFormVolMD_2::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
    Geom<double> *e, ExtData<double> *ext) const 
{
  Func<double>* E_prev_newton = u_ext[0];
  Func<double>* P_prev_newton = u_ext[2];
  return (eps_q - 1)*eps_0*eps_inf/tau * int_e_f<double, double>(n, wt, E_prev_newton, v)  
         - (1.0 / tau) * int_e_f<double, double>(n, wt, P_prev_newton, v);
}

Ord CustomWeakFormMD::VectorFormVolMD_2::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                                                 ExtData<Ord> *ext) const 
{
  Func<Ord>* E_prev_newton = u_ext[0];
  Func<Ord>* P_prev_newton = u_ext[2];
  return (eps_q - 1)*eps_0*eps_inf/tau * int_e_f<Ord, Ord>(n, wt, E_prev_newton, v)  
         - (1.0 / tau) * int_e_f<Ord, Ord>(n, wt, P_prev_newton, v);
}

VectorFormVol<double>* CustomWeakFormMD::VectorFormVolMD_2::clone() 
{
  return new VectorFormVolMD_2(*this);
}

