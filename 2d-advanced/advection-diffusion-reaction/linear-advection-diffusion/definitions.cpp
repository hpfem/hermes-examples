#include "definitions.h"

WeakFormLinearAdvectionDiffusion::WeakFormLinearAdvectionDiffusion(bool stabilization_on, bool shock_capturing_on, double b1, double b2, double epsilon) 
        : WeakForm<double>(1), stabilization_on(stabilization_on), shock_capturing_on(shock_capturing_on) 
{
  add_matrix_form(new MatrixFormVolAdvectionDiffusion(0, 0, b1, b2, epsilon));
}

template<typename Real, typename Scalar>
Scalar WeakFormLinearAdvectionDiffusion::MatrixFormVolAdvectionDiffusion::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                                                      Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = Scalar(0);
  Real h_e = Real(e->diam);
  Real s_c = Real(0.9);
  for (int i=0; i < n; i++) 
  {
    result += wt[i] * (epsilon * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i])
                               - (b1 * u->val[i] * v->dx[i] + b2 * u->val[i] * v->dy[i]));
      
    if(static_cast<WeakFormLinearAdvectionDiffusion*>(wf)->shock_capturing_on) {
      Real R_squared = Hermes::pow(b1 * u->dx[i] + b2 * u->dy[i], 2.);
      Real R = Hermes::sqrt(R_squared); //This just does fabs(b1 * u->dx[i] + b2 * u->dy[i]); but it can be parsed
      result += wt[i] * s_c * 0.5 * h_e * R *
                (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) /
                (Hermes::sqrt(Hermes::pow(u->dx[i], 2) + Hermes::pow(u->dy[i], 2)) + 1.e-8);
    }
      
    if(static_cast<WeakFormLinearAdvectionDiffusion*>(wf)->stabilization_on) 
    {
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
      double b_norm = Hermes::sqrt(b1*b1 + b2*b2);
      Real tau = 1. / Hermes::sqrt(9*Hermes::pow(4*epsilon/Hermes::pow(h_e, 2), 2) + Hermes::pow(2*b_norm/h_e, 2));
      result += wt[i] * tau * (-b1 * v->dx[i] - b2 * v->dy[i] + epsilon * v->laplace[i])
                            * (-b1 * u->dx[i] - b2 * u->dy[i] + epsilon * u->laplace[i]);
#else
      error("Define H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h if you want to use second derivatives of shape functions in weak forms.");
#endif
    }
  }
  return result;
}

double WeakFormLinearAdvectionDiffusion::MatrixFormVolAdvectionDiffusion::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                                                Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormLinearAdvectionDiffusion::MatrixFormVolAdvectionDiffusion::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                                                                           Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}


EssentialBoundaryCondition<double>::EssentialBCValueType EssentialBCNonConst::get_value_type() const 
{ 
  return EssentialBoundaryCondition<double>::BC_FUNCTION; 
}

double EssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const 
{
  return 2 - std::pow(x, 0.1) - std::pow(y, 0.1);
}

