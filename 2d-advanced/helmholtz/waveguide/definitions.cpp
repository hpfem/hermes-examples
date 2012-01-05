#include "definitions.h"

EssentialBCNonConst::EssentialBCNonConst(std::string marker) 
    : EssentialBoundaryCondition<double>(Hermes::vector<std::string>())
{
  markers.push_back(marker);
}

EssentialBoundaryCondition<double>::EssentialBCValueType EssentialBCNonConst::get_value_type() const 
{ 
  return EssentialBoundaryCondition<double>::BC_FUNCTION; 
}

double EssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
{
  return 100*std::cos(y*M_PI/0.1);
}

WeakFormHelmholtz::WeakFormHelmholtz(double eps, double mu, double omega, double sigma, double beta, 
    double E0, double h) : WeakForm<double>(2)
{
  // Jacobian.
  add_matrix_form(new MatrixFormHelmholtzEquation_real_real(0, 0, eps, omega, mu));
  add_matrix_form(new MatrixFormHelmholtzEquation_real_imag(0, 1, mu, omega, sigma));
  add_matrix_form(new MatrixFormHelmholtzEquation_imag_real(1, 0, mu, omega, sigma));
  add_matrix_form(new MatrixFormHelmholtzEquation_imag_imag(1, 1, eps, mu, omega));
  add_matrix_form_surf(new  MatrixFormSurfHelmholtz_real_imag(0, 1, "Bdy_impedance", beta));
  add_matrix_form_surf(new  MatrixFormSurfHelmholtz_imag_real(1, 0, "Bdy_impedance", beta));

  // Residual.
  add_vector_form(new VectorFormHelmholtzEquation_real(0, eps, omega, mu, sigma));
  add_vector_form(new VectorFormHelmholtzEquation_imag(1, eps, omega, mu, sigma));
  add_vector_form_surf(new VectorFormSurfHelmholtz_real(0, "Bdy_impedance", beta));
  add_vector_form_surf(new VectorFormSurfHelmholtz_imag(1, "Bdy_impedance", beta));  
}

/* Jacobian forms */

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::matrix_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) - sqr(omega) * mu * eps * int_u_v<Real, Scalar>(n, wt, u, v);
}

double WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real::clone()
{
  return new WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_real(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::matrix_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return -omega * mu * sigma * int_u_v<Real, Scalar>(n, wt, u, v);
}

double WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag::clone()
{
  return new WeakFormHelmholtz::MatrixFormHelmholtzEquation_real_imag(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::matrix_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return  omega * mu * sigma * int_u_v<Real, Scalar>(n, wt, u, v);
}

double WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real::clone()
{
  return new WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_real(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::matrix_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) - sqr(omega) * mu * eps * int_u_v<Real, Scalar>(n, wt, u, v);
}

double WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormVol<double>* WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag::clone()
{
  return new WeakFormHelmholtz::MatrixFormHelmholtzEquation_imag_imag(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::matrix_form_surf(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return beta * int_u_v<Real, Scalar>(n, wt, u, v);
}

double WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form_surf<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormSurf<double>* WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag::clone()
{
  return new WeakFormHelmholtz::MatrixFormSurfHelmholtz_real_imag(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::matrix_form_surf(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return - beta*int_u_v<Real, Scalar>(n, wt, u, v);
}

double WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form_surf<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

MatrixFormSurf<double>* WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real::clone()
{
  return new WeakFormHelmholtz::MatrixFormSurfHelmholtz_imag_real(*this);
}

/* Residual forms */

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::VectorFormHelmholtzEquation_real::vector_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return   int_grad_u_grad_v<Real, Scalar>(n, wt, u_ext[0], v) 
         - sqr(omega) * mu * eps * int_u_v<Real, Scalar>(n, wt, u_ext[0], v)
         - omega * mu * sigma * int_u_v<Real, Scalar>(n, wt, u_ext[1], v);
}

double WeakFormHelmholtz::VectorFormHelmholtzEquation_real::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord WeakFormHelmholtz::VectorFormHelmholtzEquation_real::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

VectorFormVol<double>* WeakFormHelmholtz::VectorFormHelmholtzEquation_real::clone()
{
  return new WeakFormHelmholtz::VectorFormHelmholtzEquation_real(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::VectorFormHelmholtzEquation_imag::vector_form(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return   omega * mu * sigma * int_u_v<Real, Scalar>(n, wt, u_ext[0], v) 
         - sqr(omega) * mu * eps * int_u_v<Real, Scalar>(n, wt, u_ext[1], v)
         + int_grad_u_grad_v<Real, Scalar>(n, wt, u_ext[1], v);
}

double WeakFormHelmholtz::VectorFormHelmholtzEquation_imag::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord WeakFormHelmholtz::VectorFormHelmholtzEquation_imag::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}
 
VectorFormVol<double>* WeakFormHelmholtz::VectorFormHelmholtzEquation_imag::clone()
{
  return new WeakFormHelmholtz::VectorFormHelmholtzEquation_imag(*this);
}

// Surface vector forms.

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::VectorFormSurfHelmholtz_real::vector_form_surf(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return beta * int_u_v<Real, Scalar>(n, wt, u_ext[1], v);
}

double WeakFormHelmholtz::VectorFormSurfHelmholtz_real::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form_surf<double, double>(n, wt, u_ext, v, e, ext);
}

Ord WeakFormHelmholtz::VectorFormSurfHelmholtz_real::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
}
  
VectorFormSurf<double>* WeakFormHelmholtz::VectorFormSurfHelmholtz_real::clone()
{
  return new WeakFormHelmholtz::VectorFormSurfHelmholtz_real(*this);
}

template<typename Real, typename Scalar>
Scalar WeakFormHelmholtz::VectorFormSurfHelmholtz_imag::vector_form_surf(int n, double *wt, 
    Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  return - beta*int_u_v<Real, Scalar>(n, wt, u_ext[0], v);
}

double WeakFormHelmholtz::VectorFormSurfHelmholtz_imag::value(int n, double *wt, Func<double> *u_ext[], 
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form_surf<double, double>(n, wt, u_ext, v, e, ext);
}

Ord WeakFormHelmholtz::VectorFormSurfHelmholtz_imag::ord(int n, double *wt, Func<Ord> *u_ext[], 
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
}
  
VectorFormSurf<double>* WeakFormHelmholtz::VectorFormSurfHelmholtz_imag::clone()
{
  return new WeakFormHelmholtz::VectorFormSurfHelmholtz_imag(*this);
}
