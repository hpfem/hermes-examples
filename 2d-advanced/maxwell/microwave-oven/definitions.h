#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;
using namespace Hermes::Hermes2D::WeakFormsHcurl;

typedef std::complex<double> complex;

/* Weak forms */

// Jacobian.

class CustomMatrixForm : public MatrixFormVol < ::complex >
{
public:
  CustomMatrixForm(unsigned int i, unsigned int j, double e_0, double mu_0, double mu_r, double kappa, double omega, double J, bool align_mesh)
    : MatrixFormVol<::complex>(i, j), e_0(e_0), mu_0(mu_0),
    mu_r(mu_r), kappa(kappa), omega(omega), J(J), align_mesh(align_mesh) {
    this->setSymFlag(HERMES_SYM);
  };

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
    Func<Real> *v, GeomVol<Real> *e, Func<Scalar>* *ext) const;

  virtual std::complex<double>  value(int n, double *wt, Func<::complex> *u_ext[], Func<double> *u,
    Func<double> *v, GeomVol<double> *e, Func<::complex> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
    Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const;

  // Gamma as a function of x, y.
  double gamma(int marker, double x, double y) const;

  Ord gamma(int marker, Ord x, Ord y) const;

  // Relative permittivity as a function of x, y.
  double er(int marker, double x, double y) const;

  Ord er(int marker, Ord x, Ord y) const;

  virtual MatrixFormVol<::complex>* clone() const;

  // Geometry of the load.
  bool in_load(double x, double y) const;

private:
  double e_0, mu_0, mu_r, kappa, omega, J;
  bool align_mesh;
};

// Residual.

class CustomResidualForm : public VectorFormVol < ::complex >
{
public:
  CustomResidualForm(int i, double e_0, double mu_0, double mu_r, double kappa, double omega, double J, bool align_mesh)
    : VectorFormVol<::complex>(i), e_0(e_0), mu_0(mu_0),
    mu_r(mu_r), kappa(kappa), omega(omega), J(J), align_mesh(align_mesh) {};

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
    GeomVol<Real> *e, Func<Scalar>* *ext) const;

  virtual std::complex<double>  value(int n, double *wt, Func<::complex> *u_ext[], Func<double> *v, GeomVol<double> *e,
    Func<::complex> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e,
    Func<Ord>* *ext) const;

  virtual VectorFormVol<::complex>* clone() const;

  // Gamma as a function of x, y.
  double gamma(int marker, double x, double y) const;

  Ord gamma(int marker, Ord x, Ord y) const;

  // Relative permittivity as a function of x, y.
  double er(int marker, double x, double y) const;

  Ord er(int marker, Ord x, Ord y) const;

  // Geometry of the load.
  bool in_load(double x, double y) const;

private:
  double e_0, mu_0, mu_r, kappa, omega, J;
  bool align_mesh;
};

class CustomVectorFormSurf : public VectorFormSurf < ::complex >
{
public:
  CustomVectorFormSurf(double omega, double J, std::string bnd)
    : VectorFormSurf<::complex>(0), omega(omega), J(J) { this->set_area(bnd); };

  template<typename Scalar, typename Real>
  Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
    GeomSurf<Real> *e, Func<Scalar>* *ext) const;

  virtual std::complex<double>  value(int n, double *wt, Func<::complex> *u_ext[],
    Func<double> *v, GeomSurf<double> *e, Func<::complex> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    GeomSurf<Ord> *e, Func<Ord>* *ext) const;

  virtual VectorFormSurf<::complex>* clone() const { return new CustomVectorFormSurf(omega, J, this->areas[0]); }
  double omega, J;
};

class CustomWeakForm : public WeakForm < ::complex >
{
public:
  CustomWeakForm(double e_0, double mu_0, double mu_r, double kappa, double omega,
    double J, bool align_mesh, MeshSharedPtr mesh, std::string current_bdy);
  int get_marker();

private:
  int marker;
};

// Custom error form.

class CustomErrorForm : public NormFormVol < ::complex >
{
public:
  CustomErrorForm(double kappa) : NormFormVol<::complex>(0, 0)
  {
    kappa_squared = sqr(kappa);
  };

  template<typename Real, typename Scalar>
  Scalar hcurl_form_kappa(int n, double *wt, Func<Scalar> *u, Func<Scalar> *v, GeomVol<Real> *e) const
  {
    Scalar result = Scalar(0);

    for (int i = 0; i < n; i++)
    {
      result += wt[i] * (u->curl[i] * conj(v->curl[i]));
      result += kappa_squared * wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
    }

    return result;
  }

  virtual std::complex<double> value(int n, double *wt, Func<::complex> *u, Func<::complex> *v, GeomVol<double> *e) const
  {
    return hcurl_form_kappa<double, std::complex<double> >(n, wt, u, v, e);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e) const
  {
    return hcurl_form_kappa<Ord, Ord>(n, wt, u, v, e);
  }

  double kappa_squared;
};
