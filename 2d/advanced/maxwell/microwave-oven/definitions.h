#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Weak forms */

class CustomMatrixForm : public MatrixFormVol<std::complex<double> >
{
public:
  CustomMatrixForm(int i, int j, double e_0, double mu_0, double mu_r, double kappa, double omega, double J, bool align_mesh) 
        : MatrixFormVol<std::complex<double> >(i, j, HERMES_ANY, HERMES_SYM), e_0(e_0), mu_0(mu_0), 
          mu_r(mu_r), kappa(kappa), omega(omega), J(J), align_mesh(align_mesh) {};

  template<typename Real, typename Scalar>
  Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                     Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual std::complex<double>  value(int n, double *wt, Func<std::complex<double> > *u_ext[], Func<double> *u, 
                       Func<double> *v, Geom<double> *e, ExtData<std::complex<double> > *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                  Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

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

class CustomResidualForm : public VectorFormVol<std::complex<double> >
{
public:
  CustomResidualForm(int i, double e_0, double mu_0, double mu_r, double kappa, double omega, double J, bool align_mesh) 
        : VectorFormVol<std::complex<double> >(i), e_0(e_0), mu_0(mu_0), 
          mu_r(mu_r), kappa(kappa), omega(omega), J(J), align_mesh(align_mesh) {};

  template<typename Real, typename Scalar>
  Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                     Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual std::complex<double>  value(int n, double *wt, Func<std::complex<double> > *u_ext[], Func<double> *v, Geom<double> *e,
                       ExtData<std::complex<double> > *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                  ExtData<Ord> *ext) const;

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


class CustomVectorFormSurf : public VectorFormSurf<std::complex<double> >
{
public:
  CustomVectorFormSurf(double omega, double J) 
        : VectorFormSurf<std::complex<double> >(0), omega(omega), J(J) {};

  template<typename Scalar, typename Real>
  Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                          Geom<Real> *e, ExtData<Scalar> *ext) const;

  virtual std::complex<double>  value(int n, double *wt, Func<std::complex<double> > *u_ext[], 
                       Func<double> *v, Geom<double> *e, ExtData<std::complex<double> > *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                  Geom<Ord> *e, ExtData<Ord> *ext) const;

  double omega, J;
};

class CustomWeakForm : public WeakForm<std::complex<double> >
{
public:
  CustomWeakForm(double e_0, double mu_0, double mu_r, double kappa, double omega, 
                 double J, bool align_mesh);
};

