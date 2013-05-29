#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

typedef std::complex<double> complex;

class CustomInitialCondition : public ExactSolutionScalar<complex>
{
public:
  CustomInitialCondition(MeshSharedPtr mesh) : ExactSolutionScalar<complex>(mesh) {};

  virtual void derivatives (double x, double y, std::complex<double> & dx, std::complex<double> & dy) const;

  virtual std::complex<double>  value (double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;

  virtual MeshFunction<complex>* clone() const;
};


/* Weak forms */

class CustomWeakFormGPRK : public WeakForm<complex>
{
public:
  CustomWeakFormGPRK(double h, double m, double g, double omega);

private:

  class CustomFormMatrixFormVol : public MatrixFormVol<complex>
  {
  public:
    CustomFormMatrixFormVol(int i, int j, double h, double m, double g, double omega) 
          : MatrixFormVol<complex>(i, j), h(h), m(m), g(g), omega(omega) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                          Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual std::complex<double>  value(int n, double *wt, Func<complex> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, Func<complex> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

    virtual MatrixFormVol<complex>* clone() const;

    // Members.
    double h, m, g, omega;
  };


  class CustomFormVectorFormVol : public VectorFormVol<complex>
  {
  public:
    CustomFormVectorFormVol(int i, double h, double m, double g, double omega)
          : VectorFormVol<complex>(i), h(h), m(m), g(g), omega(omega) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                            Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual std::complex<double>  value(int n, double *wt, Func<complex> *u_ext[], Func<double> *v, Geom<double> *e,
                         Func<complex> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    Func<Ord>* *ext) const;

    virtual VectorFormVol<complex>* clone() const;

    // Members.
    double h, m, g, omega;
  };
  double h, m, g, omega;
};

