#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;


class CustomInitialCondition : public ExactSolutionScalar<std::complex<double> >
{
public:
  CustomInitialCondition(Mesh* mesh) : ExactSolutionScalar<std::complex<double> >(mesh) {};

  virtual void derivatives (double x, double y, std::complex<double> & dx, std::complex<double> & dy) const;

  virtual std::complex<double>  value (double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;

  virtual MeshFunction<std::complex<double> >* clone();
};


/* Weak forms */

class CustomWeakFormGPRK : public WeakForm<std::complex<double> >
{
public:
  CustomWeakFormGPRK(double h, double m, double g, double omega);

private:

  class CustomFormMatrixFormVol : public MatrixFormVol<std::complex<double> >
  {
  public:
    CustomFormMatrixFormVol(int i, int j, double h, double m, double g, double omega) 
          : MatrixFormVol<std::complex<double> >(i, j), h(h), m(m), g(g), omega(omega) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                          Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual std::complex<double>  value(int n, double *wt, Func<std::complex<double> > *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<std::complex<double> > *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<std::complex<double> >* clone();

    // Members.
    double h, m, g, omega;
  };


  class CustomFormVectorFormVol : public VectorFormVol<std::complex<double> >
  {
  public:
    CustomFormVectorFormVol(int i, double h, double m, double g, double omega)
          : VectorFormVol<std::complex<double> >(i), h(h), m(m), g(g), omega(omega) {};

    template<typename Real, typename Scalar>
    Scalar vector_form_rk(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                            Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual std::complex<double>  value(int n, double *wt, Func<std::complex<double> > *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<std::complex<double> > *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    virtual VectorFormVol<std::complex<double> >* clone();

    // Members.
    double h, m, g, omega;
  };
  double h, m, g, omega;
};

