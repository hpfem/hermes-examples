#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Initial condition for E */

class CustomInitialConditionWave : public ExactSolutionVector < double >
{
public:
  CustomInitialConditionWave(MeshSharedPtr mesh) : ExactSolutionVector<double>(mesh) {};

  virtual Scalar2<double> value(double x, double y) const;

  virtual void derivatives(double x, double y, Scalar2<double>& dx, Scalar2<double>& dy) const;

  virtual Ord ord(double x, double y) const;

  MeshFunction<double>* clone() const;
};

/* Weak forms */

class CustomWeakFormWave : public WeakForm < double >
{
public:

  CustomWeakFormWave(double c_squared);

private:
  class MatrixFormVolWave_0_1 : public MatrixFormVol < double >
  {
  public:
    MatrixFormVolWave_0_1(double c_squared)
      : MatrixFormVol<double>(0, 1), c_squared(c_squared) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
      Func<Real> *v, GeomVol<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;

    double c_squared;
  };

  class MatrixFormVolWave_1_0 : public MatrixFormVol < double >
  {
  public:
    MatrixFormVolWave_1_0()
      : MatrixFormVol<double>(1, 0) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
      Func<Real> *v, GeomVol<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;
  };

  class VectorFormVolWave_0 : public VectorFormVol < double >
  {
  public:
    VectorFormVolWave_0(double c_squared)
      : VectorFormVol<double>(0), c_squared(c_squared) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
      GeomVol<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord>* *ext) const;

    virtual VectorFormVol<double>* clone() const;

    double c_squared;
  };

  class VectorFormVolWave_1 : public VectorFormVol < double >
  {
  public:
    VectorFormVolWave_1()
      : VectorFormVol<double>(1) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
      GeomVol<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const;

    virtual VectorFormVol<double>* clone() const;
  };
};
