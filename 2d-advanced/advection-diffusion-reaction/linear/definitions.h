#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

class WeakFormLinearAdvectionDiffusion : public WeakForm <double>
{
public:
  // Problem parameters.
  double const_f;

  WeakFormLinearAdvectionDiffusion(bool stabilization_on, bool shock_capturing_on, double b1, double b2, double epsilon);

private:
  class MatrixFormVolAdvectionDiffusion : public MatrixFormVol <double>
  {
  public:
    MatrixFormVolAdvectionDiffusion(int i, int j, double b1, double b2, double epsilon)
      : MatrixFormVol<double>(i, j), b1(b1), b2(b2), epsilon(epsilon) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, GeomVol<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;

    // Members.
    double b1, b2, epsilon;
  };

  // Members.
  bool stabilization_on;
  bool shock_capturing_on;
};

/* Essential BC */

class EssentialBCNonConst : public EssentialBoundaryCondition <double>
{
public:
  EssentialBCNonConst(std::string marker) : EssentialBoundaryCondition<double>(marker) {};

  ~EssentialBCNonConst() {};

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y) const;
};
