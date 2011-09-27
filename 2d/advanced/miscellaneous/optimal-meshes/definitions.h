#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Right-hand side for the 2D equation -Laplace u = f with Dirichlet BC.
class CustomRightHandSide
{
public:
  CustomRightHandSide(double K) : K(K) 
  {
  };

  template<typename Real>
  Real value(Real x, Real y);

  // Member.
  double K;
};

// Exact solution (needed in the Dirichlet condition).
class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(Mesh* mesh, double K) : ExactSolutionScalar<double>(mesh), K(K) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  // Members.
  double K;
};

// Weak forms.
class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(CustomRightHandSide* rhs, std::string bdy_left_right, double K);

private:
  class CustomMatrixFormVol : public MatrixFormVol<double>
  {
  public:
    CustomMatrixFormVol(int i, int j) : MatrixFormVol<double>(i, j) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
  };

  class CustomVectorFormVol : public VectorFormVol<double>
  {
  public:
    CustomVectorFormVol(int i, CustomRightHandSide* rhs) : VectorFormVol<double>(i), rhs(rhs) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Members.
    CustomRightHandSide* rhs;
  };

  class CustomVectorFormSurfRight : public VectorFormSurf<double>
  {
  public:
    CustomVectorFormSurfRight(int i, double K, std::string area) : VectorFormSurf<double>(i, area), K(K) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Members.
    double K;
  };

  class CustomVectorFormSurfLeft : public VectorFormSurf<double>
  {
  public:
    CustomVectorFormSurfLeft(int i, double K, std::string area) : VectorFormSurf(i, area), K(K) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    // Members.
    double K;
  };
};
