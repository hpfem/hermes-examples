#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Initial condition */

class CustomInitialConditionWave : public ExactSolutionScalar<double>
{
public:
  CustomInitialConditionWave(Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
    
  virtual MeshFunction<double>* clone();
};

/* Weak forms */

class CustomWeakFormWave : public WeakForm<double>
{
public:

  CustomWeakFormWave(double tau, double c_squared, Solution<double>* u_prev_sln, Solution<double>* v_prev_sln);

private:
  class VectorFormVolWave_0 : public VectorFormVol<double>
  {
  public:
    VectorFormVolWave_0() : VectorFormVol<double>(0) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                       Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const;

    virtual VectorFormVol<double>* clone();
  };

  class VectorFormVolWave_1 : public VectorFormVol<double>
  {
  public:
    VectorFormVolWave_1(double c_squared)
          : VectorFormVol<double>(1), c_squared(c_squared) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                       Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const;

    virtual VectorFormVol<double>* clone();

    double c_squared;
  };
};
