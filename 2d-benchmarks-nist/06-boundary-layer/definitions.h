#include "hermes2d.h"
#include "../NIST-util.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

/* Right-hand side */

class CustomRightHandSide : public Hermes::Hermes2DFunction<double>
{
public:
  CustomRightHandSide(double epsilon)
      : Hermes::Hermes2DFunction<double>(), epsilon(epsilon) {};

  virtual double value(double x, double y) const;

  virtual Ord value (Ord x, Ord y) const;
  
  double epsilon;
};

/* Exact solution */

class CustomBC : public EssentialBoundaryCondition<double>
{
public:
  CustomBC(Hermes::vector<std::string> markers, double amplitude = 1., double frequency = 1000.)
      : EssentialBoundaryCondition<double>(markers), amplitude(amplitude), frequency(frequency)
  {
  };

  inline typename EssentialBoundaryCondition<double>::EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<double>::BC_FUNCTION; }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    return this->amplitude * std::sin(2 * M_PI * this->frequency * this->current_time);
  }

  double amplitude;
  double frequency;
};

/* Weak forms */

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(CustomRightHandSide* f);

public:
  class CustomMatrixFormVol : public MatrixFormVol<double>
  {
  public:
    CustomMatrixFormVol(int i, int j, CustomRightHandSide* f)
        : MatrixFormVol<double>(i, j), f(f) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
        Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;
    
    CustomRightHandSide* f;
  };

  class CustomVectorFormVol : public VectorFormVol<double>
  {
  public:
    CustomVectorFormVol(int i, CustomRightHandSide* f)
        : VectorFormVol<double>(i), f(f) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
        Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[],
        Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord>* *ext) const;

    VectorFormVol<double>* clone() const;

    CustomRightHandSide* f;
  };
};
