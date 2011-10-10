#include "hermes2d.h"

using namespace Hermes::Hermes2D;
using namespace WeakFormsH1;
using Hermes::Ord;

/* Alternatively, DefaultWeakFormDiffusion may be used. This is just copied from the Kelly version of this benchmark
   for one-to-one comparison.
*/
class CustomWeakForm : public WeakForm<double>
{
  class Jacobian : public MatrixFormVol<double>
  {
  public:
    Jacobian() : MatrixFormVol<double>(0, 0, Hermes::HERMES_ANY, HERMES_SYM) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;
  };
  
  class Residual : public VectorFormVol<double>
  {
    const Hermes::Hermes2DFunction<double>* rhs;
  public:
    Residual(const Hermes::Hermes2DFunction<double>* rhs) : VectorFormVol<double>(0), rhs(rhs) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<double> *ext) const;

    virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
        Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;
  };
  
  public:
    CustomWeakForm(const Hermes::Hermes2DFunction<double>* rhs)
    {
      add_matrix_form(new Jacobian);
      add_vector_form(new Residual(rhs));
    }
};

class CustomRightHandSide : public Hermes::Hermes2DFunction<double>
{
public:
  CustomRightHandSide(double poly_deg)
      : Hermes::Hermes2DFunction<double>(), poly_deg(poly_deg) {};

  virtual double value(double x, double y) const;
  virtual Ord value (Ord x, Ord y) const { return Ord(8); }
  
  double poly_deg;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(Mesh* mesh, double poly_deg)
      : ExactSolutionScalar<double>(mesh), poly_deg(poly_deg) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord (Ord x, Ord y) const { return Ord(Ord::get_max_order()); }

  double poly_deg;
};



