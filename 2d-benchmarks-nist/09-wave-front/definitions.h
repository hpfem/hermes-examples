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
    MatrixFormVol<double>* clone();
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
    VectorFormVol<double>* clone();
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
  CustomRightHandSide(double alpha, double x_loc, double y_loc, double r_zero)
    : Hermes::Hermes2DFunction<double>(), alpha(alpha), x_loc(x_loc), y_loc(y_loc), r_zero(r_zero) 
  { };

  virtual double value(double x, double y) const;
  virtual Ord value (Ord x, Ord y) const { return Ord(8); }
  
  double alpha, x_loc, y_loc, r_zero;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(Mesh* mesh, double alpha, double x_loc, double y_loc, double r_zero)
    : ExactSolutionScalar<double>(mesh), alpha(alpha), x_loc(x_loc), y_loc(y_loc), r_zero(r_zero) 
  { };

  virtual double value(double x, double y) const {
    return atan(alpha * (sqrt(pow(x - x_loc, 2) + pow(y - y_loc, 2)) - r_zero));
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const;
  virtual Ord ord (Ord x, Ord y) const { return Ord(Ord::get_max_order()); }

  double alpha, x_loc, y_loc, r_zero;
};