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
  public:
    Residual() : VectorFormVol<double>(0) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
                            Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;
  };
  
  public:
    CustomWeakForm()
    {
      add_matrix_form(new Jacobian);
      add_vector_form(new Residual);
    }
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:

  CustomExactSolution(Mesh* mesh, double alpha)
            : ExactSolutionScalar<double>(mesh), alpha(alpha) {};

  double value(double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord (Ord x, Ord y) const { return Ord(Ord::get_max_order()); }

  double get_angle(double y, double x) const 
  {
    double theta = atan2(y, x);
    if (theta < 0)
      theta += 2 * M_PI;
    return theta;
  };

protected:
  double alpha;
};

