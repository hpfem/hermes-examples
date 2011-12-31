#include "hermes2d.h"

using namespace Hermes::Hermes2D;
using namespace WeakFormsH1;
using Hermes::Ord;

/* Right-hand side */

class CustomWeakForm : public WeakForm<double>
{
  class Jacobian : public MatrixFormVol<double>
  {
  public:
    Jacobian() : MatrixFormVol<double>(0, 0, Hermes::HERMES_ANY, HERMES_SYM) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
  };
  
  class Residual : public VectorFormVol<double>
  {
    const Hermes::Hermes2DFunction<double>* rhs;
  public:
    Residual(const Hermes::Hermes2DFunction<double>* rhs) : VectorFormVol<double>(0), rhs(rhs) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
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

/* Bilinear form inducing the energy norm */

class EnergyErrorForm : public Adapt<double>::MatrixFormVolError
{
public:
  EnergyErrorForm(WeakForm<double> *problem_wf) : Adapt<double>::MatrixFormVolError(0, 0, HERMES_UNSET_NORM)
  {
    this->form = problem_wf->get_mfvol()[0];
  }

  virtual double value(int n, double *wt, Func<double> *u_ext[],
                       Func<double> *u, Func<double> *v, Geom<double> *e,
                       ExtData<double> *ext) const
  {
    return this->form->value(n, wt, u_ext, u, v, e, ext);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                  Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                  ExtData<Ord> *ext) const
  {
    return this->form->ord(n, wt, u_ext, u, v, e, ext);
  }

private:
  const MatrixFormVol<double>* form;
};

/* Linear form for the residual error estimator */

class ResidualErrorForm : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
public:
  ResidualErrorForm(CustomRightHandSide* rhs) 
    : KellyTypeAdapt<double>::ErrorEstimatorForm(0), rhs(rhs) 
  { };

  double value(int n, double *wt, 
               Func<double> *u_ext[], Func<double> *u, 
               Geom<double> *e, ExtData<double> *ext) const;

  Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
          Geom<Ord> *e, ExtData<Ord> *ext) const;

private:
  CustomRightHandSide* rhs;

};

// Linear form for the interface error estimator.
class InterfaceErrorForm : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
public:
  InterfaceErrorForm() : KellyTypeAdapt<double>::ErrorEstimatorForm(0, Hermes::H2D_DG_INNER_EDGE) {};

  template<typename Real, typename Scalar>
  Real interface_estimator(int n, double *wt, 
                           Func<Scalar> *u_ext[], Func<Scalar> *u, 
                           Geom<Real> *e, ExtData<Scalar> *ext) const
  {
    Scalar result = Scalar(0);
    for (int i = 0; i < n; i++)
      result += wt[i] * Hermes::sqr(e->nx[i] * (u->get_dx_central(i) - u->get_dx_neighbor(i)) +
                                    e->ny[i] * (u->get_dy_central(i) - u->get_dy_neighbor(i)));
    return result * e->diam / 24.;
  }

  virtual double value(int n, double *wt, Func<double> *u_ext[],
                       Func<double> *u, Geom<double> *e,
                       ExtData<double> *ext) const
  {
    return interface_estimator<double, double>(n, wt, u_ext, u, e, ext);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                  Func<Ord> *u, Geom<Ord> *e,
                  ExtData<Ord> *ext) const
  {
    return interface_estimator<Ord, Ord>(n, wt, u_ext, u, e, ext);
  }
};
