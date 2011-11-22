#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Weak forms */

class CustomWeakFormHeatMoistureRK : public WeakForm<double>
{
public:
  CustomWeakFormHeatMoistureRK(double c_TT, double c_ww, double d_TT, double d_Tw, double d_wT, double d_ww, 
      double k_TT, double k_ww, double T_ext, double w_ext, const std::string bdy_ext);
};

/* Time-dependent Dirichlet condition for temperature */

class EssentialBCNonConst : public EssentialBoundaryCondition<double>
{
public:
  EssentialBCNonConst(std::string marker, double reactor_start_time, double temp_initial, double temp_reactor_max);
  
  ~EssentialBCNonConst() {};

  virtual EssentialBoundaryCondition<double>::EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

protected:
  double reactor_start_time;
  double temp_initial;
  double temp_reactor_max;
};

/* Custom error forms */

class CustomErrorForm : public Adapt<double>::MatrixFormVolError
{
public:
  CustomErrorForm(double d, double c) : Adapt<double>::MatrixFormVolError(), d(d), c(c) {};

  template<typename Real, typename Scalar>
  Scalar laplace_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u, 
      Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
  {
    return d / c * int_grad_u_grad_v<Scalar, Scalar>(n, wt, u, v);
  }

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const
  {
    return laplace_form<double, double>(n, wt, u_ext, u, v, e, ext);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    return laplace_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  double d, c;
};

