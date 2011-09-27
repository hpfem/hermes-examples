#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Element material markers and associated diffusion coefficients.
extern const std::string MATERIAL_1;
extern const double EPS_1;
extern const std::string MATERIAL_2;
extern const double EPS_2;

// Boundary markers and boundary condition data.
extern const std::string BDY_BOTTOM, BDY_OUTER, BDY_LEFT, BDY_INNER;
extern const double T1;                   // Prescribed temperature on Gamma_left.
extern const double T0;                   // Outer temperature on Gamma_bottom.
extern const double H;                    // Heat flux on Gamma_bottom.
extern const double CONST_GAMMA_OUTER;    // Heat flux on Gamma_inner.

// Right-hand side.
extern const double CONST_F;

class WeakFormPoisson : public WeakForm<double>
{
public:
  WeakFormPoisson(double EPS_1, double EPS_2, double H, double T0, std::string material_1, std::string material_2, std::string bdy_bottom, std::string bdy_outer, double gamma_outer) : WeakForm<double>(1) {
    add_matrix_form(new BilinearForm(0, 0, EPS_1, material_1));
    add_matrix_form(new BilinearForm(0, 0, EPS_2, material_2));

    add_matrix_form_surf(new BilinearFormSurfBottom(0, 0, H, bdy_bottom));
    add_vector_form_surf(new LinearFormSurfBottom(0, H, T0, bdy_bottom));
    add_vector_form_surf(new LinearFormSurfOuter(0, gamma_outer, bdy_outer));
  };

  class BilinearForm : public MatrixFormVol<double>
  {
  public:
    BilinearForm(int i, int j, double EPS, std::string area) : MatrixFormVol<double>(i, j, area, HERMES_SYM), EPS(EPS) {}

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return EPS * int_grad_u_grad_v<double, double>(n, wt, u, v);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return EPS * int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
    }
  protected:
    // Members.
    double EPS;
  };

  class LinearFormSurfOuter : public VectorFormSurf<double>
  {
  public:
    LinearFormSurfOuter(int i, double CONST_GAMMA_OUTER, std::string area) : VectorFormSurf<double>(i, area), CONST_GAMMA_OUTER(CONST_GAMMA_OUTER) {}

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return CONST_GAMMA_OUTER * int_v<double>(n, wt, v);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return CONST_GAMMA_OUTER * int_v<Ord>(n, wt, v);
    }
  protected:
    // Members.
    double CONST_GAMMA_OUTER;
  };

  class BilinearFormSurfBottom : public MatrixFormSurf<double>
  {
  public:
    BilinearFormSurfBottom(int i, int j, double H, std::string area) : MatrixFormSurf<double>(i, j, area), H(H) {}

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return H * int_u_v<double, double>(n, wt, u, v);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return H * int_u_v<Ord, Ord>(n, wt, u, v);
    }
  protected:
    // Members.
    double H;
  };

  class LinearFormSurfBottom : public VectorFormSurf<double>
  {
  public:
    LinearFormSurfBottom(int i, double H, double T0, std::string area) : VectorFormSurf<double>(i, area), H(H), T0(T0) {}

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return T0 * H * int_v<double>(n, wt, v);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return T0 * H * int_v<Ord>(n, wt, v);
    }
  protected:
    // Members.
    double H;
    double T0;
  };
};


class ErrorEstimatorFormInterface : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
public:
  ErrorEstimatorFormInterface(int i) : ErrorEstimatorForm(i, H2D_DG_INNER_EDGE) {}

  virtual double value(int n, double *wt, Func<double> *u_ext[],
             Func<double> *u, Geom<double> *e,
             ExtData<double> *ext)
  {
    return kelly_interface_estimator<double, double>(n, wt, u_ext, u, e, ext);
  }
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                   Func<Ord> *u, Geom<Ord> *e,
                   ExtData<Ord> *ext)
  {
    return kelly_interface_estimator<Ord, Ord>(n, wt, u_ext, u, e, ext);
  }

private:
  template<typename Real, typename Scalar>
  Scalar kelly_interface_estimator(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                   Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = Scalar(0);

    double EPS_C = adapt->get_element_markers_conversion()->get_user_marker(e->elem_marker).marker == MATERIAL_1 ? EPS_1 : EPS_2;
    double EPS_N = adapt->get_element_markers_conversion()->get_user_marker(e->get_neighbor_marker()).marker == MATERIAL_1 ? EPS_1 : EPS_2;

    for (int i = 0; i < n; i++)
      result += wt[i] * sqr( e->nx[i] * (EPS_C*u->get_dx_central(i) - EPS_N*u->get_dx_neighbor(i)) +
                             e->ny[i] * (EPS_C*u->get_dy_central(i) - EPS_N*u->get_dy_neighbor(i))  );

    return result;  // Multiplication by element diameter will be done automatically by the KellyTypeAdapt class.
                    // This allows to call this function only once for each interface and add the correctly scaled
                    // result to the total error estimate for both elements sharing the interface.
  }
};

class ErrorEstimatorFormZeroNeumann : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
public:
  ErrorEstimatorFormZeroNeumann(int i, std::string area) : ErrorEstimatorForm(i, area) {}

  virtual double value(int n, double *wt, Func<double> *u_ext[],
             Func<double> *u, Geom<double> *e,
             ExtData<double> *ext)
  {
    return kelly_zero_neumann_boundary_estimator<double, double>(n, wt, u_ext, u, e, ext);
  }
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                   Func<Ord> *u, Geom<Ord> *e,
                   ExtData<Ord> *ext)
  {
    return kelly_zero_neumann_boundary_estimator<Ord, Ord>(n, wt, u_ext, u, e, ext);
  }

private:
  template<typename Real, typename Scalar>
  Scalar kelly_zero_neumann_boundary_estimator(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                               Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = Scalar(0);

    double EPS_C = adapt->get_element_markers_conversion()->get_user_marker(e->elem_marker).marker == MATERIAL_1 ? EPS_1 : EPS_2;

    for (int i = 0; i < n; i++)
      result += wt[i] * sqr( e->nx[i] * EPS_C * u->dx[i] + e->ny[i] * EPS_C * u->dy[i] );

    return e->diam * result;
  }
};

class ErrorEstimatorFormNeumann : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
public:
  ErrorEstimatorFormNeumann(int i, std::string area) : ErrorEstimatorForm(i, area) {}

  virtual double value(int n, double *wt, Func<double> *u_ext[],
             Func<double> *u, Geom<double> *e,
             ExtData<double> *ext)
  {
    return kelly_neumann_boundary_estimator<double, double>(n, wt, u_ext, u, e, ext);
  }
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                   Func<Ord> *u, Geom<Ord> *e,
                   ExtData<Ord> *ext)
  {
    return kelly_neumann_boundary_estimator<Ord, Ord>(n, wt, u_ext, u, e, ext);
  }

private:
  template<typename Real, typename Scalar>
Scalar kelly_neumann_boundary_estimator(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = Scalar(0);

  double EPS_C = adapt->get_element_markers_conversion()->get_user_marker(e->elem_marker).marker == MATERIAL_1 ? EPS_1 : EPS_2;

  for (int i = 0; i < n; i++)
    result += wt[i] * sqr( CONST_GAMMA_OUTER - e->nx[i] * EPS_C * u->dx[i] - e->ny[i] * EPS_C * u->dy[i] );

  return e->diam * result;
}
};

class ErrorEstimatorFormNewton : public KellyTypeAdapt<double>::ErrorEstimatorForm
{
public:
  ErrorEstimatorFormNewton(int i, std::string area) : ErrorEstimatorForm(i, area) {}

  virtual double value(int n, double *wt, Func<double> *u_ext[],
             Func<double> *u, Geom<double> *e,
             ExtData<double> *ext)
  {
    return kelly_newton_boundary_estimator<double, double>(n, wt, u_ext, u, e, ext);
  }
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                   Func<Ord> *u, Geom<Ord> *e,
                   ExtData<Ord> *ext)
  {
    return kelly_newton_boundary_estimator<Ord, Ord>(n, wt, u_ext, u, e, ext);
  }

private:
  template<typename Real, typename Scalar>
  Scalar kelly_newton_boundary_estimator(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u,
                                         Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = Scalar(0);

    double EPS_C = adapt->get_element_markers_conversion()->get_user_marker(e->elem_marker).marker == MATERIAL_1 ? EPS_1 : EPS_2;

    for (int i = 0; i < n; i++)
      result += wt[i] * sqr( H*u->val[i] - H*T0 - e->nx[i] * EPS_C * u->dx[i] - e->ny[i] * EPS_C * u->dy[i] );

    return e->diam * result;
  }
};
