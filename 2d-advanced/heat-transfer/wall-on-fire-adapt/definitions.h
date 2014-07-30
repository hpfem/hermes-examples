#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Space-dependent thermal conductivity */
double lambda(double x, double y);

/* Time-dependent fire temperature */
template<typename Real>
Real T_fire_x(Real x);

template<typename Real>
Real T_fire_t(Real t);

/* Weak forms */

class CustomWeakFormHeatRK : public WeakForm < double >
{
public:
  CustomWeakFormHeatRK(std::string bdy_fire, std::string bdy_air,
    double alpha_fire, double alpha_air, double rho, double heatcap,
    double temp_ext_air, double temp_init, double* current_time_ptr);

private:
  // This form is custom since it contains space-dependent thermal conductivity.
  class CustomJacobianVol : public MatrixFormVol < double >
  {
  public:
    CustomJacobianVol(int i, int j, double rho, double heatcap)
      : MatrixFormVol<double>(i, j), rho(rho), heatcap(heatcap) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, GeomVol<double> *e,
      Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const;

    // This is needed for the rk_time_step_newton() method.
    virtual MatrixFormVol<double>* clone() const;

    double rho, heatcap;
  };

  // This form is custom since it contains space-dependent thermal conductivity.
  class CustomFormResidualVol : public VectorFormVol < double >
  {
  public:
    CustomFormResidualVol(int i, double rho, double heatcap)
      : VectorFormVol<double>(i), rho(rho), heatcap(heatcap) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomVol<double> *e,
      Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const;

    // Needed for the rk_time_step_newton() method.
    virtual VectorFormVol<double>* clone() const;

    double rho, heatcap;
  };

  // Custom due to time-dependent exterior temperature.
  class CustomFormResidualSurfFire : public VectorFormSurf < double >
  {
  public:
    CustomFormResidualSurfFire(int i, std::string area, double alpha_fire, double rho,
      double heatcap, double* current_time_ptr)
      : VectorFormSurf<double>(i), alpha_fire(alpha_fire), rho(rho),
      heatcap(heatcap), current_time_ptr(current_time_ptr) {
      this->set_area(area);
    };

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
      GeomSurf<Real> *e, Func<Scalar>* *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomSurf<double> *e,
      Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomSurf<Ord> *e, Func<Ord>* *ext) const;

    // Needed for the rk_time_step_newton() method.
    virtual VectorFormSurf<double>* clone() const;

    // Fire temperature as function of x and t.
    template<typename Real>
    Real T_fire(Real x, Real t) const;

    double alpha_fire, rho, heatcap, *current_time_ptr;
  };
};
