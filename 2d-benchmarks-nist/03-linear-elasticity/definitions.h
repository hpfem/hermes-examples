#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

/* Exact solution */

class CustomExactSolutionU : public ExactSolutionScalar<double>
{
public:
  CustomExactSolutionU(Mesh* mesh, double E, double nu, double lambda, double Q) 
      : ExactSolutionScalar<double>(mesh), E(E), nu(nu), lambda(lambda), Q(Q) {};

  double get_angle(double y, double x) const;

  double d_theta_dx(double x, double y) const;

  double d_theta_dxd_theta_dx(double x, double y) const;

  double d_theta_dy(double x, double y) const;

  double d_theta_dyd_theta_dy(double x, double y) const;

  double d_theta_dxd_theta_dy(double x, double y) const;

  double r(double x, double y) const;

  double drdx(double x, double y) const;

  double drdxdrdx(double x, double y) const;

  double drdy(double x, double y) const;

  double drdydrdy(double x, double y) const;

  double drdxdrdy(double x, double y) const;

  double u_r(double x, double y) const;

  double du_rdx(double x, double y) const;

  double du_rdxdu_rdx(double x, double y) const;

  double du_rdy(double x, double y) const;

  double du_rdydu_rdy(double x, double y) const;

  double du_rdxdu_rdy(double x, double y) const;

  double v_r(double x, double y) const;

  double dv_rdx(double x, double y) const;

  double dv_rdxdv_rdx(double x, double y) const;

  double dv_rdy(double x, double y) const;

  double dv_rdydv_rdy(double x, double y) const;

  double dv_rdxdv_rdy(double x, double y) const;

  // \partial^2 u / \partial x^2 
  double dudxdudx(double x, double y) const;

  double dudydudy(double x, double y) const;

  double dudxdudy(double x, double y) const;

  // \partial^2 v / \partial x^2
  double dvdxdvdx(double x, double y) const;

  double dvdydvdy(double x, double y) const;

  double dvdxdvdy(double x, double y) const;

  // Exact solution u(x,y) and its derivatives.
  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord (Ord x, Ord y) const;

  double E, nu, lambda, Q;
};

class CustomExactSolutionV : public ExactSolutionScalar<double>
{
public:
  CustomExactSolutionV(Mesh* mesh, double E, double nu, double lambda, double Q) 
      : ExactSolutionScalar<double>(mesh), E(E), nu(nu), lambda(lambda), Q(Q) {};

  double get_angle(double y, double x) const;

  double d_theta_dx(double x, double y) const;

  double d_theta_dxd_theta_dx(double x, double y) const;

  double d_theta_dy(double x, double y) const;

  double d_theta_dyd_theta_dy(double x, double y) const;

  double d_theta_dxd_theta_dy(double x, double y) const;

  double r(double x, double y) const;

  double drdx(double x, double y) const;

  double drdxdrdx(double x, double y) const;

  double drdy(double x, double y) const;

  double drdydrdy(double x, double y) const;

  double drdxdrdy(double x, double y) const;

  double u_r(double x, double y) const;

  double du_rdx(double x, double y) const;

  double du_rdxdu_rdx(double x, double y) const;

  double du_rdy(double x, double y) const;

  double du_rdydu_rdy(double x, double y) const;

  double du_rdxdu_rdy(double x, double y) const;

  double v_r(double x, double y) const;

  double dv_rdx(double x, double y) const;

  double dv_rdxdv_rdx(double x, double y) const;

  double dv_rdy(double x, double y) const;

  double dv_rdydv_rdy(double x, double y) const;

  double dv_rdxdv_rdy(double x, double y) const;

  // \partial^2 u / \partial x^2 
  double dudxdudx(double x, double y) const;

  double dudydudy(double x, double y) const;

  double dudxdudy(double x, double y) const;

  // \partial^2 v / \partial x^2
  double dvdxdvdx(double x, double y) const;

  double dvdydvdy(double x, double y) const;

  double dvdxdvdy(double x, double y) const;

  // Exact solution u(x,y) and its derivatives.
  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord (Ord x, Ord y) const;

  double E, nu, lambda, Q;
};

/* Weak forms */

class CustomWeakFormElasticityNIST : public WeakForm<double>
{
public:
  CustomWeakFormElasticityNIST(double E, double nu, double mu, double lambda);

public:
  class CustomMatrixFormVolElasticityNIST_0_0 : public MatrixFormVol<double>
  {
  public:
    CustomMatrixFormVolElasticityNIST_0_0(int i, int j, double E, double nu)
        : MatrixFormVol<double>(i, j), E(E), nu(nu) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    double E, nu;
  };

  class CustomMatrixFormVolElasticityNIST_0_1 : public MatrixFormVol<double>
  {
  public:
    CustomMatrixFormVolElasticityNIST_0_1(int i, int j, double E, double nu)
        : MatrixFormVol<double>(i, j), E(E), nu(nu) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    double E, nu;
  };

  class CustomMatrixFormVolElasticityNIST_1_1 : public MatrixFormVol<double>
  {
  public:
    CustomMatrixFormVolElasticityNIST_1_1(int i, int j, double E, double nu)
        : MatrixFormVol<double>(i, j), E(E), nu(nu) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    double E, nu;
  };


  class CustomVectorFormVolElasticityNIST_0 : public VectorFormVol<double>
  {
  public:
    CustomVectorFormVolElasticityNIST_0(int i, double E, double nu)
        : VectorFormVol<double>(i), E(E), nu(nu) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[],
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    double E, nu;
  };

  class CustomVectorFormVolElasticityNIST_1 : public VectorFormVol<double>
  {
  public:
    CustomVectorFormVolElasticityNIST_1(int i, double E, double nu)
        : VectorFormVol<double>(i), E(E), nu(nu) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[],
        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[],
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    double E, nu;
  };
};
