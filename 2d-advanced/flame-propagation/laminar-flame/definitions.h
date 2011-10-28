#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(double Le, double alpha, double beta, double kappa, double x1, 
      Hermes::Hermes2D::Filter<double>* omega, Hermes::Hermes2D::Filter<double>* d_omega_dT, Hermes::Hermes2D::Filter<double>* d_omega_dC);

  ~CustomWeakForm() {};

private:
  class JacobianFormVol_0_0 : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol_0_0() 
            : MatrixFormVol<double>(0, 0, HERMES_ANY, HERMES_SYM) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();
  };
  
  class JacobianFormVol_0_1 : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol_0_1() 
            : MatrixFormVol<double>(0, 1, HERMES_ANY, HERMES_SYM) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();
  };

  class JacobianFormVol_1_0 : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol_1_0() 
            : MatrixFormVol<double>(1, 0, HERMES_ANY, HERMES_SYM) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();
  };

  class JacobianFormVol_1_1 : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol_1_1(double Le) 
            : MatrixFormVol<double>(1, 1, HERMES_ANY, HERMES_SYM), Le(Le) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();

    double Le;
  };

  class JacobianFormSurf_0_0 : public MatrixFormSurf<double>
  {
  public:
    JacobianFormSurf_0_0(std::string bnd_marker, double kappa) 
            : MatrixFormSurf<double>(0, 0, bnd_marker), kappa(kappa) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormSurf<double>* clone();

    double kappa;
  };

  class ResidualFormVol_0 : public VectorFormVol<double>
  {
  public:
    ResidualFormVol_0() 
            : VectorFormVol<double>(0)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual VectorFormVol<double>* clone();
  };

  class ResidualFormVol_1 : public VectorFormVol<double>
  {
  public:
    ResidualFormVol_1(double Le) 
            : VectorFormVol<double>(1), Le(Le)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual VectorFormVol<double>* clone();

  private:
    double Le;
  };
  
  class ResidualFormSurf_0 : public VectorFormSurf<double>
  {
  public:
    ResidualFormSurf_0(std::string bnd_marker, double kappa) 
            : VectorFormSurf<double>(0, bnd_marker), kappa(kappa)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual VectorFormSurf<double>* clone();

  private:
    double kappa;
  };

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
};

class InitialSolutionTemperature : public ExactSolutionScalar<double>
{
public:
  InitialSolutionTemperature(Mesh* mesh, double x1) : ExactSolutionScalar<double>(mesh), x1(x1) {};

  virtual double value (double x, double y) const {
    return (x <= x1) ? 1.0 : exp(x1 - x);
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = (x <= x1) ? 0.0 : -exp(x1 - x);
    dy = 0.0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return -exp(x1 - x);
  }

  // Value.
  double x1;
};

class InitialSolutionConcentration : public ExactSolutionScalar<double>
{
public:
  InitialSolutionConcentration(Mesh* mesh, double x1, double Le) : ExactSolutionScalar<double>(mesh), x1(x1), Le(Le) {};

  virtual double value (double x, double y) const {
    return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x));
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = (x <= x1) ? 0.0 : Le * exp(x1 - x);
    dy = 0.0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return exp(Le*(x1 - x));
  }

  // Value.
  double x1, Le;
};

class CustomFilter : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilter(Hermes::vector<Solution<double>*> solutions, double Le, double alpha, double beta, double kappa, double x1) : Hermes::Hermes2D::DXDYFilter<double>(solutions), Le(Le), alpha(alpha), beta(beta), kappa(kappa), x1(x1)
  {
  }

private:
  virtual void filter_fn (int n, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
};

class CustomFilterDc : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterDc(Hermes::vector<Solution<double>*> solutions, double Le, double alpha, double beta, double kappa, double x1) : Hermes::Hermes2D::DXDYFilter<double>(solutions), Le(Le), alpha(alpha), beta(beta), kappa(kappa), x1(x1)
  {
  }

private:
  virtual void filter_fn (int n, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
};

class CustomFilterDt : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterDt(Hermes::vector<Solution<double>*> solutions, double Le, double alpha, double beta, double kappa, double x1) : Hermes::Hermes2D::DXDYFilter<double>(solutions), Le(Le), alpha(alpha), beta(beta), kappa(kappa), x1(x1)
  {
  }

private:
  virtual void filter_fn (int n, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);

  double Le;
  double alpha;
  double beta;
  double kappa;
  double x1;
};
