#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Global function alpha */

double alpha(double omega, double k);

/* Initial condition for E */

class CustomInitialConditionE : public ExactSolutionVector<double>
{
public:
 CustomInitialConditionE(Mesh* mesh, double time, double omega, double k_x, double k_y) 
     : ExactSolutionVector<double>(mesh), time(time), omega(omega), k_x(k_x), k_y(k_y) {};

  virtual Scalar2<double> value (double x, double y) const;

  virtual void derivatives (double x, double y, Scalar2<double>& dx, Scalar2<double>& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
  
  double time, omega, k_x, k_y;
};

/* Initial condition for H */

class CustomInitialConditionH : public ExactSolutionScalar<double>
{
public:
  CustomInitialConditionH(Mesh* mesh, double time, double omega, double k_x, double k_y) 
      : ExactSolutionScalar<double>(mesh), time(time), omega(omega), k_x(k_x), k_y(k_y) {};

  virtual double value (double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double time, omega, k_x, k_y;
};

/* Initial condition for P */

class CustomInitialConditionP : public ExactSolutionVector<double>
{
public:
  CustomInitialConditionP(Mesh* mesh, double time, double omega, double k_x, double k_y) 
      : ExactSolutionVector<double>(mesh), time(time), omega(omega), k_x(k_x), k_y(k_y) {};

  virtual Scalar2<double> value (double x, double y) const;

  virtual void derivatives (double x, double y, Scalar2<double>& dx, Scalar2<double>& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  double time, omega, k_x, k_y;
};

/* Weak forms */

class CustomWeakFormMD : public WeakForm<double>
{
public:

  CustomWeakFormMD(double omega, double k_x, double k_y, double mu_0, 
      double eps_0, double eps_inf, double eps_q, double tau);

private:
  class MatrixFormVolMD_0_0 : public MatrixFormVol<double>
  {
  public:
  MatrixFormVolMD_0_0(double eps_q, double tau) 
      : MatrixFormVol<double>(0, 0, HERMES_ANY, HERMES_NONSYM), eps_q(eps_q), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();

    double eps_q, tau;
  };

  class MatrixFormVolMD_0_1 : public MatrixFormVol<double>
  {
  public:
  MatrixFormVolMD_0_1(double eps_0, double eps_inf) : MatrixFormVol<double>(0, 1, HERMES_ANY, HERMES_NONSYM), eps_0(eps_0), eps_inf(eps_inf) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();

    double eps_0, eps_inf;
  };

  class MatrixFormVolMD_0_2 : public MatrixFormVol<double>
  {
  public:
  MatrixFormVolMD_0_2(double eps_0, double eps_inf, double tau) : MatrixFormVol<double>(0, 2, HERMES_ANY, HERMES_NONSYM), eps_0(eps_0), eps_inf(eps_inf), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();

    double eps_0, eps_inf, tau;
  };

  class MatrixFormVolMD_1_0 : public MatrixFormVol<double>
  {
  public:
  MatrixFormVolMD_1_0(double mu_0) : MatrixFormVol<double>(1, 0, HERMES_ANY, HERMES_NONSYM), mu_0(mu_0) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();

    double mu_0;
  };

  class MatrixFormVolMD_2_0 : public MatrixFormVol<double>
  {
  public:
  MatrixFormVolMD_2_0(double eps_0, double eps_inf, double eps_q, double tau) : MatrixFormVol<double>(2, 0, HERMES_ANY, HERMES_NONSYM), eps_0(eps_0), eps_inf(eps_inf), eps_q(eps_q), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();

    double eps_0, eps_inf, eps_q, tau;
  };

  class MatrixFormVolMD_2_2 : public MatrixFormVol<double>
  {
  public:
  MatrixFormVolMD_2_2(double tau) : MatrixFormVol<double>(2, 2, HERMES_ANY, HERMES_NONSYM), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
        Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
        Geom<Ord> *e, ExtData<Ord> *ext) const;

    virtual MatrixFormVol<double>* clone();

    double tau;
  };

  class VectorFormVolMD_0 : public VectorFormVol<double>
  {
  public:
  VectorFormVolMD_0(double eps_0, double eps_inf, double eps_q, double tau) : VectorFormVol<double>(0), eps_0(eps_0), eps_inf(eps_inf), eps_q(eps_q), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
        ExtData<Ord> *ext) const;

    virtual VectorFormVol<double>* clone();

    double eps_0, eps_inf, eps_q, tau;
  };

  class VectorFormVolMD_1 : public VectorFormVol<double>
  {
  public:
  VectorFormVolMD_1(double mu_0) : VectorFormVol<double>(1), mu_0(mu_0) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
        ExtData<Ord> *ext) const;

    virtual VectorFormVol<double>* clone();

    double mu_0;
  };

  class VectorFormVolMD_2 : public VectorFormVol<double>
  {
  public:
  VectorFormVolMD_2(double eps_0, double eps_inf, double eps_q, double tau) : VectorFormVol<double>(2), eps_0(eps_0), eps_inf(eps_inf), eps_q(eps_q), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
        Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
        ExtData<Ord> *ext) const;

    virtual VectorFormVol<double>* clone();

    double eps_0, eps_inf, eps_q, tau;
  };
};
