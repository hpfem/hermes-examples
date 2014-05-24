#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::WeakFormsElasticity;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

class CustomWeakFormLinearElasticity : public WeakForm<double>
{
public:
  CustomWeakFormLinearElasticity(double E, double nu, Hermes2DFunction<double>* Eo11d, Hermes2DFunction<double>* Eo12d, Hermes2DFunction<double>* Eo22d, Hermes2DFunction<double>* Eo33d, double Eo11, double Eo12, double Eo22, double Eo33, double rho_g,
      std::string surface_force_bdy, double f0, double f1); //
private:
  class CustomJacobianElast00 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianElast00(int i, int j, double lambda, double mu, double kappa, double Eo11, double Eo12, double Eo22, double Eo33)
    : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), kappa(kappa), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33) {};
    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, Func<Ord>* *ext) const;
    virtual MatrixFormVol<double>* clone() const;
    protected:
      double lambda;
      double mu;
      double kappa;
      double Eo11;
      double Eo12;
      double Eo22;
      double Eo33;
  };
  class CustomJacobianElast01 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianElast01(int i, int j, double lambda, double mu, double kappa, double Eo11, double Eo12, double Eo22, double Eo33)
    : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), kappa(kappa), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33) {};
    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, Func<Ord>* *ext) const;
    virtual MatrixFormVol<double>* clone() const;
    protected:
      double lambda;
      double mu;
      double kappa;
      double Eo11;
      double Eo12;
      double Eo22;
      double Eo33;
  };
  class CustomJacobianElast02 : public MatrixFormVol<double>
    {
    public:
      CustomJacobianElast02(int i, int j, double lambda, double mu, double kappa, double Eo11, double Eo12, double Eo22, double Eo33)
      : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), kappa(kappa), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33) {};
      virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, Func<Ord>* *ext) const;
      virtual MatrixFormVol<double>* clone() const;
      protected:
        double lambda;
        double mu;
        double kappa;
        double Eo11;
        double Eo12;
        double Eo22;
        double Eo33;
    };
  class CustomJacobianElast10 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianElast10(int i, int j, double lambda, double mu, double kappa, double Eo11, double Eo12, double Eo22, double Eo33)
    : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), kappa(kappa), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33) {};
    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, Func<Ord>* *ext) const;
    virtual MatrixFormVol<double>* clone() const;
    protected:
      double lambda;
      double mu;
      double kappa;
      double Eo11;
      double Eo12;
      double Eo22;
      double Eo33;
  };
  class CustomJacobianElast11 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianElast11(int i, int j, double lambda, double mu, double kappa, double Eo11, double Eo12, double Eo22, double Eo33)
    : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), kappa(kappa), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33) {};
    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, Func<Ord>* *ext) const;
    virtual MatrixFormVol<double>* clone() const;
    protected:
      double lambda;
      double mu;
      double kappa;
      double Eo11;
      double Eo12;
      double Eo22;
      double Eo33;
  };
  class CustomJacobianElast12 : public MatrixFormVol<double>
    {
    public:
      CustomJacobianElast12(int i, int j, double lambda, double mu, double kappa, double Eo11, double Eo12, double Eo22, double Eo33)
      : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), kappa(kappa), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33) {};
      virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, Func<Ord>* *ext) const;
      virtual MatrixFormVol<double>* clone() const;
      protected:
        double lambda;
        double mu;
        double kappa;
        double Eo11;
        double Eo12;
        double Eo22;
        double Eo33;
  };
  class CustomJacobianElast20 : public MatrixFormVol<double>
    {
    public:
      CustomJacobianElast20(int i, int j, double lambda, double mu, double kappa, double Eo11, double Eo12, double Eo22, double Eo33)
      : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), kappa(kappa), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33) {};
      virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, Func<Ord>* *ext) const;
      virtual MatrixFormVol<double>* clone() const;
      protected:
        double lambda;
        double mu;
        double kappa;
        double Eo11;
        double Eo12;
        double Eo22;
        double Eo33;
  };
  class CustomJacobianElast21 : public MatrixFormVol<double>
    {
    public:
      CustomJacobianElast21(int i, int j, double lambda, double mu, double kappa, double Eo11, double Eo12, double Eo22, double Eo33)
      : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), kappa(kappa), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33) {};
      virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, Func<Ord>* *ext) const;
      virtual MatrixFormVol<double>* clone() const;
      protected:
        double lambda;
        double mu;
        double kappa;
        double Eo11;
        double Eo12;
        double Eo22;
        double Eo33;
  };
  class CustomJacobianElast22 : public MatrixFormVol<double>
    {
    public:
      CustomJacobianElast22(int i, int j, double lambda, double mu, double kappa, double Eo11, double Eo12, double Eo22, double Eo33)
      : MatrixFormVol<double>(i, j), lambda(lambda), mu(mu), kappa(kappa), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33) {};
      virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                           Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, Func<Ord>* *ext) const;
      virtual MatrixFormVol<double>* clone() const;
      protected:
        double lambda;
        double mu;
        double kappa;
        double Eo11;
        double Eo12;
        double Eo22;
        double Eo33;
  };
  class CustomVectorRes0 : public VectorFormVol<double>
  {
  public:
	  CustomVectorRes0(int i, double lambda, double mu, double kappa, Hermes2DFunction<double>* Eo11d, Hermes2DFunction<double>* Eo12d, Hermes2DFunction<double>* Eo22d, Hermes2DFunction<double>* Eo33d, double Eo11, double Eo12, double Eo22, double Eo33)
      : VectorFormVol<double>(i), lambda(lambda), mu(mu), kappa(kappa), Eo11d(Eo11d), Eo12d(Eo12d), Eo22d(Eo22d), Eo33d(Eo33d), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[],
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, Func<Ord>* *ext) const;
    virtual VectorFormVol<double>* clone() const;
    private:
    double lambda;
    double mu;
    double kappa;
    Hermes2DFunction<double>* Eo11d;
    Hermes2DFunction<double>* Eo12d;
    Hermes2DFunction<double>* Eo22d;
    Hermes2DFunction<double>* Eo33d;
    double Eo11;
    double Eo12;
    double Eo22;
    double Eo33;
  };
  class CustomVectorRes1 : public VectorFormVol<double>
  {
  public:
	  CustomVectorRes1(int i, double lambda, double mu, double kappa, Hermes2DFunction<double>* Eo11d, Hermes2DFunction<double>* Eo12d, Hermes2DFunction<double>* Eo22d, Hermes2DFunction<double>* Eo33d, double Eo11, double Eo12, double Eo22, double Eo33)
      : VectorFormVol<double>(i), lambda(lambda), mu(mu), kappa(kappa), Eo11d(Eo11d), Eo12d(Eo12d), Eo22d(Eo22d), Eo33d(Eo33d), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[],
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, Func<Ord>* *ext) const;
    virtual VectorFormVol<double>* clone() const;
    private:
    double lambda;
    double mu;
    double kappa;
    Hermes2DFunction<double>* Eo11d;
    Hermes2DFunction<double>* Eo12d;
    Hermes2DFunction<double>* Eo22d;
    Hermes2DFunction<double>* Eo33d;
    double Eo11;
    double Eo12;
    double Eo22;
    double Eo33;
  };
  class CustomVectorRes2 : public VectorFormVol<double>
  {
  public:
	  CustomVectorRes2(int i, double lambda, double mu, double kappa, Hermes2DFunction<double>* Eo11d, Hermes2DFunction<double>* Eo12d, Hermes2DFunction<double>* Eo22d, Hermes2DFunction<double>* Eo33d, double Eo11, double Eo12, double Eo22, double Eo33)
      : VectorFormVol<double>(i), lambda(lambda), mu(mu), kappa(kappa), Eo11d(Eo11d), Eo12d(Eo12d), Eo22d(Eo22d), Eo33d(Eo33d), Eo11(Eo11), Eo12(Eo12), Eo22(Eo22), Eo33(Eo33)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[],
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;
    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, Func<Ord>* *ext) const;
    virtual VectorFormVol<double>* clone() const;
    private:
    double lambda;
    double mu;
    double kappa;
    Hermes2DFunction<double>* Eo11d;
    Hermes2DFunction<double>* Eo12d;
    Hermes2DFunction<double>* Eo22d;
    Hermes2DFunction<double>* Eo33d;
    double Eo11;
    double Eo12;
    double Eo22;
    double Eo33;
  };
};

// Custom function Eo11
class CustomExactFunctionEo11
{
public:
  CustomExactFunctionEo11(double dist) : dist(dist) {};
  double val(double x, double y);
  double dx(double x, double y);
  double ddxx(double x, double y);
protected :
  double dist;
};

// Custom function Eo12
class CustomExactFunctionEo12
{
public:
  CustomExactFunctionEo12(double dist) : dist(dist) {};
  double val(double x, double y);
  double dx(double x, double y);
protected :
  double dist;
};

// Custom function Eo22
class CustomExactFunctionEo22
{
public:
  CustomExactFunctionEo22(double dist) : dist(dist) {};
  double val(double x, double y);
  double dx(double x, double y);
protected :
  double dist;
};

// Custom function Eo33
class CustomExactFunctionEo33
{
public:
  CustomExactFunctionEo33(double dist) : dist(dist) {};
  double val(double x, double y);
  double dx(double x, double y);
protected :
  double dist;
};

// Eo11 Hermes2DFunction
class CustomEo11: public Hermes2DFunction<double>
{
public:
  CustomEo11(double dist);
  virtual double value(double x, double y) const;
  virtual Ord value(Ord x, Ord y) const;
  ~CustomEo11();
  CustomExactFunctionEo11* cefEo11;
protected :
  double dist;
};

// Eo12 Hermes2DFunction
class CustomEo12: public Hermes2DFunction<double>
{
public:
  CustomEo12(double dist);
  virtual double value(double x, double y) const;
  virtual Ord value(Ord x, Ord y) const;
  ~CustomEo12();
  CustomExactFunctionEo12* cefEo12;
protected :
  double dist;
};

// Eo22 Hermes2DFunction
class CustomEo22: public Hermes2DFunction<double>
{
public:
  CustomEo22(double dist);
  virtual double value(double x, double y) const;
  virtual Ord value(Ord x, Ord y) const;
  ~CustomEo22();
  CustomExactFunctionEo22* cefEo22;
protected :
  double dist;
};

// Eo33 Hermes2DFunction
class CustomEo33: public Hermes2DFunction<double>
{
public:
  CustomEo33(double dist);
  virtual double value(double x, double y) const;
  virtual Ord value(Ord x, Ord y) const;
  ~CustomEo33();
  CustomExactFunctionEo33* cefEo33;
protected :
  double dist;
};

// Exact solution Eo11
class ExactSolutionEo11 : public ExactSolutionScalar<double>
{
public:
  ExactSolutionEo11(MeshSharedPtr mesh, double dist);
  virtual double value(double x, double y) const;
  virtual void derivatives(double x, double y, double& dx, double& dy) const;
  virtual Ord ord(double x, double y) const;
  ~ExactSolutionEo11();
  virtual MeshFunction<double>* clone() const;
  CustomExactFunctionEo11* cefEo11;
protected :
  double dist;
};

// Exact solution Eo12
class ExactSolutionEo12 : public ExactSolutionScalar<double>
{
public:
  ExactSolutionEo12(MeshSharedPtr mesh, double dist);
  virtual double value(double x, double y) const;
  virtual void derivatives(double x, double y, double& dx, double& dy) const;
  virtual Ord ord(double x, double y) const;
  ~ExactSolutionEo12();
  virtual MeshFunction<double>* clone() const;
  CustomExactFunctionEo12* cefEo12;
protected :
  double dist;
};

// Exact solution Eo22
class ExactSolutionEo22 : public ExactSolutionScalar<double>
{
public:
  ExactSolutionEo22(MeshSharedPtr mesh, double dist);
  virtual double value(double x, double y) const;
  virtual void derivatives(double x, double y, double& dx, double& dy) const;
  virtual Ord ord(double x, double y) const;
  ~ExactSolutionEo22();
  virtual MeshFunction<double>* clone() const;
  CustomExactFunctionEo22* cefEo22;
protected :
  double dist;
};

// Exact solution Eo33
class ExactSolutionEo33 : public ExactSolutionScalar<double>
{
public:
  ExactSolutionEo33(MeshSharedPtr mesh, double dist);
  virtual double value(double x, double y) const;
  virtual void derivatives(double x, double y, double& dx, double& dy) const;
  virtual Ord ord(double x, double y) const;
  ~ExactSolutionEo33();
  virtual MeshFunction<double>* clone() const;
  CustomExactFunctionEo33* cefEo33;
protected :
  double dist;
};

// Exact solution TrEo
class ExactSolutionTrEo : public ExactSolutionScalar<double>
{
public:
  ExactSolutionTrEo(MeshSharedPtr mesh, double Eo11, double Eo22, double Eo33);
  virtual double value(double x, double y) const;
  virtual void derivatives(double x, double y, double& dx, double& dy) const;
  virtual Ord ord(double x, double y) const;
  ~ExactSolutionTrEo();
  virtual MeshFunction<double>* clone() const;
  CustomExactFunctionEo11* cefEo11;
  CustomExactFunctionEo22* cefEo22;
  CustomExactFunctionEo33* cefEo33;
protected :
  double Eo11;
  double Eo22;
  double Eo33;
};

// Custom filter S11
class CustomFilterS11 : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterS11(Hermes::vector<MeshFunctionSharedPtr<double>> solutions, double lambda, double mu, CustomEo11* Eo11d, CustomEo12* Eo12d, CustomEo22* Eo22d, CustomEo33* Eo33d) : Hermes::Hermes2D::DXDYFilter<double>(solutions), lambda(lambda), mu(mu), Eo11d(Eo11d), Eo12d(Eo12d), Eo22d(Eo22d), Eo33d(Eo33d)
{
}
virtual MeshFunction<double>* clone() const
{
  Hermes::vector<MeshFunctionSharedPtr<double> > slns;
  for(int i = 0; i < this->num; i++)
  {
    slns.push_back(this->sln[i]->clone());
  }
  CustomFilterS11* filter = new CustomFilterS11(slns, lambda, mu, Eo11d, Eo12d, Eo22d, Eo33d);
  return filter;
}
private:
virtual void filter_fn(int n, double* x, double* y, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);
double mu, lambda;
CustomEo11* Eo11d;
CustomEo12* Eo12d;
CustomEo22* Eo22d;
CustomEo33* Eo33d;
};

// Custom filter S12
class CustomFilterS12 : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterS12(Hermes::vector<MeshFunctionSharedPtr<double>> solutions, double lambda, double mu, CustomEo11* Eo11d, CustomEo12* Eo12d, CustomEo22* Eo22d, CustomEo33* Eo33d) : Hermes::Hermes2D::DXDYFilter<double>(solutions), lambda(lambda), mu(mu), Eo11d(Eo11d), Eo12d(Eo12d), Eo22d(Eo22d), Eo33d(Eo33d)
{
}
virtual MeshFunction<double>* clone() const
{
  Hermes::vector<MeshFunctionSharedPtr<double> > slns;
  for(int i = 0; i < this->num; i++)
  {
    slns.push_back(this->sln[i]->clone());
  }
  CustomFilterS12* filter = new CustomFilterS12(slns, lambda, mu, Eo11d, Eo12d, Eo22d, Eo33d);
  return filter;
}
private:
virtual void filter_fn(int n, double* x, double* y, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);
double mu, lambda;
CustomEo11* Eo11d;
CustomEo12* Eo12d;
CustomEo22* Eo22d;
CustomEo33* Eo33d;
};

// Custom filter S22
class CustomFilterS22 : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterS22(Hermes::vector<MeshFunctionSharedPtr<double>> solutions, double lambda, double mu, CustomEo11* Eo11d, CustomEo12* Eo12d, CustomEo22* Eo22d, CustomEo33* Eo33d) : Hermes::Hermes2D::DXDYFilter<double>(solutions), lambda(lambda), mu(mu), Eo11d(Eo11d), Eo12d(Eo12d), Eo22d(Eo22d), Eo33d(Eo33d)
{
}
virtual MeshFunction<double>* clone() const
{
  Hermes::vector<MeshFunctionSharedPtr<double> > slns;
  for(int i = 0; i < this->num; i++)
  {
    slns.push_back(this->sln[i]->clone());
  }
  CustomFilterS22* filter = new CustomFilterS22(slns, lambda, mu, Eo11d, Eo12d, Eo22d, Eo33d);
  return filter;
}
private:
virtual void filter_fn(int n, double* x, double* y, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);
double mu, lambda;
CustomEo11* Eo11d;
CustomEo12* Eo12d;
CustomEo22* Eo22d;
CustomEo33* Eo33d;
};

// Custom filter S33
class CustomFilterS33 : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilterS33(Hermes::vector<MeshFunctionSharedPtr<double>> solutions, double lambda, double mu, CustomEo11* Eo11d, CustomEo12* Eo12d, CustomEo22* Eo22d, CustomEo33* Eo33d) : Hermes::Hermes2D::DXDYFilter<double>(solutions), lambda(lambda), mu(mu), Eo11d(Eo11d), Eo12d(Eo12d), Eo22d(Eo22d), Eo33d(Eo33d)
{
}
virtual MeshFunction<double>* clone() const
{
  Hermes::vector<MeshFunctionSharedPtr<double> > slns;
  for(int i = 0; i < this->num; i++)
  {
    slns.push_back(this->sln[i]->clone());
  }
  CustomFilterS33* filter = new CustomFilterS33(slns, lambda, mu, Eo11d, Eo12d, Eo22d, Eo33d);
  return filter;
}
private:
virtual void filter_fn(int n, double* x, double* y, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);
double mu, lambda;
CustomEo11* Eo11d;
CustomEo12* Eo12d;
CustomEo22* Eo22d;
CustomEo33* Eo33d;
};

// Custom filter von Mises (with distortion)
class CustomFilter_vM : public Hermes::Hermes2D::DXDYFilter<double>
{
public:
  CustomFilter_vM(Hermes::vector<MeshFunctionSharedPtr<double>> solutions, double lambda, double mu, CustomEo11* Eo11d, CustomEo12* Eo12d, CustomEo22* Eo22d, CustomEo33* Eo33d) : Hermes::Hermes2D::DXDYFilter<double>(solutions), lambda(lambda), mu(mu), Eo11d(Eo11d), Eo12d(Eo12d), Eo22d(Eo22d), Eo33d(Eo33d)
{
}
virtual MeshFunction<double>* clone() const
{
  Hermes::vector<MeshFunctionSharedPtr<double> > slns;
  for(int i = 0; i < this->num; i++)
  {
    slns.push_back(this->sln[i]->clone());
  }
  CustomFilter_vM* filter = new CustomFilter_vM(slns, lambda, mu, Eo11d, Eo12d, Eo22d, Eo33d);
  return filter;
}
private:
virtual void filter_fn(int n, double* x, double* y, Hermes::vector<double *> values, Hermes::vector<double *> dx, Hermes::vector<double *> dy, double* rslt, double* rslt_dx, double* rslt_dy);
double mu, lambda;
CustomEo11* Eo11d;
CustomEo12* Eo12d;
CustomEo22* Eo22d;
CustomEo33* Eo33d;
};
