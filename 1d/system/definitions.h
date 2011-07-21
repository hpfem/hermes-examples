#include "hermes2d.h"

using namespace Hermes::Hermes2D;

/* Custom function that is used in the exact solution and in right-hand side */

class CustomExactFunction1
{
public:
  CustomExactFunction1() {};

  double val(double x);
  
  double dx(double x);
  
  double ddxx(double x);
};

class CustomExactFunction2
{
public:
  CustomExactFunction2(double K) : K(K) {};

  double val(double x);
  
  double dx(double x);
  
  double ddxx(double x);

  double K;
};

/* Right-hand side */

class CustomRightHandSide1: public HermesFunction<double>
{
public:
  CustomRightHandSide1(double K, double d_u, double sigma);

  virtual double value(double x, double y) const;

  virtual Hermes::Ord ord(Hermes::Ord x, Hermes::Ord y) const;

  ~CustomRightHandSide1();

  CustomExactFunction1* cef1;
  CustomExactFunction2* cef2;
  double d_u, sigma;
};

class CustomRightHandSide2: public HermesFunction<double>
{
public:
  CustomRightHandSide2(double K, double d_v);

  virtual double value(double x, double y) const;

  virtual Hermes::Ord ord(Hermes::Ord x, Hermes::Ord y) const;

  ~CustomRightHandSide2();

  CustomExactFunction1* cef1;
  CustomExactFunction2* cef2;
  double d_v;
};

/* Exact solution */

class ExactSolutionFitzHughNagumo1 : public ExactSolutionScalar<double>
{
public:
  ExactSolutionFitzHughNagumo1(Mesh* mesh);

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Hermes::Ord ord(Hermes::Ord x, Hermes::Ord y) const;

  ~ExactSolutionFitzHughNagumo1();

  CustomExactFunction1* cef1;
};

class ExactSolutionFitzHughNagumo2 : public ExactSolutionScalar<double>
{
public:
  ExactSolutionFitzHughNagumo2(Mesh* mesh, double K);

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Hermes::Ord ord(Hermes::Ord x, Hermes::Ord y) const;

  ~ExactSolutionFitzHughNagumo2();

  CustomExactFunction2* cef2;
};

/* Weak forms */

class CustomResidual1 : public VectorFormVol<double>
{
public:
  CustomResidual1(double d_u, double sigma, CustomRightHandSide1* g1)
    : VectorFormVol<double>(0, Hermes::HERMES_ANY), d_u(d_u), sigma(sigma), g1(g1) {};

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<double> *ext) const;

  virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
                  Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;

  virtual VectorFormVol<double>* clone();

private:
  double d_u;
  double sigma;
  CustomRightHandSide1* g1;
};

class CustomResidual2 : public VectorFormVol<double>
{
public:
  CustomResidual2(double d_v, CustomRightHandSide2* g2)
    : VectorFormVol<double>(1, Hermes::HERMES_ANY), d_v(d_v), g2(g2) {};

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<double> *ext) const;
  
  virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
                  Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const;
  
  virtual VectorFormVol<double>* clone();
  
private:
  double d_v;
  CustomRightHandSide2* g2;
};

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(CustomRightHandSide1* g1, CustomRightHandSide2* g2);
};
