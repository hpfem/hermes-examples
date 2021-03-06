#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(MeshSharedPtr mesh, double x0, double x1, double y0, double y1, 
                      double* t_ptr, double s, double c)
    : ExactSolutionScalar<double>(mesh), x0(x0), x1(x1), y0(y0), y1(y1), t_ptr(t_ptr), s(s), c(c) {};

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord (double x, double y) const;

  double x0, x1, y0, y1, *t_ptr, s, c;
};

/* Custom function f */

class CustomFunction: public Hermes::Hermes2DFunction<double>
{
public:
  CustomFunction(double x0, double x1, double y0, double y1,
    double s, double c)
    : Hermes::Hermes2DFunction<double>(), x0(x0), x1(x1), y0(y0), y1(y1), s(s), c(c) {};

  virtual double value(double x, double y, double t) const;

  virtual Ord value(Ord x, Ord y) const;

  double x0, x1, y0, y1, s, c, *t_ptr;
};

class CustomVectorFormVol : public VectorFormVol<double>
{
public:
  CustomVectorFormVol(int i = 0, 
    Hermes::Hermes2DFunction<double>* coeff = nullptr, std::string area = HERMES_ANY,
    GeomType gt = HERMES_PLANAR);

  CustomVectorFormVol(int i, 
    Hermes::Hermes2DFunction<double>* coeff, std::vector<std::string> areas,
    GeomType gt = HERMES_PLANAR);

  ~CustomVectorFormVol();

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                       GeomVol<double> *e, Func<double>* *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    GeomVol<Ord> *e, Func<Ord>* *ext) const;

  virtual VectorFormVol<double>* clone() const;

private:
  Hermes::Hermes2DFunction<double>* coeff;
  GeomType gt;
};

class CustomWeakFormPoisson : public WeakFormsH1::DefaultWeakFormPoisson<double>
{
public:
  CustomWeakFormPoisson(
    Hermes::Hermes1DFunction<double>* coeff = nullptr,
    Hermes::Hermes2DFunction<double>* f = nullptr,
    std::string area = HERMES_ANY, 
    GeomType gt = HERMES_PLANAR);
};

