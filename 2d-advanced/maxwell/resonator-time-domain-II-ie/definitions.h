#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Weak forms */

class CustomWeakFormWaveIE : public WeakForm<double>
{
public:

  CustomWeakFormWaveIE(double tau, double c_squared, Solution<double>* E_prev_sln, Solution<double>* F_prev_sln);

private:
  class MatrixFormVolWave_0_0 : public MatrixFormVol<double>
  {
  public:
    MatrixFormVolWave_0_0(double tau) : MatrixFormVol<double>(0, 0), tau(tau) { this->setSymFlag(HERMES_SYM);};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;

    double tau;
  };

  class MatrixFormVolWave_0_1 : public MatrixFormVol<double>
  {
  public:
    MatrixFormVolWave_0_1() : MatrixFormVol<double>(0, 1) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, Func<Ord>* *ext) const;
    MatrixFormVol<double>* clone() const;
  };

  class MatrixFormVolWave_1_0 : public MatrixFormVol<double>
  {
  public:
    MatrixFormVolWave_1_0(double c_squared) 
          : MatrixFormVol<double>(1, 0), c_squared(c_squared) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;
    double c_squared;
  };

  class MatrixFormVolWave_1_1 : public MatrixFormVol<double>
  {
  public:
    MatrixFormVolWave_1_1(double tau) : MatrixFormVol<double>(1, 1), tau(tau) { this->setSymFlag(HERMES_SYM); };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, Func<Ord>* *ext) const;

    MatrixFormVol<double>* clone() const;
    double tau;
  };

  class VectorFormVolWave_0 : public VectorFormVol<double>
  {
  public:
    VectorFormVolWave_0(double tau) : VectorFormVol<double>(0), tau(tau) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
                    Func<Ord>* *ext) const;
    
    VectorFormVol<double>* clone() const;
    double tau;
  };

  class VectorFormVolWave_1 : public VectorFormVol<double>
  {
  public:
  VectorFormVolWave_1(double tau, double c_squared) 
    : VectorFormVol<double>(1), tau(tau), c_squared(c_squared) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

    VectorFormVol<double>* clone() const;
    double tau, c_squared;
  };
};

/* Initial condition for E */

class CustomInitialConditionWave : public ExactSolutionVector<double>
{
public:
  CustomInitialConditionWave(const Mesh* mesh) : ExactSolutionVector<double>(mesh) {};

  virtual Scalar2<double> value (double x, double y) const;

  virtual void derivatives (double x, double y, Scalar2<double>& dx, Scalar2<double>& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
  MeshFunction<double>* clone() const;
};

