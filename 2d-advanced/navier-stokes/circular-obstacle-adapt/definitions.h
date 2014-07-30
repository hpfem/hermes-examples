#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

class WeakFormNSSimpleLinearization : public WeakForm < double >
{
public:
  WeakFormNSSimpleLinearization(bool Stokes, double Reynolds, double time_step, MeshFunctionSharedPtr<double>  x_vel_previous_time,
    MeshFunctionSharedPtr<double>  y_vel_previous_time);

  class BilinearFormSymVel : public MatrixFormVol < double >
  {
  public:
    BilinearFormSymVel(unsigned int i, unsigned int j, bool Stokes, double Reynolds, double time_step);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord>* *ext);

    virtual MatrixFormVol<double>* clone() const;

  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class BilinearFormNonsymVel : public MatrixFormVol < double >
  {
  public:
    BilinearFormNonsymVel(unsigned int i, unsigned int j, bool Stokes);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;

    virtual MatrixFormVol<double>* clone() const;
  protected:
    bool Stokes;
  };

  class BilinearFormNonsymXVelPressure : public MatrixFormVol < double >
  {
  public:
    BilinearFormNonsymXVelPressure(int i, int j);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;
    virtual MatrixFormVol<double>* clone() const;
  };

  class BilinearFormNonsymYVelPressure : public MatrixFormVol < double >
  {
  public:
    BilinearFormNonsymYVelPressure(int i, int j);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;
    virtual MatrixFormVol<double>* clone() const;
  };

  class VectorFormVolVel : public VectorFormVol < double >
  {
  public:
    VectorFormVolVel(int i, bool Stokes, double time_step);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomVol<double> *e,
      Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const;

    virtual VectorFormVol<double>* clone() const;
  protected:
    bool Stokes;
    double time_step;
  };

protected:
  bool Stokes;
  double Reynolds;
  double time_step;
  MeshFunctionSharedPtr<double>  x_vel_previous_time;
  MeshFunctionSharedPtr<double>  y_vel_previous_time;
};

class WeakFormNSNewton : public WeakForm < double >
{
public:
  WeakFormNSNewton(bool Stokes, double Reynolds, double time_step, MeshFunctionSharedPtr<double>  x_vel_previous_time,
    MeshFunctionSharedPtr<double>  y_vel_previous_time);

  class BilinearFormSymVel : public MatrixFormVol < double >
  {
  public:
    BilinearFormSymVel(unsigned int i, unsigned int j, bool Stokes, double Reynolds, double time_step);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormSymVel(this->i, this->j, this->Stokes, this->Reynolds, this->time_step);
    }
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class BilinearFormNonsymVel_0_0 : public MatrixFormVol < double >
  {
  public:
    BilinearFormNonsymVel_0_0(unsigned int i, unsigned int j, bool Stokes);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormNonsymVel_0_0(this->i, this->j, this->Stokes);
    }
  protected:
    bool Stokes;
  };

  class BilinearFormNonsymVel_0_1 : public MatrixFormVol < double >
  {
  public:
    BilinearFormNonsymVel_0_1(unsigned int i, unsigned int j, bool Stokes);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormNonsymVel_0_1(this->i, this->j, this->Stokes);
    }
  protected:
    bool Stokes;
  };

  class BilinearFormNonsymVel_1_0 : public MatrixFormVol < double >
  {
  public:
    BilinearFormNonsymVel_1_0(unsigned int i, unsigned int j, bool Stokes);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormNonsymVel_1_0(this->i, this->j, this->Stokes);
    }
  protected:
    bool Stokes;
  };

  class BilinearFormNonsymVel_1_1 : public MatrixFormVol < double >
  {
  public:
    BilinearFormNonsymVel_1_1(unsigned int i, unsigned int j, bool Stokes);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormNonsymVel_1_1(this->i, this->j, this->Stokes);
    }
  protected:
    bool Stokes;
  };

  class BilinearFormNonsymXVelPressure : public MatrixFormVol < double >
  {
  public:
    BilinearFormNonsymXVelPressure(int i, int j);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormNonsymXVelPressure(this->i, this->j);
    }
  };

  class BilinearFormNonsymYVelPressure : public MatrixFormVol < double >
  {
  public:
    BilinearFormNonsymYVelPressure(int i, int j);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
      GeomVol<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e,
      Func<Ord>* *ext) const;
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormNonsymYVelPressure(this->i, this->j);
    }
  };

  class VectorFormNS_0 : public VectorFormVol < double >
  {
  public:
    VectorFormNS_0(int i, bool Stokes, double Reynolds, double time_step);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomVol<double> *e,
      Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const;

    VectorFormVol<double>* clone() const
    {
      return new VectorFormNS_0(this->i, this->Stokes, this->Reynolds, this->time_step);
    }
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_1 : public VectorFormVol < double >
  {
  public:
    VectorFormNS_1(int i, bool Stokes, double Reynolds, double time_step);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomVol<double> *e,
      Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const;

    VectorFormVol<double>* clone() const
    {
      return new VectorFormNS_1(this->i, this->Stokes, this->Reynolds, this->time_step);
    }
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_2 : public VectorFormVol < double >
  {
  public:
    VectorFormNS_2(int i);

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomVol<double> *e,
      Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const;

    VectorFormVol<double>* clone() const
    {
      return new VectorFormNS_2(this->i);
    }
  };

protected:
  bool Stokes;
  double Reynolds;
  double time_step;
  MeshFunctionSharedPtr<double>  x_vel_previous_time;
  MeshFunctionSharedPtr<double>  y_vel_previous_time;
};

class EssentialBCNonConst : public EssentialBoundaryCondition < double >
{
public:
  EssentialBCNonConst(std::vector<std::string> markers, double vel_inlet, double H, double startup_time)
    : EssentialBoundaryCondition<double>(markers), startup_time(startup_time), vel_inlet(vel_inlet), H(H)  {};

  EssentialBCNonConst(std::string marker, double vel_inlet, double H, double startup_time);

  ~EssentialBCNonConst() {};

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y) const;

protected:
  double startup_time;
  double vel_inlet;
  double H;
};
