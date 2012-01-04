#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;
using namespace Hermes::Hermes2D::WeakFormsH1;

class WeakFormNSSimpleLinearization : public WeakForm<double>
{
public:
  WeakFormNSSimpleLinearization(bool Stokes, double Reynolds, double time_step, Solution<double>* x_vel_previous_time, 
                                Solution<double>* y_vel_previous_time);

  class BilinearFormSymVel : public MatrixFormVol<double>
  {
  public:
    BilinearFormSymVel(int i, int j, bool Stokes, double Reynolds, double time_step) 
            : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_SYM), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext);

    MatrixFormVol<double>* clone();
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class BilinearFormNonsymVel : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymVel(int i, int j, bool Stokes) 
            : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;
    
    MatrixFormVol<double>* clone();
  protected:
    bool Stokes;
  };

  class BilinearFormNonsymXVelPressure : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymXVelPressure(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_ANTISYM) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;
    
    MatrixFormVol<double>* clone();
  };

  class BilinearFormNonsymYVelPressure : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymYVelPressure(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_ANTISYM) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;
    
    MatrixFormVol<double>* clone();
  };

  class VectorFormVolVel : public VectorFormVol<double>
  {
  public:
    VectorFormVolVel(int i, bool Stokes, double time_step) 
          : VectorFormVol<double>(i), Stokes(Stokes), time_step(time_step) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    VectorFormVol<double>* clone();
  protected:
    bool Stokes;
    double time_step;
  };

protected:
  bool Stokes;
  double Reynolds;
  double time_step;
  Solution<double>* x_vel_previous_time;
  Solution<double>* y_vel_previous_time;
};

class WeakFormNSNewton : public WeakForm<double>
{
public:
  WeakFormNSNewton(bool Stokes, double Reynolds, double time_step, Solution<double>* x_vel_previous_time, 
                   Solution<double>* y_vel_previous_time);

  class BilinearFormSymVel : public MatrixFormVol<double>
  {
  public:
    BilinearFormSymVel(int i, int j, bool Stokes, double Reynolds, double time_step) 
            : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_SYM), Stokes(Stokes), 
                        Reynolds(Reynolds), time_step(time_step) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    MatrixFormVol<double>* clone();
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class BilinearFormNonsymVel_0_0 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymVel_0_0(int i, int j, bool Stokes) 
            : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    MatrixFormVol<double>* clone();
  protected:
    bool Stokes;
  };

  class BilinearFormNonsymVel_0_1 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymVel_0_1(int i, int j, bool Stokes) 
            : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    MatrixFormVol<double>* clone();
  protected:
    bool Stokes;
  };

  class BilinearFormNonsymVel_1_0 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymVel_1_0(int i, int j, bool Stokes) 
            : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    MatrixFormVol<double>* clone();
  protected:
    bool Stokes;
  };

  class BilinearFormNonsymVel_1_1 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymVel_1_1(int i, int j, bool Stokes) 
            : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), Stokes(Stokes) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;

    MatrixFormVol<double>* clone();
  protected:
    bool Stokes;
  };

  class BilinearFormNonsymXVelPressure : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymXVelPressure(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_ANTISYM) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;
    MatrixFormVol<double>* clone();
  };

  class BilinearFormNonsymYVelPressure : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymYVelPressure(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_ANTISYM) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                    ExtData<Ord> *ext) const;
    MatrixFormVol<double>* clone();
  };

  class VectorFormNS_0 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_0(int i, bool Stokes, double Reynolds, double time_step) : VectorFormVol<double>(i), 
                   Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    VectorFormVol<double>* clone();
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_1 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_1(int i, bool Stokes, double Reynolds, double time_step) 
          : VectorFormVol<double>(i), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;

    VectorFormVol<double>* clone();
  protected:
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_2 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_2(int i) : VectorFormVol<double>(i) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, 
                         ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
    
    VectorFormVol<double>* clone();
  };

protected:
  bool Stokes;
  double Reynolds;
  double time_step;
  Solution<double>* x_vel_previous_time;
  Solution<double>* y_vel_previous_time;
};

class EssentialBCNonConst : public EssentialBoundaryCondition<double>
{
public:
  EssentialBCNonConst(Hermes::vector<std::string> markers, double vel_inlet, double H, double startup_time) 
             : EssentialBoundaryCondition<double>(markers), startup_time(startup_time), vel_inlet(vel_inlet), H(H)  {};

  EssentialBCNonConst(std::string marker, double vel_inlet, double H, double startup_time);
  
  ~EssentialBCNonConst() {};

  virtual EssentialBoundaryCondition<double>::EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

protected:
  double startup_time;
  double vel_inlet;
  double H;
};