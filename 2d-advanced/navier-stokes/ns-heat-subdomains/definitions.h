#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* These numbers must be compatible with mesh file */

// These numbers must be compatible with mesh file.
const double H = 1.0;                               
const double OBSTACLE_DIAMETER = 0.3 * 1.4142136;    
const double HOLE_MID_X = 0.5;
const double HOLE_MID_Y = 0.5;

/* Custom initial condition for temperature*/

class CustomInitialConditionTemperature : public ExactSolutionScalar<double>
{
public:
  CustomInitialConditionTemperature(MeshSharedPtr mesh, double mid_x, double mid_y, double radius, double temp_fluid, double temp_graphite);

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(double x, double y) const;

  MeshFunction<double>* clone() const;

  // Members.
  double mid_x, mid_y, radius, temp_fluid, temp_graphite;
};

/* Weak forms */

class CustomWeakFormHeatAndFlow : public WeakForm<double>
{
public:
  CustomWeakFormHeatAndFlow(bool Stokes, double Reynolds, double time_step, MeshFunctionSharedPtr<double> x_vel_previous_time, 
    MeshFunctionSharedPtr<double> y_vel_previous_time, MeshFunctionSharedPtr<double> T_prev_time, double heat_source, double specific_heat_graphite, 
    double specific_heat_fluid, double rho_graphite, double rho_fluid, double thermal_conductivity_graphite, 
    double thermal_conductivity_fluid, bool simple_temp_advection);

  class BilinearFormTime: public MatrixFormVol<double>
  {
  public:
    BilinearFormTime(int i, int j, std::string area, double time_step) : MatrixFormVol<double>(i, j), time_step(time_step) {
      this->set_area(area);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const
    {
      double result = int_u_v<double, double>(n, wt, u, v) / time_step;
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
    {
      Ord result = int_u_v<Ord, Ord>(n, wt, u, v) / time_step;
      return result;
    }
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormTime(*this);
    }
  protected:
    // Members.
    double time_step;
  };

  class BilinearFormSymVel : public MatrixFormVol<double>
  {
  public:
    BilinearFormSymVel(int i, int j, bool Stokes, double Reynolds, double time_step) : MatrixFormVol<double>(i, j), Stokes(Stokes), 
      Reynolds(Reynolds), time_step(time_step) {
      this->setSymFlag(sym);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      double result = int_grad_u_grad_v<double, double>(n, wt, u, v) / Reynolds;
      if(!Stokes)
        result += int_u_v<double, double>(n, wt, u, v) / time_step;
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
    {
      Ord result = int_grad_u_grad_v<Ord, Ord>(n, wt, u, v) / Reynolds;
      if(!Stokes)
        result += int_u_v<Ord, Ord>(n, wt, u, v) / time_step;
      return result;
    }
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormSymVel(*this);
    }
  protected:
    // Members.
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class BilinearFormUnSymVel_0_0 : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymVel_0_0(int i, int j, bool Stokes) : MatrixFormVol<double>(i, j), Stokes(Stokes) {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      double result = 0;
      if(!Stokes) {
        Func<double>* xvel_prev_newton = u_ext[0];
        Func<double>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i]
        * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * xvel_prev_newton->dx[i]);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const{
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* xvel_prev_newton = u_ext[0];
        Func<Ord>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i]
        * u->dy[i]) * v->val[i] + u->val[i] * v->val[i] * xvel_prev_newton->dx[i]);
      }
      return result;
    }
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormUnSymVel_0_0(*this);
    }
  protected:
    // Members.
    bool Stokes;
  };

  class BilinearFormUnSymVel_0_1 : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymVel_0_1(int i, int j, bool Stokes) : MatrixFormVol<double>(i, j), Stokes(Stokes) {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      double result = 0;
      if(!Stokes) {
        Func<double>* xvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* xvel_prev_newton = u_ext[0];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);
      }
      return result;
    }
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormUnSymVel_0_1(*this);
    }
  protected:
    // Members.
    bool Stokes;
  };

  class BilinearFormUnSymVel_1_0 : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymVel_1_0(int i, int j, bool Stokes) : MatrixFormVol<double>(i, j), Stokes(Stokes) {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      double result = 0;
      if(!Stokes) {
        Func<double>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const{
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * (u->val[i] * v->val[i] * yvel_prev_newton->dx[i]);
      }
      return result;
    }
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormUnSymVel_1_0(*this);
    }
  protected:
    // Members.
    bool Stokes;
  };

  class BilinearFormUnSymVel_1_1 : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymVel_1_1(int i, int j, bool Stokes) : MatrixFormVol<double>(i, j), Stokes(Stokes) {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      double result = 0;
      if(!Stokes) {
        Func<double>* xvel_prev_newton = u_ext[0];
        Func<double>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] + u->val[i]
        * v->val[i] * yvel_prev_newton->dy[i]);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      Ord result = Ord(0);
      if(!Stokes) {
        Func<Ord>* xvel_prev_newton = u_ext[0];
        Func<Ord>* yvel_prev_newton = u_ext[1];
        for (int i = 0; i < n; i++)
          result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i] + u->val[i]
        * v->val[i] * yvel_prev_newton->dy[i]);
      }
      return result;
    }
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormUnSymVel_1_1(*this);
    }
  protected:
    // Members.
    bool Stokes;
  };

  class BilinearFormUnSymXVelPressure : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymXVelPressure(int i, int j) : MatrixFormVol<double>(i, j) {
      this->setSymFlag(HERMES_ANTISYM);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      return - int_u_dvdx<double, double>(n, wt, u, v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      return - int_u_dvdx<Ord, Ord>(n, wt, u, v);
    }
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormUnSymXVelPressure(*this);
    }
  };

  class BilinearFormUnSymYVelPressure : public MatrixFormVol<double>
  {
  public:
    BilinearFormUnSymYVelPressure(int i, int j) : MatrixFormVol<double>(i, j) {
      this->setSymFlag(HERMES_ANTISYM);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      return - int_u_dvdy<double, double>(n, wt, u, v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      return - int_u_dvdy<Ord, Ord>(n, wt, u, v);
    }
    MatrixFormVol<double>* clone() const
    {
      return new BilinearFormUnSymYVelPressure(*this);
    }
  };

  class CustomJacobianTempAdvection_3_0 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianTempAdvection_3_0(int i, int j, std::string area) : MatrixFormVol<double>(i, j) 
    {
      this->set_area(area);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* T_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++) 
      {
        result += wt[i] * u->val[i] * T_prev_newton->dx[i] * v->val[i];
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
    {
      Ord result = Ord(0);
      Func<Ord>* T_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++) 
      {
        result += wt[i] * u->val[i] * T_prev_newton->dx[i] * v->val[i];
      }
      return result;
    }
    MatrixFormVol<double>* clone() const
    {
      return new CustomJacobianTempAdvection_3_0(*this);
    }
  };

  class CustomJacobianTempAdvection_3_3_simple : public MatrixFormVol<double>
  {
  public:
    CustomJacobianTempAdvection_3_3_simple(int i, int j, std::string area) : MatrixFormVol<double>(i, j) 
    {
      this->set_area(area);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* xvel_prev_time = ext[0];
      Func<double>* yvel_prev_time = ext[1];
      for (int i = 0; i < n; i++) 
      {
        result += wt[i] * (xvel_prev_time->val[i] * u->dx[i] + yvel_prev_time->val[i] * u->dy[i]) * v->val[i];
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
    {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_time = ext[0];
      Func<Ord>* yvel_prev_time = ext[1];
      for (int i = 0; i < n; i++) 
      {
        result += wt[i] * (xvel_prev_time->val[i] * u->dx[i] + yvel_prev_time->val[i] * u->dy[i]) * v->val[i];
      }
      return result;
    }
    MatrixFormVol<double>* clone() const
    {
      return new CustomJacobianTempAdvection_3_3_simple(*this);
    }
  };

  class CustomJacobianTempAdvection_3_1 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianTempAdvection_3_1(int i, int j, std::string area) : MatrixFormVol<double>(i, j) 
    {
      this->set_area(area);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      double result = 0;
      Func<double>* T_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++) 
      {
        result += wt[i] * u->val[i] * T_prev_newton->dy[i] * v->val[i];
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
    {
      Ord result = Ord(0);
      Func<Ord>* T_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
      {
        result += wt[i] * u->val[i] * T_prev_newton->dy[i] * v->val[i];
      }
      return result;
    }
    MatrixFormVol<double>* clone() const
    {
      return new CustomJacobianTempAdvection_3_1(*this);
    }
  };

  class CustomJacobianTempAdvection_3_3 : public MatrixFormVol<double>
  {
  public:
    CustomJacobianTempAdvection_3_3(int i, int j, std::string area) : MatrixFormVol<double>(i, j) 
    {
      this->set_area(area);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++) 
      {
        result += wt[i] * (  xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i] ) * v->val[i];
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
    {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
      {
        result += wt[i] * (  xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i] ) * v->val[i];
      }
      return result;
    }
    MatrixFormVol<double>* clone() const
    {
      return new CustomJacobianTempAdvection_3_3(*this);
    }
  };

  class VectorFormTime: public VectorFormVol<double>
  {
  public:
    VectorFormTime(int i, std::string area, double time_step) : VectorFormVol<double>(i), time_step(time_step) 
    {
      this->set_area(area);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const
    {
      Func<double>* func_prev_time = ext[0];
      double result = (int_u_v<double, double>(n, wt, u_ext[3], v) - int_u_v<double, double>(n, wt, func_prev_time, v)) / time_step;
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
    {
      Func<Ord>* func_prev_time = ext[0];
      Ord result = (int_u_v<Ord, Ord>(n, wt, u_ext[3], v) - int_u_v<Ord, Ord>(n, wt, func_prev_time, v)) / time_step;
      return result;
    }
    VectorFormVol<double>* clone() const
    {
      return new VectorFormTime(*this);
    }
  protected:
    // Members.
    double time_step;
  };

  class CustomResidualTempAdvection : public VectorFormVol<double>
  {
  public:
    CustomResidualTempAdvection(int i, std::string area) : VectorFormVol<double>(i) 
    {
      this->set_area(area);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];
      Func<double>* T_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++) 
      {
        result += wt[i] * ( xvel_prev_newton->val[i] * T_prev_newton->dx[i] + yvel_prev_newton->val[i] * T_prev_newton->dy[i] ) * v->val[i]; 
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
    {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      Func<Ord>* T_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++) 
      {
        result += wt[i] * ( xvel_prev_newton->val[i] * T_prev_newton->dx[i] + yvel_prev_newton->val[i] * T_prev_newton->dy[i] ) * v->val[i]; 
      }
      return result;
    }
    VectorFormVol<double>* clone() const
    {
      return new CustomResidualTempAdvection(*this);
    }
  };

  class CustomResidualTempAdvection_simple : public VectorFormVol<double>
  {
  public:
    CustomResidualTempAdvection_simple(int i, std::string area) : VectorFormVol<double>(i) 
    {
      this->set_area(area);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* xvel_prev_time = ext[0];
      Func<double>* yvel_prev_time = ext[1];
      Func<double>* T_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) 
      {
        result += wt[i] * ( xvel_prev_time->val[i] * T_prev_newton->dx[i] + yvel_prev_time->val[i] * T_prev_newton->dy[i] ) * v->val[i]; 
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
    {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_time = ext[0];
      Func<Ord>* yvel_prev_time = ext[1];
      Func<Ord>* T_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) 
      {
        result += wt[i] * ( xvel_prev_time->val[i] * T_prev_newton->dx[i] + yvel_prev_time->val[i] * T_prev_newton->dy[i] ) * v->val[i]; 
      }
      return result;
    }
    VectorFormVol<double>* clone() const
    {
      return new CustomResidualTempAdvection_simple(*this);
    }
  };

  class VectorFormNS_0 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_0(int i, bool Stokes, double Reynolds, double time_step) : VectorFormVol<double>(i), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) 
    {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      double result = 0;
      Func<double>* xvel_prev_time = ext[0];
      Func<double>* yvel_prev_time = ext[1];
      Func<double>* xvel_prev_newton = u_ext[0];  
      Func<double>* yvel_prev_newton = u_ext[1];  
      Func<double>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Reynolds - (p_prev_newton->val[i] * v->dx[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step )
          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_time = ext[0];  
      Func<Ord>* yvel_prev_time = ext[1];
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  
      Func<Ord>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Reynolds - (p_prev_newton->val[i] * v->dx[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step)
          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
    VectorFormVol<double>* clone() const
    {
      return new VectorFormNS_0(*this);
    }
  protected:
    // Members.
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_1 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_1(int i, bool Stokes, double Reynolds, double time_step) : VectorFormVol<double>(i), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step) {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      double result = 0;
      Func<double>* xvel_prev_time = ext[0];  
      Func<double>* yvel_prev_time = ext[1];
      Func<double>* xvel_prev_newton = u_ext[0];  
      Func<double>* yvel_prev_newton = u_ext[1];  
      Func<double>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((yvel_prev_newton->dx[i] * v->dx[i] + yvel_prev_newton->dy[i] * v->dy[i]) / Reynolds - (p_prev_newton->val[i] * v->dy[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((yvel_prev_newton->val[i] - yvel_prev_time->val[i]) * v->val[i] / time_step )
          + ((xvel_prev_newton->val[i] * yvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * yvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const{
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_time = ext[0];  
      Func<Ord>* yvel_prev_time = ext[1];
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  
      Func<Ord>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Reynolds - (p_prev_newton->val[i] * v->dx[i]));
      if(!Stokes)
        for (int i = 0; i < n; i++)
          result += wt[i] * (((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step )
          + ((xvel_prev_newton->val[i] * yvel_prev_newton->dx[i] + yvel_prev_newton->val[i] * yvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
    VectorFormVol<double>* clone() const
    {
      return new VectorFormNS_1(*this);
    }
  protected:
    // Members.
    bool Stokes;
    double Reynolds;
    double time_step;
  };

  class VectorFormNS_2 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_2(int i) : VectorFormVol<double>(i) {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const{
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];  
      Func<double>* yvel_prev_newton = u_ext[1];  

      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->dx[i] * v->val[i] + yvel_prev_newton->dy[i] * v->val[i]);
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_newton = u_ext[0];  
      Func<Ord>* yvel_prev_newton = u_ext[1];  

      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->dx[i] * v->val[i] + yvel_prev_newton->dy[i] * v->val[i]);
      return result;
    }
    VectorFormVol<double>* clone() const
    {
      return new VectorFormNS_2(*this);
    }
  };

protected:
  // Members.
  bool Stokes;
  double Reynolds;
  double time_step;
  MeshFunctionSharedPtr<double> x_vel_previous_time;
  MeshFunctionSharedPtr<double> y_vel_previous_time;
};

class EssentialBCNonConst : public EssentialBoundaryCondition<double>
{
public:
  EssentialBCNonConst(Hermes::vector<std::string> markers, double vel_inlet, double H, double startup_time) : 
      EssentialBoundaryCondition<double>(markers), vel_inlet(vel_inlet), H(H), startup_time(startup_time) {};
      EssentialBCNonConst(std::string marker, double vel_inlet, double H, double startup_time) : 
      EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), vel_inlet(vel_inlet), H(H), startup_time(startup_time) {
        markers.push_back(marker);
      };

      ~EssentialBCNonConst() {};

      virtual EssentialBCValueType get_value_type() const { 
        return BC_FUNCTION; 
      };

      virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
        double val_y = vel_inlet * y*(H-y) / (H/2.)/(H/2.);  // Parabolic profile.
        //double val_y = vel_inlet;                            // Constant profile.
        if (current_time <= startup_time) 
          return val_y * current_time/startup_time;
        else 
          return val_y;
      };

protected:
  // Members.
  double vel_inlet;
  double H;
  double startup_time;
};

bool point_in_graphite(double x, double y);
int element_in_graphite(Element* e);
int element_in_fluid(Element* e);
