#include "definitions.h"

/* Weak forms */

class WeakFormRayleighBenard : public WeakForm<double>
{
public:
  WeakFormRayleighBenard(double Pr, double Ra, std::string bdy_top, double temp_ext, double alpha_air,
                         double time_step, Solution<double>* x_vel_previous_time,
                         Solution<double>* y_vel_previous_time, Solution<double>* temp_previous_time)
    : WeakForm<double>(4), Pr(Pr), Ra(Ra), time_step(time_step), x_vel_previous_time(x_vel_previous_time),
                y_vel_previous_time(y_vel_previous_time), temp_previous_time(temp_previous_time) {
    /* Jacobian terms - first velocity equation */
    // Time derivative in the first velocity equation.
    add_matrix_form(new DefaultMatrixFormVol<double>(0, 0, HERMES_ANY, new Hermes2DFunction<double>(1./time_step)));
    // Laplacian divided by Pr in the first velocity equation.
    add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, HERMES_ANY, new Hermes1DFunction<double>(1./Pr)));
    // First part of the convective term in the first velocity equation.
    BilinearFormNonsymVel_0_0* nonsym_vel_form_0_0 = new BilinearFormNonsymVel_0_0(0, 0);
    add_matrix_form(nonsym_vel_form_0_0);
    // Second part of the convective term in the first velocity equation.
    BilinearFormNonsymVel_0_1* nonsym_vel_form_0_1 = new BilinearFormNonsymVel_0_1(0, 1);
    add_matrix_form(nonsym_vel_form_0_1);
    // Pressure term in the first velocity equation.
    BilinearFormNonsymXVelPressure* nonsym_velx_pressure_form = new BilinearFormNonsymXVelPressure(0, 2);
    add_matrix_form(nonsym_velx_pressure_form);

    /* Jacobian terms - second velocity equation, continuity equation */
    // Time derivative in the second velocity equation.
    add_matrix_form(new DefaultMatrixFormVol<double>(1, 1, HERMES_ANY, new Hermes2DFunction<double>(1./time_step)));
    // Laplacian divided by Pr in the second velocity equation.
    add_matrix_form(new DefaultJacobianDiffusion<double>(1, 1, HERMES_ANY, new Hermes1DFunction<double>(1./Pr)));
    // First part of the convective term in the second velocity equation.
    BilinearFormNonsymVel_1_0* nonsym_vel_form_1_0 = new BilinearFormNonsymVel_1_0(1, 0);
    add_matrix_form(nonsym_vel_form_1_0);
    // Second part of the convective term in the second velocity equation.
    BilinearFormNonsymVel_1_1* nonsym_vel_form_1_1 = new BilinearFormNonsymVel_1_1(1, 1);
    add_matrix_form(nonsym_vel_form_1_1);
    // Pressure term in the second velocity equation.
    BilinearFormNonsymYVelPressure* nonsym_vely_pressure_form = new BilinearFormNonsymYVelPressure(1, 2);
    add_matrix_form(nonsym_vely_pressure_form);
    // Temperature term in the second velocity equation.
    add_matrix_form(new DefaultMatrixFormVol<double>(1, 3, HERMES_ANY, new Hermes2DFunction<double>(Ra * Pr)));

    /* Jacobian terms - temperature equation */
    // Time derivative in the temperature equation.
    add_matrix_form(new DefaultMatrixFormVol<double>(3, 3, HERMES_ANY, new Hermes2DFunction<double>(1./time_step)));
    // Laplacian in the temperature equation.
    add_matrix_form(new DefaultJacobianDiffusion<double>(3, 3));
    // First part of temperature advection term.
    add_matrix_form(new BilinearFormNonsymTemp_3_0(3, 0));
    // Second part of temperature advection term.
    add_matrix_form(new BilinearFormNonsymTemp_3_1(3, 1));
    // Third part of temperature advection term.
    add_matrix_form(new BilinearFormNonsymTemp_3_3(3, 3));
    // Surface term generated by the Newton condition on top edge.
    add_matrix_form_surf(new DefaultMatrixFormSurf<double>(3, 3, bdy_top, new Hermes2DFunction<double>(alpha_air)));

    /* Residual - volumetric */
    // First velocity equation.
    VectorFormNS_0* F_0 = new VectorFormNS_0(0, Pr, time_step);
    F_0->ext.push_back(x_vel_previous_time);
    add_vector_form(F_0);
    // Second velocity equation.
    VectorFormNS_1* F_1 = new VectorFormNS_1(1, Pr, Ra, time_step);
    F_1->ext.push_back(y_vel_previous_time);
    add_vector_form(F_1);
    // Continuity equation.
    VectorFormNS_2* F_2 = new VectorFormNS_2(2);
    add_vector_form(F_2);
    // Temperature equation.
    VectorFormNS_3* F_3 = new VectorFormNS_3(3, time_step);
    F_3->ext.push_back(temp_previous_time);
    add_vector_form(F_3);
    add_vector_form_surf(new DefaultVectorFormSurf<double>(3, bdy_top, new Hermes2DFunction<double>(-alpha_air * temp_ext)));
    add_vector_form_surf(new CustomResidualSurfConst(3, bdy_top, alpha_air));
  };

  class BilinearFormNonsymVel_0_0 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymVel_0_0(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM) {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const {
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i]
         + u->val[i] * xvel_prev_newton->dx[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
    Ord result = Ord(0);
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i]
         + u->val[i] * xvel_prev_newton->dx[i] * v->val[i]);
      return result;
    }
  };

  class BilinearFormNonsymVel_0_1 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymVel_0_1(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const {
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * xvel_prev_newton->dy[i] * v->val[i] ;

      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
      Ord result = Ord(0);

      Func<Ord>* xvel_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] * xvel_prev_newton->dy[i]);

      return result;
    }
  };
  class BilinearFormNonsymVel_1_0 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymVel_1_0(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const {
      double result = 0;
      Func<double>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * yvel_prev_newton->dx[i] * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * yvel_prev_newton->dx[i] * v->val[i];
      return result;
    }
  };
  class BilinearFormNonsymVel_1_1 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymVel_1_1(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const {
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i]
                           + u->val[i] * yvel_prev_newton->dy[i] * v->val[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] * u->dx[i] + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i]
                           + u->val[i] * yvel_prev_newton->dy[i] * v->val[i]);
      return result;
    }
  };

  class BilinearFormNonsymXVelPressure : public MatrixFormVol<double>
  {
  public:
    // The antisym flag is used here to generate a term in the continuity equation.
    BilinearFormNonsymXVelPressure(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_ANTISYM) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const {
      return - int_u_dvdx<double, double>(n, wt, u, v);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
      return - int_u_dvdx<Ord, Ord>(n, wt, u, v);
    }
  };
  class BilinearFormNonsymYVelPressure : public MatrixFormVol<double>
  {
  public:
    // The antisym flag is used here to generate a term in the continuity equation.
    BilinearFormNonsymYVelPressure(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_ANTISYM) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const {
      return - int_u_dvdy<double, double>(n, wt, u, v);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
      return - int_u_dvdy<Ord, Ord>(n, wt, u, v);
    }
  };

  class VectorFormNS_0 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_0(int i, double Pr, double time_step) : VectorFormVol<double>(i), Pr(Pr), time_step(time_step) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<double> *ext) const {
      double result = 0;
      Func<double>* xvel_prev_time = ext->fn[0];
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];
      Func<double>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Pr
                          - (p_prev_newton->val[i] * v->dx[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step
                          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i]
                          + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_time = ext->fn[0];
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      Func<Ord>* p_prev_newton = u_ext[2];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Pr
                          - (p_prev_newton->val[i] * v->dx[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->val[i] - xvel_prev_time->val[i]) * v->val[i] / time_step
                          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i]
                          + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i]));
      return result;
    }
  protected:
    double Pr, time_step;
  };

  class VectorFormNS_1 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_1(int i, double Pr, double Ra, double time_step)
      : VectorFormVol<double>(i), Pr(Pr), Ra(Ra), time_step(time_step) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<double> *ext) const {
      double result = 0;
      Func<double>* yvel_prev_time = ext->fn[0];
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];
      Func<double>* p_prev_newton = u_ext[2];
      Func<double>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((yvel_prev_newton->dx[i] * v->dx[i] + yvel_prev_newton->dy[i] * v->dy[i]) / Pr
                          - (p_prev_newton->val[i] * v->dy[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * ((yvel_prev_newton->val[i] - yvel_prev_time->val[i]) * v->val[i] / time_step
                          + ((xvel_prev_newton->val[i] * yvel_prev_newton->dx[i]
            + yvel_prev_newton->val[i] * yvel_prev_newton->dy[i]) * v->val[i])
         + Pr*Ra*temp_prev_newton->val[i]*v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* yvel_prev_time = ext->fn[0];
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      Func<Ord>* p_prev_newton = u_ext[2];
      Func<Ord>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * ((xvel_prev_newton->dx[i] * v->dx[i] + xvel_prev_newton->dy[i] * v->dy[i]) / Pr
                  - (p_prev_newton->val[i] * v->dx[i]));
      for (int i = 0; i < n; i++)
        result += wt[i] * ((yvel_prev_newton->val[i] - yvel_prev_time->val[i]) * v->val[i] / time_step
                          + ((xvel_prev_newton->val[i] * xvel_prev_newton->dx[i]
                          + yvel_prev_newton->val[i] * xvel_prev_newton->dy[i]) * v->val[i])
                           + Pr*Ra*temp_prev_newton->val[i]*v->val[i]);
      return result;
    }
  protected:
    double Pr, Ra, time_step;
  };

  class VectorFormNS_2 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_2(int i) : VectorFormVol<double>(i) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<double> *ext) const {
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];

      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->dx[i] + yvel_prev_newton->dy[i]) * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];

      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->dx[i] + yvel_prev_newton->dy[i]) * v->val[i];
      return result;
    }
  };

  class VectorFormNS_3 : public VectorFormVol<double>
  {
  public:
    VectorFormNS_3(int i, double time_step) : VectorFormVol<double>(i), time_step(time_step) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
                         ExtData<double> *ext) const {
      double result = 0;
      Func<double>* temp_prev_time = ext->fn[0];
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];
      Func<double>* temp_prev_newton = u_ext[3];

      for (int i = 0; i < n; i++)
        result += wt[i] * (((temp_prev_newton->val[i] - temp_prev_time->val[i]) / time_step
                           + xvel_prev_newton->val[i] * temp_prev_newton->dx[i]
         + yvel_prev_newton->val[i] * temp_prev_newton->dy[i]) * v->val[i]
                     + temp_prev_newton->dx[i] * v->dx[i] + temp_prev_newton->dy[i] * v->dy[i]
         );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* temp_prev_time = ext->fn[0];
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      Func<Ord>* temp_prev_newton = u_ext[3];

      for (int i = 0; i < n; i++)
        result += wt[i] * (((temp_prev_newton->val[i] - temp_prev_time->val[i]) / time_step
                           + xvel_prev_newton->val[i] * temp_prev_newton->dx[i]
         + yvel_prev_newton->val[i] * temp_prev_newton->dy[i]) * v->val[i]
                     + temp_prev_newton->dx[i] * v->dx[i] + temp_prev_newton->dy[i] * v->dy[i]
         );
      return result;
    }

    private:
      double time_step;
  };

  class BilinearFormNonsymTemp_3_0 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymTemp_3_0(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const {
      double result = 0;
      Func<double>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * temp_prev_newton->dx[i] * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * temp_prev_newton->dx[i] * v->val[i];
      return result;
    }
  };

  class BilinearFormNonsymTemp_3_1 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymTemp_3_1(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const {
      double result = 0;
      Func<double>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * temp_prev_newton->dy[i] * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* temp_prev_newton = u_ext[3];
      for (int i = 0; i < n; i++)
        result += wt[i] * u->val[i] * temp_prev_newton->dy[i] * v->val[i];
      return result;
    }
  };

  class BilinearFormNonsymTemp_3_3 : public MatrixFormVol<double>
  {
  public:
    BilinearFormNonsymTemp_3_3(int i, int j) : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM) {
      
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const {
      double result = 0;
      Func<double>* xvel_prev_newton = u_ext[0];
      Func<double>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->val[i] * u->dx[i]
                           + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i];
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
      Ord result = Ord(0);
      Func<Ord>* xvel_prev_newton = u_ext[0];
      Func<Ord>* yvel_prev_newton = u_ext[1];
      for (int i = 0; i < n; i++)
        result += wt[i] * (xvel_prev_newton->val[i] * u->dx[i]
                           + yvel_prev_newton->val[i] * u->dy[i]) * v->val[i];
      return result;
    }
  };

  class CustomResidualSurfConst : public VectorFormSurf<double>
  {
  public:
    CustomResidualSurfConst(int i, double coeff = 1.0,
                             GeomType gt = HERMES_PLANAR)
           : VectorFormSurf<double>(i), coeff(coeff), gt(gt) { }
    CustomResidualSurfConst(int i, std::string area, double coeff = 1.0,
                             GeomType gt = HERMES_PLANAR)
           : VectorFormSurf<double>(i, area), coeff(coeff), gt(gt) { }

    template<typename Real, typename Scalar>
    Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = Scalar(0);
      for (int i = 0; i < n; i++) {
        result += wt[i] * u_ext[3]->val[i] * v->val[i];
      }
      return coeff * result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<double> *ext) const {
      return vector_form_surf<double, double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // This is to make the form usable in rk_time_step_newton().
    virtual VectorFormSurf* clone() {
      return new CustomResidualSurfConst(*this);
    }

    private:
      double coeff;
      GeomType gt;
  };
protected:
  double Pr, Ra, time_step;
  Solution<double>* x_vel_previous_time;
  Solution<double>* y_vel_previous_time;
  Solution<double>* temp_previous_time;
};
