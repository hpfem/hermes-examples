#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Global variables for forms.
double K_S, ALPHA, THETA_R, THETA_S, N, M;

// Problem parameters.
const double TAU = 5e-3;                          // Time step.
const double STARTUP_TIME = 1.1e-2;               // Start-up time for time-dependent Dirichlet boundary condition.
const double T_FINAL = 5.0;                       // Time interval length.
double TIME = 0;                                  // Global time variable initialized with first time step.
double H_INIT = -9.5;                             // Initial pressure head.
double H_ELEVATION = 5.2;

double K_S_1 = 0.108;
double K_S_3 = 0.0048;
double K_S_2 = 0.0168;
double K_S_4 = 1.061;

double ALPHA_1 = 0.01;
double ALPHA_3 = 0.005;
double ALPHA_2 = 0.01;
double ALPHA_4 = 0.05;

double THETA_R_1 = 0.1020;
double THETA_R_2 = 0.09849;
double THETA_R_3 = 0.08590;
double THETA_R_4 = 0.08590;

double THETA_S_1 = 0.4570;
double THETA_S_2 = 0.4510;
double THETA_S_3 = 0.4650;
double THETA_S_4 = 0.5650;

double N_1 = 1.982;
double N_2 = 1.632; 
double N_3 = 5.0;
double N_4 = 5.0;

double M_1 = 0.49546;
double M_2 = 0.38726;
double M_3 = 0.8;
double M_4 = 0.8;

// Boundary markers.
const std::string BDY_1 = "1";
const std::string BDY_2 = "2";
const std::string BDY_3 = "3";
const std::string BDY_4 = "4";
const std::string BDY_5 = "5";
const std::string BDY_6 = "6";

class CustomWeakForm : public WeakForm<double>
{
  CustomWeakForm(int TIME_INTEGRATION, MeshFunctionSharedPtr<double> prev_sln)
  {
    this->set_ext(prev_sln);
    if (TIME_INTEGRATION == 1) {
      this->add_matrix_form(new JacobianFormVolEuler(0, 0));
      this->add_matrix_form_surf(new JacobianFormSurf1Euler(0, 0));
      this->mfsurf.back()->set_area(BDY_1);
      this->add_matrix_form_surf(new JacobianFormSurf4Euler(0, 0));
      this->mfsurf.back()->set_area(BDY_4);
      this->add_matrix_form_surf(new JacobianFormSurf6Euler(0, 0));
      this->mfsurf.back()->set_area(BDY_6);

      this->add_vector_form(new ResidualFormVolEuler(0));
      this->add_vector_form_surf(new ResidualFormSurf1Euler(0));
      this->vfsurf.back()->set_area(BDY_1);
      this->add_vector_form_surf(new ResidualFormSurf4Euler(0));
      this->vfsurf.back()->set_area(BDY_4);
      this->add_vector_form_surf(new ResidualFormSurf6Euler(0));
      this->vfsurf.back()->set_area(BDY_6);
    }
    else
    {
      this->add_matrix_form(new JacobianFormVolCrankNicolson(0, 0));
      this->add_matrix_form_surf(new JacobianFormSurf1CrankNicolson(0, 0));
      this->mfsurf.back()->set_area(BDY_1);
      this->add_matrix_form_surf(new JacobianFormSurf4CrankNicolson(0, 0));
      this->mfsurf.back()->set_area(BDY_4);
      this->add_matrix_form_surf(new JacobianFormSurf6CrankNicolson(0, 0));
      this->mfsurf.back()->set_area(BDY_6);

      this->add_vector_form(new ResidualFormVolCrankNicolson(0));
      this->add_vector_form_surf(new ResidualFormSurf1CrankNicolson(0));
      this->vfsurf.back()->set_area(BDY_1);
      this->add_vector_form_surf(new ResidualFormSurf4CrankNicolson(0));
      this->vfsurf.back()->set_area(BDY_4);
      this->add_vector_form_surf(new ResidualFormSurf6CrankNicolson(0));
      this->vfsurf.back()->set_area(BDY_6);
    }
  }

  class JacobianFormVolEuler : public MatrixFormVol<double>
  {
  public:
    JacobianFormVolEuler(int i, int j) : MatrixFormVol<double>(i, j)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const
    {
      double x = e->x[0];
      double y = e->x[1];
      if (is_in_mat_1(x,y)) {
        K_S = K_S_1;
        ALPHA = ALPHA_1;
        THETA_R = THETA_R_1;
        THETA_S = THETA_R_1;
        N = N_1;
        M = M_1;
      }
      if (is_in_mat_2(x,y)) {
        K_S = K_S_2;
        ALPHA = ALPHA_2;
        THETA_R = THETA_R_2;
        THETA_S = THETA_R_2;
        N = N_2;
        M = M_2;
      }
      if (is_in_mat_3(x,y)) {
        K_S = K_S_3;
        ALPHA = ALPHA_3;
        THETA_R = THETA_R_3;
        THETA_S = THETA_R_3;
        N = N_3;
        M = M_3;
      }
      if (is_in_mat_4(x,y)) {
        K_S = K_S_4;
        ALPHA = ALPHA_4;
        THETA_R = THETA_R_4;
        THETA_S = THETA_R_4;
        N = N_4;
        M = M_4;
      }

      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * (
        C(h_prev_newton->val[i]) * u->val[i] * v->val[i] / TAU
        + dCdh(h_prev_newton->val[i]) * u->val[i] * h_prev_newton->val[i] * v->val[i] / TAU
        - dCdh(h_prev_newton->val[i]) * u->val[i] * h_prev_time->val[i] * v->val[i] / TAU
        + K(h_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
        + dKdh(h_prev_newton->val[i]) * u->val[i] * 
        (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
        - dKdh(h_prev_newton->val[i]) * u->dy[i] * v->val[i]
      - ddKdhh(h_prev_newton->val[i]) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
      );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomVol<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    MatrixFormVol<double>* clone() const
    {
      return new JacobianFormVolEuler();
    }
  };

  class JacobianFormVolCrankNicolson : public MatrixFormVol<double>
  {
  public:
    JacobianFormVolCrankNicolson(int i, int j) : MatrixFormVol<double>(i, j)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const
    {
      double x = e->x[0];
      double y = e->x[1];
      if (is_in_mat_1(x,y)) {
        K_S = K_S_1;
        ALPHA = ALPHA_1;
        THETA_R = THETA_R_1;
        THETA_S = THETA_R_1;
        N = N_1;
        M = M_1;
      }
      if (is_in_mat_2(x,y)) {
        K_S = K_S_2;
        ALPHA = ALPHA_2;
        THETA_R = THETA_R_2;
        THETA_S = THETA_R_2;
        N = N_2;
        M = M_2;
      }
      if (is_in_mat_3(x,y)) {
        K_S = K_S_3;
        ALPHA = ALPHA_3;
        THETA_R = THETA_R_3;
        THETA_S = THETA_R_3;
        N = N_3;
        M = M_3;
      }
      if (is_in_mat_4(x,y)) {
        K_S = K_S_4;
        ALPHA = ALPHA_4;
        THETA_R = THETA_R_4;
        THETA_S = THETA_R_4;
        N = N_4;
        M = M_4;
      }

      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * 0.5 * ( // implicit Euler part:
        C(h_prev_newton->val[i]) * u->val[i] * v->val[i] / TAU
        + dCdh(h_prev_newton->val[i]) * u->val[i] * h_prev_newton->val[i] * v->val[i] / TAU
        - dCdh(h_prev_newton->val[i]) * u->val[i] * h_prev_time->val[i] * v->val[i] / TAU
        + K(h_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
        + dKdh(h_prev_newton->val[i]) * u->val[i] * 
        (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
        - dKdh(h_prev_newton->val[i]) * u->dy[i] * v->val[i]
      - ddKdhh(h_prev_newton->val[i]) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
      )
        + wt[i] * 0.5 * ( // explicit Euler part, 
        C(h_prev_time->val[i]) * u->val[i] * v->val[i] / TAU
        );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomVol<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    MatrixFormVol<double>* clone() const
    {
      return new JacobianFormVolEuler(this->i, this->j);
    }
  };

  class JacobianFormSurf1Euler : public MatrixFormSurf<double>
  {
  public:
    JacobianFormSurf1Euler(int i, int j) : MatrixFormSurf<double>(i, j)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
      }

      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* clone() const
    {
      return new JacobianFormSurf1Euler(this->i, this->j);
    }
  };

  class JacobianFormSurf1CrankNicolson : public MatrixFormSurf<double>
  {
  public:
    JacobianFormSurf1CrankNicolson(int i, int j) : MatrixFormSurf<double>(i, j)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) {
        // Just the implicit Euler contributes:
        result += wt[i] * 0.5 * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
      }

      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* clone() const
    {
      return new JacobianFormSurf1CrankNicolson(this->i, this->j);
    }
  };

  class JacobianFormSurf4Euler : public MatrixFormSurf<double>
  {
  public:
    JacobianFormSurf4Euler(int i, int j) : MatrixFormSurf<double>(i, j)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) {
        result -= wt[i] * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
      }

      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* clone() const
    {
      return new JacobianFormSurf4Euler(this->i, this->j);
    }
  };

  class JacobianFormSurf4CrankNicolson : public MatrixFormSurf<double>
  {
  public:
    JacobianFormSurf4CrankNicolson(int i, int j) : MatrixFormSurf<double>(i, j)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) {
        // Just the implicit Euler contributes:
        result -= wt[i] * 0.5 * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
      }

      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* clone() const
    {
      return new JacobianFormSurf4CrankNicolson(this->i, this->j);
    }
  };

  class JacobianFormSurf6Euler : public MatrixFormSurf<double>
  {
  public:
    JacobianFormSurf6Euler(int i, int j) : MatrixFormSurf<double>(i, j)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
      }

      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* clone() const
    {
      return new JacobianFormSurf6Euler(this->i, this->j);
    }
  };

  class JacobianFormSurf6CrankNicolson : public MatrixFormSurf<double>
  {
  public:
    JacobianFormSurf6CrankNicolson(int i, int j) : MatrixFormSurf<double>(i, j)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      // Just the implicit Euler contributes:
      for (int i = 0; i < n; i++) {
        result += wt[i] * 0.5 * dKdh(h_prev_newton->val[i]) * u->val[i] * v->val[i];
      }

      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    MatrixFormSurf<double>* clone() const
    {
      return new JacobianFormSurf6CrankNicolson(this->i, this->j);
    }
  };

  class ResidualFormVolEuler : public VectorFormVol<double>
  {
  public:
    ResidualFormVolEuler(int i) : VectorFormVol<double>(i)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * (
          C(h_prev_newton->val[i]) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / TAU
          + K(h_prev_newton->val[i]) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
          - dKdh(h_prev_newton->val[i]) * h_prev_newton->dy[i] * v->val[i]
        );
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomVol<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    VectorFormVol<double>* clone() const
    {
      return new ResidualFormVolEuler(this->i);
    }
  };

  class ResidualFormVolCrankNicolson : public VectorFormVol<double>
  {
  public:
    ResidualFormVolCrankNicolson(int i) : VectorFormVol<double>(i)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * 0.5 * ( // implicit Euler part
          C(h_prev_newton->val[i]) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / TAU
          + K(h_prev_newton->val[i]) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
          - dKdh(h_prev_newton->val[i]) * h_prev_newton->dy[i] * v->val[i]
        )
          + wt[i] * 0.5 * ( // explicit Euler part
          C(h_prev_time->val[i]) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / TAU
          + K(h_prev_time->val[i]) * (h_prev_time->dx[i] * v->dx[i] + h_prev_time->dy[i] * v->dy[i])
          - dKdh(h_prev_time->val[i]) * h_prev_time->dy[i] * v->val[i]
        );
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomVol<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    VectorFormVol<double>* clone() const
    {
      return new ResidualFormVolEuler(this->i);
    }
  };

  class ResidualFormSurf1Euler : public VectorFormSurf<double>
  {
  public:
    ResidualFormSurf1Euler(int i) : VectorFormSurf<double>(i)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * K(h_prev_newton->val[i]) * v->val[i];
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    VectorFormSurf<double>* clone() const
    {
      return new ResidualFormSurf1Euler(this->i);
    }
  };

  class ResidualFormSurf1CrankNicolson : public VectorFormSurf<double>
  {
  public:
    ResidualFormSurf1CrankNicolson(int i) : VectorFormSurf<double>(i)
    {
    }

    virtual double value(GeomSurf n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * 0.5 * (K(h_prev_newton->val[i]) + K(h_prev_time->val[i])) * v->val[i];
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    VectorFormSurf<double>* clone() const
    {
      return new ResidualFormSurf1CrankNicolson(this->i);
    }
  };

  class ResidualFormSurf4Euler : public VectorFormSurf<double>
  {
  public:
    ResidualFormSurf4Euler(int i) : VectorFormSurf<double>(i)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) {
        result -= wt[i] * K(h_prev_newton->val[i]) * v->val[i];
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    VectorFormSurf<double>* clone() const
    {
      return new ResidualFormSurf4Euler(this->i);
    }
  };

  class ResidualFormSurf4CrankNicolson : public VectorFormSurf<double>
  {
  public:
    ResidualFormSurf4CrankNicolson(int i) : VectorFormSurf<double>(i)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];
      for (int i = 0; i < n; i++) {
        result -= wt[i] * 0.5 * (K(h_prev_newton->val[i]) + K(h_prev_time->val[i]))* v->val[i];
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    VectorFormSurf<double>* clone() const
    {
      return new ResidualFormSurf4CrankNicolson(this->i, this->j);
    }
  };

  class ResidualFormSurf6Euler : public VectorFormSurf<double>
  {
  public:
    ResidualFormSurf6Euler(int i) : VectorFormSurf<double>(i)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * (q_function() + K(h_prev_newton->val[i])) * v->val[i];
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    VectorFormSurf<double>* clone() const
    {
      return new ResidualFormSurf6Euler(this->i);
    }
  };

  class ResidualFormSurf6CrankNicolson : public VectorFormSurf<double>
  {
  public:
    ResidualFormSurf6CrankNicolson(int i) : VectorFormSurf<double>(i)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, GeomSurf<double> *e, Func<double>* *ext) const
    {
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * (q_function() + 0.5 * (K(h_prev_newton->val[i]) + K(h_prev_time->val[i]))) * v->val[i];
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      GeomSurf<Ord> *e, Func<Ord>* *ext) const
    {
      return Ord(20);
    }

    VectorFormSurf<double>* clone() const
    {
      return new ResidualFormSurf6CrankNicolson(this->i);
    }
  };
};