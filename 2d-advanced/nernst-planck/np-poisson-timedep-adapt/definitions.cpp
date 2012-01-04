#include "definitions.h"

class ScaledWeakFormPNPCranic : public WeakForm<double> {
public:
  ScaledWeakFormPNPCranic(double* tau, double epsilon,
        Solution<double>* C_prev_time, Solution<double>* phi_prev_time) : WeakForm<double>(2) {
      for(unsigned int i = 0; i < 2; i++) {
        ScaledWeakFormPNPCranic::Residual* vector_form =
            new ScaledWeakFormPNPCranic::Residual(i, tau, epsilon);
        if(i == 0) {
          vector_form->ext.push_back(C_prev_time);
          vector_form->ext.push_back(phi_prev_time);
        }
        add_vector_form(vector_form);
        for(unsigned int j = 0; j < 2; j++)
          add_matrix_form(new ScaledWeakFormPNPCranic::Jacobian(i, j, tau, epsilon));
      }
    };

private:
  class Jacobian : public MatrixFormVol<double> {
  public:
    Jacobian(int i, int j, double* tau, double epsilon) : MatrixFormVol<double>(i, j),
          i(i), j(j), tau(tau), epsilon(epsilon) {}

    template<typename Real, typename Scalar>
    Real matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Real result = Real(0);
      Func<Scalar>* prev_newton;
      switch(i * 10 + j) {
        case 0:
          prev_newton = u_ext[1];
          for (int i = 0; i < n; i++) {

            result += wt[i] * (u->val[i] * v->val[i] / *(this->tau) +
                this->epsilon * 0.5 * ((u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) +
                    u->val[i] * (prev_newton->dx[i] * v->dx[i] + prev_newton->dy[i] * v->dy[i])));
          }
          return result;
          break;
        case 1:
          prev_newton = u_ext[0];
          for (int i = 0; i < n; i++) {
            result += wt[i] * (0.5 * this->epsilon * prev_newton->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
          }
          return result;
          break;
        case 10:
          for (int i = 0; i < n; i++) {
            result += wt[i] * ( -1.0/(2 * this->epsilon * this->epsilon) * u->val[i] * v->val[i]);
          }
          return result;
          break;
        case 11:
          for (int i = 0; i < n; i++) {
            result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          }
          return result;
          break;
        default:

          return result;
      }
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormVol<double>* clone()
    {
      return new Jacobian(*this);
    }

    // Members.
    int i, j;
    double* tau;
    double epsilon;
  };

  class Residual : public VectorFormVol<double>
      {
      public:
        Residual(int i, double* tau, double epsilon)
          : VectorFormVol<double>(i), i(i), tau(tau), epsilon(epsilon) {}

        template<typename Real, typename Scalar>
        Real vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
          Real result = Real(0);
          Func<Scalar>* C_prev_time;
          Func<Scalar>* phi_prev_time;
          Func<Scalar>* C_prev_newton;
          Func<Scalar>* phi_prev_newton;
          switch(i) {
            case 0:
              C_prev_time = ext->fn[0];
              phi_prev_time = ext->fn[1];
              C_prev_newton = u_ext[0];
              phi_prev_newton = u_ext[1];
              for (int i = 0; i < n; i++) {
                result += wt[i] * ((C_prev_newton->val[i] - C_prev_time->val[i]) * v->val[i] / *(this->tau) +
                    0.5 * this->epsilon * ((C_prev_newton->dx[i] * v->dx[i] + C_prev_newton->dy[i] * v->dy[i]) +
                    (C_prev_time->dx[i] * v->dx[i] + C_prev_time->dy[i] * v->dy[i]) +
                      C_prev_newton->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
                      C_prev_time->val[i] * (phi_prev_time->dx[i] * v->dx[i] + phi_prev_time->dy[i] * v->dy[i])));
              }
              return result;
              break;
            case 1:
              C_prev_newton = u_ext[0];
              phi_prev_newton = u_ext[1];
              for (int i = 0; i < n; i++) {
                result += wt[i] * ((phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
                    v->val[i] * 1 / (2 * this->epsilon * this->epsilon) * (1 - C_prev_newton->val[i]));
              }
              return result;
              break;
            default:
              return result;
          }
        }

        virtual double value(int n, double *wt, Func<double> *u_ext[],
                     Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
          return vector_form<double, double>(n, wt, u_ext, v, e, ext);
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                Geom<Ord> *e, ExtData<Ord> *ext) const {
          return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
        }
        
        VectorFormVol<double>* clone()
        {
          return new Residual(*this);
        }

        // Members.
        int i;
        double* tau;
        double epsilon;
      };
};


class WeakFormPNPCranic : public WeakForm<double> {
public:

  WeakFormPNPCranic(double* tau, double C0, double K, double L, double D,
      Solution<double>* C_prev_time, Solution<double>* phi_prev_time) : WeakForm<double>(2) {
    for(unsigned int i = 0; i < 2; i++) {
      WeakFormPNPCranic::Residual* vector_form =
          new WeakFormPNPCranic::Residual(i, tau, C0, K, L, D);
      if(i == 0) {
        vector_form->ext.push_back(C_prev_time);
        vector_form->ext.push_back(phi_prev_time);
      }
      add_vector_form(vector_form);
      for(unsigned int j = 0; j < 2; j++)
        add_matrix_form(new WeakFormPNPCranic::Jacobian(
            i, j, tau, C0, K, L, D));
    }
  };

private:
  class Jacobian : public MatrixFormVol<double> {
  public:
    Jacobian(int i, int j, double* tau, double C0, double K, double L, double D) : MatrixFormVol<double>(i, j),
          i(i), j(j), tau(tau), C0(C0), L(L), D(D) {}

    template<typename Real, typename Scalar>
    Real matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Real result = Real(0);
      Func<Scalar>* prev_newton;
      switch(i * 10 + j) {
        case 0:
          prev_newton = u_ext[1];
          for (int i = 0; i < n; i++) {

            result += wt[i] * (u->val[i] * v->val[i] / *(this->tau) +
                this->D * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) +
                this->K * u->val[i] * (prev_newton->dx[i] * v->dx[i] + prev_newton->dy[i] * v->dy[i]));
          }
          return result;
          break;
        case 1:
          prev_newton = u_ext[0];
          for (int i = 0; i < n; i++) {
            result += wt[i] * (0.5 * this->K * prev_newton->val[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
          }
          return result;
          break;
        case 10:
          for (int i = 0; i < n; i++) {
            result += wt[i] * ( -this->L * u->val[i] * v->val[i]);
          }
          return result;
          break;
        case 11:
          for (int i = 0; i < n; i++) {
            result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          }
          return result;
          break;
        default:

          return result;
      }
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
    
    MatrixFormVol<double>* clone()
    {
      return new Jacobian(*this);
    }

    // Members.
    int i, j;
    double* tau;
    double C0;
    double K;
    double L;
    double D;
  };

  class Residual : public VectorFormVol<double>
    {
    public:
      Residual(int i, double* tau, double C0, double K, double L, double D)
        : VectorFormVol<double>(i), i(i), tau(tau), C0(C0), K(K), L(L), D(D) {}

      template<typename Real, typename Scalar>
      Real vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                          Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Real result = Real(0);
        Func<Scalar>* C_prev_time;
        Func<Scalar>* phi_prev_time;
        Func<Scalar>* C_prev_newton;
        Func<Scalar>* phi_prev_newton;
        Func<Scalar>* u1_prev_newton;
        Func<Scalar>* u2_prev_newton;
        switch(i) {
          case 0:
            C_prev_time = ext->fn[0];
            phi_prev_time = ext->fn[1];
            C_prev_newton = u_ext[0];
            phi_prev_newton = u_ext[1];
            for (int i = 0; i < n; i++) {
              result += wt[i] * ((C_prev_newton->val[i] - C_prev_time->val[i]) * v->val[i] / *(this->tau) +
                  0.5 * this->D * (C_prev_newton->dx[i] * v->dx[i] + C_prev_newton->dy[i] * v->dy[i]) +
                  0.5 * this->D * (C_prev_time->dx[i] * v->dx[i] + C_prev_time->dy[i] * v->dy[i]) +
                  0.5 * this->K * C_prev_newton->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
                  0.5 * this->K * C_prev_time->val[i] * (phi_prev_time->dx[i] * v->dx[i] + phi_prev_time->dy[i] * v->dy[i]));
            }
            return result;
            break;
          case 1:
            C_prev_newton = u_ext[0];
            phi_prev_newton = u_ext[1];
            for (int i = 0; i < n; i++) {
              result += wt[i] * ((phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
                    this->L * v->val[i] * (this->C0 - C_prev_newton->val[i]));
            }
            return result;
            break;
          default:
            return result;
        }
      }

      virtual double value(int n, double *wt, Func<double> *u_ext[],
                   Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
        return vector_form<double, double>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
              Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }
      
      VectorFormVol<double>* clone()
      {
        return new Residual(*this);
      }

      // Members.
      int i;
      double* tau;
      double C0;
      double K;
      double L;
      double D;
    };

};

class WeakFormPNPEuler : public WeakForm<double>
{
public:
  WeakFormPNPEuler(double* tau, double C0, double K, double L, double D, Solution<double>* C_prev_time)
    : WeakForm<double>(2) {
    for(unsigned int i = 0; i < 2; i++) {
      WeakFormPNPEuler::Residual* vector_form =
          new WeakFormPNPEuler::Residual(i, tau, C0, K, L, D);
      if(i == 0)
        vector_form->ext.push_back(C_prev_time);
      add_vector_form(vector_form);
      for(unsigned int j = 0; j < 2; j++)
        add_matrix_form(new WeakFormPNPEuler::Jacobian(i, j, tau, C0, K, L, D));
    }
  };

private:
  class Jacobian : public MatrixFormVol<double>
  {
  public:
    Jacobian(int i, int j, double* tau, double C0, double K, double L, double D)
      : MatrixFormVol<double>(i, j), i(i), j(j), tau(tau), C0(C0), K(K), L(L), D(D) {}

    template<typename Real, typename Scalar>
    Real matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Real result = Real(0);
	    Func<Scalar>* prev_newton;
      switch(i * 10 + j) {
        case 0:
	        prev_newton = u_ext[1];
	        for (int i = 0; i < n; i++) {
	          result += wt[i] * (u->val[i] * v->val[i] / *(this->tau) +
	              0.5 * this->D * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) +
	              0.5 * this->K * u->val[i] * (prev_newton->dx[i] * v->dx[i] + prev_newton->dy[i] * v->dy[i]));
	        }
	        return result;
          break;
        case 1:
	        prev_newton = u_ext[0];
	        for (int i = 0; i < n; i++) {
		        result += wt[i] * this->K * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) * prev_newton->val[i];
	        }
	        return result;
          break;
        case 10:
	        for (int i = 0; i < n; i++) {
		        result += wt[i] * ( -this->L * u->val[i] * v->val[i]);
	        }
	        return result;
          break;
        case 11:
	        for (int i = 0; i < n; i++) {
		        result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
	        }
	        return result;
          break;
        default:
          
          return result;
      }
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
    
    MatrixFormVol<double>* clone()
    {
      return new Jacobian(*this);
    }

    // Members.
    int i, j;
    double* tau;
    double C0;
    double K;
    double L;
    double D;
  };

  class Residual : public VectorFormVol<double>
  {
  public:
    Residual(int i, double* tau, double C0, double K, double L, double D)
      : VectorFormVol<double>(i), i(i), tau(tau), C0(C0), K(K), L(L), D(D) {}

    template<typename Real, typename Scalar>
    Real vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                        Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Real result = Real(0);
      Func<Scalar>* C_prev_time;
	    Func<Scalar>* C_prev_newton;
	    Func<Scalar>* phi_prev_newton;
      Func<Scalar>* u1_prev_newton;
      Func<Scalar>* u2_prev_newton;
      switch(i) {
        case 0:
	        C_prev_time = ext->fn[0];
	        C_prev_newton = u_ext[0];
	        phi_prev_newton = u_ext[1];
	        for (int i = 0; i < n; i++) {
		        result += wt[i] * ((C_prev_newton->val[i] - C_prev_time->val[i]) * v->val[i] / *(this->tau) +
				        this->D * (C_prev_newton->dx[i] * v->dx[i] + C_prev_newton->dy[i] * v->dy[i]) +
				        this->K * C_prev_newton->val[i] * (phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]));
	        }
	        return result;
          break;
        case 1:
	        C_prev_newton = u_ext[0];
	        phi_prev_newton = u_ext[1];
	        for (int i = 0; i < n; i++) {
	          result += wt[i] * ((phi_prev_newton->dx[i] * v->dx[i] + phi_prev_newton->dy[i] * v->dy[i]) +
					        this->L * v->val[i] * (this->C0 - C_prev_newton->val[i]));
	        }
	        return result;
          break;
        default:
          return result;
      }
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[],
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
    
    VectorFormVol<double>* clone()
    {
      return new Residual(*this);
    }

    // Members.
    int i;
    double* tau;
    double C0;
    double K;
    double L;
    double D;
  };
};
