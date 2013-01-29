#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

class CustomInitialConditionWave : public ExactSolutionScalar<double>
{
public:
  CustomInitialConditionWave(const Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual double value (double x, double y) const {
    return exp(-x*x - y*y);
  }

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = exp(-x*x - y*y) * (-2*x);
    dy = exp(-x*x - y*y) * (-2*y);
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(10);
  }

  MeshFunction<double>* clone() const
  {
    return new CustomInitialConditionWave(mesh);
  }
};

/* Weak forms */

class CustomWeakFormWave : public WeakForm<double>
{
public:

  CustomWeakFormWave(double tau, double c_squared, Solution<double>* u_prev_sln, Solution<double>* v_prev_sln) : WeakForm(2) {
    add_matrix_form(new MatrixFormVolWave_0_1());
    add_matrix_form(new MatrixFormVolWave_1_0(c_squared));

    VectorFormVolWave_0* vector_form_0 = new VectorFormVolWave_0();
    vector_form_0->set_ext(v_prev_sln);
    add_vector_form(vector_form_0);
    
    VectorFormVolWave_1* vector_form_1 = new VectorFormVolWave_1(c_squared);
    vector_form_1->set_ext(u_prev_sln);
    add_vector_form(vector_form_1);
  };

private:
  // This form is custom because of the clone() method.
  class MatrixFormVolWave_0_1 : public MatrixFormVol<double>
  {
  public:
    MatrixFormVolWave_0_1() : MatrixFormVol<double>(0, 1) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    virtual MatrixFormVol<double>* clone() const 
    {
      return new MatrixFormVolWave_0_1(*this);
    }
  };

  // This form is custom because of the clone() method.
  class MatrixFormVolWave_1_0 : public MatrixFormVol<double>
  {
  public:
    MatrixFormVolWave_1_0(double c_squared) 
          : MatrixFormVol<double>(1, 0), c_squared(c_squared) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar> **ext) const {
      return - c_squared * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double> **ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, Func<Ord> **ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    virtual MatrixFormVol<double>* clone() const {
      return new MatrixFormVolWave_1_0(*this);
    }

    double c_squared;
  };

  class VectorFormVolWave_0 : public VectorFormVol<double>
  {
  public:
    VectorFormVolWave_0() : VectorFormVol<double>(0) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, Func<Scalar> **ext) const {
      Scalar result = Scalar(0);
      Func<Scalar>* K_sln = u_ext[0];
      Func<Scalar>* sln_prev_time = ext[0];

      for (int i = 0; i < n; i++) {
        Scalar sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
        result += wt[i] * sln_val_i * v->val[i];
      }
      return result;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                 Geom<double> *e, Func<double> **ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            Func<Ord> **ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    virtual VectorFormVol<double>* clone() const {
      return new VectorFormVolWave_0(*this);
    }
  };

  class VectorFormVolWave_1 : public VectorFormVol<double>
  {
  public:
    VectorFormVolWave_1(double c_squared) 
          : VectorFormVol<double>(1), c_squared(c_squared) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, Func<Scalar> **ext) const {
      Scalar result = Scalar(0);
      Func<Scalar>* K_sln = u_ext[0];
      Func<Scalar>* sln_prev_time = ext[0];
      
      for (int i = 0; i < n; i++) {
        Scalar sln_dx_i = sln_prev_time->dx[i] + K_sln->dx[i];
        Scalar sln_dy_i = sln_prev_time->dy[i] + K_sln->dy[i];
        result += wt[i] * c_squared * (sln_dx_i * v->dx[i] 
                  + sln_dy_i * v->dy[i]);
      }
      return - result;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                 Geom<double> *e, Func<double> **ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            Func<Ord> **ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    virtual VectorFormVol<double>* clone() const {
      return new VectorFormVolWave_1(*this);
    }

    double c_squared;
  };
};
