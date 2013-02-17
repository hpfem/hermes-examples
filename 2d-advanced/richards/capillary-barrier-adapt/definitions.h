#include "hermes2d.h"

#include "../constitutive.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// The first part of the file dontains forms for the Newton's
// method. Forms for Picard are in the second part.

/*** INITIAL CONDITION ***/

class InitialSolutionRichards : public ExactSolutionScalar<double>
{
public:
  InitialSolutionRichards(const Mesh* mesh, double constant) 
         : ExactSolutionScalar<double>(mesh), constant(constant) {};

  virtual double value (double x, double y) const {
    return constant;
  }

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0.0;
    dy = 0.0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(0);
  }

  virtual MeshFunction<double>* clone() const
  {
    return new InitialSolutionRichards(mesh, constant);
  }

  // Value.
  double constant;
};

class ExactSolutionPoisson : public ExactSolutionScalar<double>
{
public:
  ExactSolutionPoisson(const Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual double value (double x, double y) const {
    return x*x +y*y;
  }

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 2*x;
    dy = 2*y;
  };
  
  virtual MeshFunction<double>* clone() const
  {
    return new ExactSolutionPoisson(mesh);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return x*x +y*y;
  }
};

/*** NEWTON ***/

class WeakFormRichardsNewtonEuler : public WeakForm<double>
{
public:
  WeakFormRichardsNewtonEuler(ConstitutiveRelationsGenuchtenWithLayer* relations, double tau, Solution<double>* prev_time_sln, Mesh* mesh) 
               : WeakForm<double>(1), mesh(mesh), relations(relations) {
    JacobianFormNewtonEuler* jac_form = new JacobianFormNewtonEuler(0, 0, relations, tau);
    jac_form->set_ext(prev_time_sln);
    add_matrix_form(jac_form);

    ResidualFormNewtonEuler* res_form = new ResidualFormNewtonEuler(0, relations, tau);
    res_form->set_ext(prev_time_sln);
    add_vector_form(res_form);
  }

private:
  class JacobianFormNewtonEuler : public MatrixFormVol<double>
  {
  public:
    JacobianFormNewtonEuler(int i, int j, ConstitutiveRelationsGenuchtenWithLayer* relations, double tau) 
      : MatrixFormVol<double>(i, j), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const {
                         std::string elem_marker = static_cast<WeakFormRichardsNewtonEuler*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];

      for (int i = 0; i < n; i++)
        result += wt[i] * (
        relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
		             + relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * h_prev_newton->val[i] * v->val[i] / tau
		             - relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * h_prev_time->val[i] * v->val[i] / tau
			     + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                             + relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * 
                               (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                             - relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->dy[i] * v->val[i]
                             - relations->ddKdhh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                          );
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, Func<double>* *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
            Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      return Ord(30);
    }

    MatrixFormVol<double>* clone() const
    {
      JacobianFormNewtonEuler* form = new JacobianFormNewtonEuler(i, j, relations, tau);
      form->wf = this->wf;
      return form;
    }

    // Members.
    double tau;
    ConstitutiveRelationsGenuchtenWithLayer* relations;
  };

  class ResidualFormNewtonEuler : public VectorFormVol<double>
  {
  public:
    ResidualFormNewtonEuler(int i, ConstitutiveRelationsGenuchtenWithLayer* relations, double tau) 
               : VectorFormVol<double>(i), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const {
      std::string elem_marker = static_cast<WeakFormRichardsNewtonEuler*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * (
		           relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                           * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / tau
                           + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                           * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                           -relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * h_prev_newton->dy[i] * v->val[i]
                          );
      }
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      return Ord(30);
    }
    
    VectorFormVol<double>* clone() const
    {
      ResidualFormNewtonEuler* form = new ResidualFormNewtonEuler(i, relations, tau);
      form->wf = this->wf;
      return form;
    }

    // Members.
    double tau;
    ConstitutiveRelationsGenuchtenWithLayer* relations;
  };

  Mesh* mesh;

  ConstitutiveRelationsGenuchtenWithLayer* relations;

  WeakForm<double>* clone() const
  {
    WeakFormRichardsNewtonEuler* wf = new WeakFormRichardsNewtonEuler(*this);
    wf->relations = this->relations;
    return wf;
  }
};


class WeakFormRichardsNewtonCrankNicolson : public WeakForm<double>
{
public:
  WeakFormRichardsNewtonCrankNicolson(ConstitutiveRelationsGenuchtenWithLayer* relations, double tau, Solution<double>* prev_time_sln, Mesh* mesh) : WeakForm<double>(1), relations(relations), mesh(mesh) {
    JacobianFormNewtonCrankNicolson* jac_form = new JacobianFormNewtonCrankNicolson(0, 0, relations, tau);
    jac_form->set_ext(prev_time_sln);
    add_matrix_form(jac_form);

    ResidualFormNewtonCrankNicolson* res_form = new ResidualFormNewtonCrankNicolson(0, relations, tau);
    res_form->set_ext(prev_time_sln);
    add_vector_form(res_form);
  }

private:
  class JacobianFormNewtonCrankNicolson : public MatrixFormVol<double>
  {
  public:
    JacobianFormNewtonCrankNicolson(int i, int j, ConstitutiveRelationsGenuchtenWithLayer* relations, double tau) 
      : MatrixFormVol<double>(i, j), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, Func<Scalar>* *ext) const {
      std::string elem_marker = static_cast<WeakFormRichardsNewtonCrankNicolson*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * 0.5 * ( // implicit Euler part:
		             relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
		             + relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * h_prev_newton->val[i] * v->val[i] / tau
		             - relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * h_prev_time->val[i] * v->val[i] / tau
			     + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                             + relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * 
                               (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                             - relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->dy[i] * v->val[i]
                             - relations->ddKdhh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                           )
                + wt[i] * 0.5 * ( // explicit Euler part, 
		             relations->C(h_prev_time->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
                           );
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }
    
    MatrixFormVol<double>* clone() const
    {
      JacobianFormNewtonCrankNicolson* form = new JacobianFormNewtonCrankNicolson(i, j, relations, tau);
      form->wf = this->wf;
      return form;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      return Ord(30);
    }

    // Members.
    double tau;
    ConstitutiveRelationsGenuchtenWithLayer* relations;
  };

  class ResidualFormNewtonCrankNicolson : public VectorFormVol<double>
  {
  public:
    ResidualFormNewtonCrankNicolson(int i, ConstitutiveRelationsGenuchtenWithLayer* relations, double tau) : VectorFormVol<double>(i), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const {
      double result = 0;
      std::string elem_marker = static_cast<WeakFormRichardsNewtonCrankNicolson*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * 0.5 * ( // implicit Euler part
		           relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / tau
                           + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                           - relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * h_prev_newton->dy[i] * v->val[i]
                          )
                + wt[i] * 0.5 * ( // explicit Euler part
		           relations->C(h_prev_time->val[i], atoi(elem_marker.c_str())) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / tau
                           + relations->K(h_prev_time->val[i], atoi(elem_marker.c_str())) * (h_prev_time->dx[i] * v->dx[i] + h_prev_time->dy[i] * v->dy[i])
                           - relations->dKdh(h_prev_time->val[i], atoi(elem_marker.c_str())) * h_prev_time->dy[i] * v->val[i]
		           );
      }
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      return Ord(30);
    }
    
    VectorFormVol<double>* clone() const
    {
      ResidualFormNewtonCrankNicolson* form = new ResidualFormNewtonCrankNicolson(i, relations, tau);
      form->wf = this->wf;
      return form;
    }

    // Members.
    double tau;
    ConstitutiveRelationsGenuchtenWithLayer* relations;
  };

  Mesh* mesh;

  ConstitutiveRelationsGenuchtenWithLayer* relations;

  WeakForm<double>* clone() const
  {
    WeakFormRichardsNewtonCrankNicolson* wf = new WeakFormRichardsNewtonCrankNicolson(*this);
    wf->relations = this->relations;
    return wf;
  }
};


class WeakFormRichardsPicardEuler : public WeakForm<double>
{
public:
  WeakFormRichardsPicardEuler(ConstitutiveRelationsGenuchtenWithLayer* relations, double tau, Solution<double>* prev_picard_sln, Solution<double>* prev_time_sln, Mesh* mesh) : WeakForm<double>(1), relations(relations), mesh(mesh) {
    JacobianFormPicardEuler* jac_form = new JacobianFormPicardEuler(0, 0, relations, tau);
    jac_form->set_ext(prev_picard_sln);
    add_matrix_form(jac_form);

    ResidualFormPicardEuler* res_form = new ResidualFormPicardEuler(0, relations, tau);
    res_form->set_ext(prev_picard_sln);
    res_form->set_ext(prev_time_sln);
    add_vector_form(res_form);
  }

private:
  class JacobianFormPicardEuler : public MatrixFormVol<double>
  {
  public:
    JacobianFormPicardEuler(int i, int j, ConstitutiveRelationsGenuchtenWithLayer* relations, double tau) 
      : MatrixFormVol<double>(i, j), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, Func<Scalar>* *ext) const {
      std::string elem_marker = static_cast<WeakFormRichardsPicardEuler*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      double result = 0;
      Func<double>* h_prev_picard = ext[0];

      for (int i = 0; i < n; i++) {
        result += wt[i] * (  relations->C(h_prev_picard->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
                             + relations->K(h_prev_picard->val[i], atoi(elem_marker.c_str())) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                             - relations->dKdh(h_prev_picard->val[i], atoi(elem_marker.c_str())) * u->dy[i] * v->val[i]);
      }
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, Func<double>* *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }
    
    MatrixFormVol<double>* clone() const
    {
      JacobianFormPicardEuler* form = new JacobianFormPicardEuler(i, j, relations, tau);
      form->wf = this->wf;
      return form;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      return Ord(30);
    }

    // Members.
    double tau;
    ConstitutiveRelationsGenuchtenWithLayer* relations;
  };

  class ResidualFormPicardEuler : public VectorFormVol<double>
  {
  public:
    ResidualFormPicardEuler(int i, ConstitutiveRelationsGenuchtenWithLayer* relations, double tau) : VectorFormVol<double>(i), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, Func<Scalar>* *ext) const {
      std::string elem_marker = static_cast<WeakFormRichardsPicardEuler*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      double result = 0;
      Func<double>* h_prev_picard = ext[0];
      Func<double>* h_prev_time = ext[1];
      for (int i = 0; i < n; i++) 
        result += wt[i] * relations->C(h_prev_picard->val[i], atoi(elem_marker.c_str())) * h_prev_time->val[i] * v->val[i] / tau;
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, Func<double>* *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }
    
    VectorFormVol<double>* clone() const
    {
      ResidualFormPicardEuler* form = new ResidualFormPicardEuler(i, relations, tau);
      form->wf = this->wf;
      return form;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const {
      return Ord(30);
    }
    
    // Members.
    double tau;
    ConstitutiveRelationsGenuchtenWithLayer* relations;
  };

  Mesh* mesh;

  ConstitutiveRelationsGenuchtenWithLayer* relations;

  WeakForm<double>* clone() const
  {
    WeakFormRichardsPicardEuler* wf = new WeakFormRichardsPicardEuler(*this);
    wf->relations = this->relations;
    return wf;
  }
};

class RichardsEssentialBC : public EssentialBoundaryCondition<double> {
public:

  RichardsEssentialBC(std::string marker, double h_elevation, double pulse_end_time, double h_init, double startup_time) :
  EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), h_elevation(h_elevation), pulse_end_time(pulse_end_time), h_init(h_init), startup_time(startup_time)
  {
    markers.push_back(marker);
  }

  ~RichardsEssentialBC() {}

  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<double>::BC_FUNCTION; }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
    if (current_time < startup_time)
      return h_init + current_time/startup_time*(h_elevation-h_init);
    else if (current_time > pulse_end_time)
      return h_init;
    else
      return h_elevation;
  }

  // Member.
  double h_elevation;
  double pulse_end_time;
  double h_init;
  double startup_time;
};