#include "hermes2d.h"

// Numerical fluxes.
#include "numerical_flux.h"

// Utility functions for the Euler equations.
#include "euler_util.h"

class EulerEquationsWeakFormStabilization : public WeakForm<double>
{
public:
  EulerEquationsWeakFormStabilization(Solution<double>* prev_rho) : WeakForm<double>()
  {
    DGVectorFormIndicator* form = new DGVectorFormIndicator();
    form->setExt(prev_rho);
    add_vector_form_DG(form);
  }

  class DGVectorFormIndicator : public VectorFormDG<double>
  {
  public:
    DGVectorFormIndicator() 
      : VectorFormDG<double>(0)
    {
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0;
      double w_L[4], w_R[4];
      for (int i = 0;i < n;i++)
        result += wt[i] * v->val[i] * (ext->fn[0]->get_val_central(i) - ext->fn[0]->get_val_neighbor(i)) * (ext->fn[0]->get_val_central(i) - ext->fn[0]->get_val_neighbor(i)) / (e->diam * std::pow(e->area, 0.75));

      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return v->val[0] * v->val[0] * Ord(6);
    }

    VectorFormDG<double>* clone() { return new DGVectorFormIndicator; }
  };
};

class EulerEquationsWeakFormExplicit : public WeakForm<double>
{
public:
  // Constructor.
  EulerEquationsWeakFormExplicit(double kappa, double rho_ext, double v1_ext, double v2_ext, double pressure_ext, 
    std::string solid_wall_bottom_marker, std::string solid_wall_top_marker, std::string inlet_marker, std::string outlet_marker, 
    Solution<double>* prev_density, Solution<double>* prev_density_vel_x, Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, bool fvm_only = false, int num_of_equations = 4) :
  WeakForm<double>(num_of_equations), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), pressure_ext(pressure_ext), 
    energy_ext(QuantityCalculator::calc_energy(rho_ext, rho_ext * v1_ext, rho_ext * v2_ext, pressure_ext, kappa)), euler_fluxes(new EulerFluxes(kappa)) 
  {
    add_matrix_form(new EulerEquationsBilinearFormTime(0));
    add_matrix_form(new EulerEquationsBilinearFormTime(1));
    add_matrix_form(new EulerEquationsBilinearFormTime(2));
    add_matrix_form(new EulerEquationsBilinearFormTime(3));
    add_vector_form(new EulerEquationsLinearFormTime(0));
    add_vector_form(new EulerEquationsLinearFormTime(1));
    add_vector_form(new EulerEquationsLinearFormTime(2));
    add_vector_form(new EulerEquationsLinearFormTime(3));

    if(!fvm_only) 
    {
      add_vector_form(new EulerEquationsLinearFormDensity());
      add_vector_form(new EulerEquationsLinearFormDensityVelX(kappa));
      add_vector_form(new EulerEquationsLinearFormDensityVelY(kappa));
      add_vector_form(new EulerEquationsLinearFormEnergy(kappa));
    }

    add_vector_form_DG(new EulerEquationsLinearFormInterface(0, kappa));
    add_vector_form_DG(new EulerEquationsLinearFormInterface(1, kappa));
    add_vector_form_DG(new EulerEquationsLinearFormInterface(2, kappa));
    add_vector_form_DG(new EulerEquationsLinearFormInterface(3, kappa));

    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(0, solid_wall_bottom_marker, kappa));
    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(1, solid_wall_bottom_marker, kappa));
    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(2, solid_wall_bottom_marker, kappa));
    add_vector_form_surf(new EulerEquationsLinearFormSolidWall(3, solid_wall_bottom_marker, kappa));

    if(solid_wall_bottom_marker != solid_wall_top_marker) 
    {
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(0, solid_wall_top_marker, kappa));
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(1, solid_wall_top_marker, kappa));
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(2, solid_wall_top_marker, kappa));
      add_vector_form_surf(new EulerEquationsLinearFormSolidWall(3, solid_wall_top_marker, kappa));
    }
    else
      Hermes::Mixins::Loggable::Static::warn("Are you sure that solid wall top and bottom markers should coincide?");

    add_vector_form_surf(new EulerEquationsLinearFormInlet(0, inlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsLinearFormInlet(1, inlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsLinearFormInlet(2, inlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsLinearFormInlet(3, inlet_marker, kappa));

    add_vector_form_surf(new EulerEquationsLinearFormOutlet(0, outlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsLinearFormOutlet(1, outlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsLinearFormOutlet(2, outlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsLinearFormOutlet(3, outlet_marker, kappa));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfvol.size();vector_form_i++) 
      vfvol.at(vector_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfsurf.size();vector_form_i++) 
      vfsurf.at(vector_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfDG.size();vector_form_i++) 
      vfDG.at(vector_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
  };

  void set_time_step(double tau) 
  {
    this->tau = tau;
  }

  double get_tau() const 
  {
    return tau;
  }

  // Destructor.
  ~EulerEquationsWeakFormExplicit() {};
protected:
  class EulerEquationsBilinearFormTime : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormTime(int i) : MatrixFormVol<double>(i, i) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormVol<double>* clone() { return new EulerEquationsBilinearFormTime(this->i); }
  };

  class EulerEquationsLinearFormDensity : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormDensity() : VectorFormVol<double>(0) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, 
      ExtData<Scalar> *ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * ext->fn[0]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_0_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_0_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_0_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_0_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_0_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_0_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_0_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_0_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, 
      ExtData<double> *ext) const 
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }
    VectorFormVol<double>* clone()
    {
      EulerEquationsLinearFormDensity* form = new EulerEquationsLinearFormDensity(*this);
      form->wf = this->wf;
      return form;
    }
  };

  class EulerEquationsLinearFormDensityVelX : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormDensityVelX(double kappa)
      : VectorFormVol<double>(1), kappa(kappa) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * ext->fn[0]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_1_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_1_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_1_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_1_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    double kappa;
    VectorFormVol<double>* clone() { return new EulerEquationsLinearFormDensityVelX(*this); }
  };

  class EulerEquationsLinearFormDensityVelY : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormDensityVelY(double kappa) 
      : VectorFormVol<double>(2), kappa(kappa) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * ext->fn[0]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_2_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_2_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_2_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_2_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    VectorFormVol<double>* clone() { return new EulerEquationsLinearFormDensityVelY(*this); }

    double kappa;
  };

  class EulerEquationsLinearFormEnergy : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormEnergy(double kappa) 
      : VectorFormVol<double>(3), kappa(kappa) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++) 
      {
        result += wt[i] * ext->fn[0]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
          * v->dx[i];
        result += wt[i] * ext->fn[0]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_3_0<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
          * v->dy[i];
        result += wt[i] * ext->fn[1]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
          * v->dx[i];
        result += wt[i] * ext->fn[1]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_3_1<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
        result += wt[i] * ext->fn[2]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[2]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_3_2<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], ext->fn[3]->val[i]) 
          * v->dy[i];
        result += wt[i] * ext->fn[3]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_1_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dx[i];
        result += wt[i] * ext->fn[3]->val[i] 
        * (static_cast<EulerEquationsWeakFormExplicit*>(wf))->euler_fluxes->A_2_3_3<Scalar>(ext->fn[0]->val[i], ext->fn[1]->val[i], ext->fn[2]->val[i], Scalar(0)) 
          * v->dy[i];
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    VectorFormVol<double>* clone() { return new EulerEquationsLinearFormEnergy(*this); }
    double kappa;
  };

  class EulerEquationsLinearFormTime : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormTime(int i) 
      : VectorFormVol<double>(i), component_i(i) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      return int_u_v<Real, Scalar>(n, wt, ext->fn[component_i], v);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    VectorFormVol<double>* clone() { return new EulerEquationsLinearFormTime(*this); }
    // Member.
    int component_i;
  };

  class EulerEquationsLinearFormInterface : public VectorFormDG<double>
  {
  public:
    EulerEquationsLinearFormInterface(int i, double kappa) 
      : VectorFormDG<double>(i), element(i), num_flux(new StegerWarmingNumericalFlux(kappa)) {}

    ~EulerEquationsLinearFormInterface()
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0;
      double w_L[4], w_R[4];
      for (int point_i = 0; point_i < n;point_i++) 
      {
        w_L[0] = ext->fn[0]->get_val_central(point_i);
        w_R[0] = ext->fn[0]->get_val_neighbor(point_i);

        w_L[1] = ext->fn[1]->get_val_central(point_i);
        w_R[1] = ext->fn[1]->get_val_neighbor(point_i);

        w_L[2] = ext->fn[2]->get_val_central(point_i);
        w_R[2] = ext->fn[2]->get_val_neighbor(point_i);

        w_L[3] = ext->fn[3]->get_val_central(point_i);
        w_R[3] = ext->fn[3]->get_val_neighbor(point_i);

        result -= wt[point_i] * v->val[point_i] 
        * num_flux->numerical_flux_i(element, w_L, w_R, e->nx[point_i], e->ny[point_i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    VectorFormDG<double>* clone()
    {
      EulerEquationsLinearFormInterface* form = new EulerEquationsLinearFormInterface(this->i, this->num_flux->kappa);
      form->wf = this->wf;
      return form;
    }

    // Members.
    int element;
    NumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormSolidWall : public VectorFormSurf<double>
  {
  public:
    EulerEquationsLinearFormSolidWall(int i, std::string marker, double kappa) 
      : VectorFormSurf<double>(i), element(i), num_flux(new StegerWarmingNumericalFlux(kappa)) { setArea(marker); }

    ~EulerEquationsLinearFormSolidWall()
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0;
      for (int i = 0;i < n;i++) 
      {
        double w_L[4];
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];

        result -= wt[i] * v->val[i] * num_flux->numerical_flux_solid_wall_i(element, 
          w_L, e->nx[i], e->ny[i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    VectorFormSurf<double>* clone()
    {
      EulerEquationsLinearFormSolidWall* form = new EulerEquationsLinearFormSolidWall(i, areas[0], num_flux->kappa);
      form->wf = this->wf;
      return form;
    }

    // Members.
    int element;
    NumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormInlet : public VectorFormSurf<double>
  {
  public:
    EulerEquationsLinearFormInlet(int i, std::string marker, double kappa) 
      : VectorFormSurf<double>(i), element(i), num_flux(new StegerWarmingNumericalFlux(kappa)) { setArea(marker); }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, 
      ExtData<double> *ext) const 
    {
      double result = 0;
      double w_L[4], w_B[4];

      for (int i = 0;i < n;i++) 
      {
        // Left (inner) state from the previous time level solution.
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];

        w_B[0] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->rho_ext;
        w_B[1] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->rho_ext 
          * static_cast<EulerEquationsWeakFormExplicit*>(wf)->v1_ext;
        w_B[2] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->rho_ext 
          * static_cast<EulerEquationsWeakFormExplicit*>(wf)->v2_ext;
        w_B[3] = static_cast<EulerEquationsWeakFormExplicit*>(wf)->energy_ext;

        result -= wt[i] * v->val[i] * num_flux->numerical_flux_inlet_i(element, 
          w_L, w_B, e->nx[i], e->ny[i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    VectorFormSurf<double>* clone()
    {
      EulerEquationsLinearFormInlet* form = new EulerEquationsLinearFormInlet(i, areas[0], num_flux->kappa);
      form->wf = this->wf;
      return form;
    }

    // Members.
    int element;
    NumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormOutlet : public VectorFormSurf<double>
  {
  public:
    EulerEquationsLinearFormOutlet(int i, std::string marker, double kappa) : VectorFormSurf<double>(i), element(i), num_flux(new StegerWarmingNumericalFlux(kappa)) 
    { setArea(marker); }

    ~EulerEquationsLinearFormOutlet()
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0;
      double w_L[4];
      for (int i = 0;i < n;i++) 
      {
        w_L[0] = ext->fn[0]->val[i];
        w_L[1] = ext->fn[1]->val[i];
        w_L[2] = ext->fn[2]->val[i];
        w_L[3] = ext->fn[3]->val[i];
        result -= wt[i] * v->val[i] 
        * num_flux->numerical_flux_outlet_i(element, w_L, 
          static_cast<EulerEquationsWeakFormExplicit*>(wf)->pressure_ext, 
          e->nx[i], e->ny[i]);
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(20);
    }

    VectorFormSurf<double>* clone()
    {
      EulerEquationsLinearFormOutlet* form = new EulerEquationsLinearFormOutlet(i, areas[0], num_flux->kappa);
      form->wf = this->wf;
      return form;
    }

    // Members.
    int element;
    NumericalFlux* num_flux;
  };
  // Members.
  double rho_ext;
  double v1_ext;
  double v2_ext;
  double pressure_ext;
  double energy_ext;
  double tau;

  EulerFluxes* euler_fluxes;

  friend class EulerEquationsWeakFormExplicitCoupled;
  friend class EulerEquationsWeakFormExplicitMultiComponent;
  friend class EulerEquationsWeakFormSemiImplicit;
  friend class EulerEquationsWeakFormSemiImplicitCoupled;
};

class EulerEquationsWeakFormSemiImplicit : public WeakForm<double>
{
public:
  // Constructor.
  EulerEquationsWeakFormSemiImplicit(double kappa, double rho_ext, double v1_ext, double v2_ext, double pressure_ext, 
    std::string solid_wall_bottom_marker, std::string solid_wall_top_marker, std::string inlet_marker, std::string outlet_marker, 
    Solution<double>* prev_density, Solution<double>* prev_density_vel_x, Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, bool fvm_only = false, int num_of_equations = 4) :
  WeakForm<double>(num_of_equations), rho_ext(rho_ext), v1_ext(v1_ext), v2_ext(v2_ext), pressure_ext(pressure_ext), 
    energy_ext(QuantityCalculator::calc_energy(rho_ext, rho_ext * v1_ext, rho_ext * v2_ext, pressure_ext, kappa)), euler_fluxes(new EulerFluxes(kappa)), centralElemId(-1), neighbElemId(-1), uTrn(-1), vTrn(-1), iSurf(-1)
  {
    P_plus_cache = new double**[Hermes::Hermes2D::Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
    for(int i = 0; i < Hermes::Hermes2D::Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
    {
      P_plus_cache[i] = new double*[13];
      for(int j = 0; j < 13; j++)
        P_plus_cache[i][j] = new double[16];
    }

    P_minus_cache = new double**[Hermes::Hermes2D::Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
    for(int i = 0; i < Hermes::Hermes2D::Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
    {
      P_minus_cache[i] = new double*[13];
      for(int j = 0; j < 13; j++)
        P_minus_cache[i][j] = new double[16];
    }

    add_matrix_form(new EulerEquationsBilinearFormTime(0));
    add_matrix_form(new EulerEquationsBilinearFormTime(1));
    add_matrix_form(new EulerEquationsBilinearFormTime(2));
    add_matrix_form(new EulerEquationsBilinearFormTime(3));

    add_vector_form(new EulerEquationsLinearFormTime(0));
    add_vector_form(new EulerEquationsLinearFormTime(1));
    add_vector_form(new EulerEquationsLinearFormTime(2));
    add_vector_form(new EulerEquationsLinearFormTime(3));

    if(!fvm_only) 
    {
      add_matrix_form(new EulerEquationsBilinearForm(0, 0));
      add_matrix_form(new EulerEquationsBilinearForm(0, 1));
      add_matrix_form(new EulerEquationsBilinearForm(0, 2));
      add_matrix_form(new EulerEquationsBilinearForm(0, 3));
      add_matrix_form(new EulerEquationsBilinearForm(1, 0));
      add_matrix_form(new EulerEquationsBilinearForm(1, 1));
      add_matrix_form(new EulerEquationsBilinearForm(1, 2));
      add_matrix_form(new EulerEquationsBilinearForm(1, 3));
      add_matrix_form(new EulerEquationsBilinearForm(2, 0));
      add_matrix_form(new EulerEquationsBilinearForm(2, 1));
      add_matrix_form(new EulerEquationsBilinearForm(2, 2));
      add_matrix_form(new EulerEquationsBilinearForm(2, 3));
      add_matrix_form(new EulerEquationsBilinearForm(3, 0));
      add_matrix_form(new EulerEquationsBilinearForm(3, 1));
      add_matrix_form(new EulerEquationsBilinearForm(3, 2));
      add_matrix_form(new EulerEquationsBilinearForm(3, 3));
    }

    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(0, 0, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(1, 0, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(2, 0, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(3, 0, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(0, 1, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(1, 1, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(2, 1, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(3, 1, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(0, 2, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(1, 2, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(2, 2, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(3, 2, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(0, 3, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(1, 3, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(2, 3, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(3, 3, kappa));

    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 0, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 0, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 0, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 0, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 1, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 1, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 1, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 1, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 2, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 2, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 2, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 2, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 3, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 3, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 3, solid_wall_bottom_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 3, solid_wall_bottom_marker, kappa));

    if(solid_wall_bottom_marker != solid_wall_top_marker) 
    {
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 0, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 0, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 0, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 0, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 1, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 1, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 1, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 1, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 2, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 2, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 2, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 2, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 3, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 3, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 3, solid_wall_top_marker, kappa));
      add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 3, solid_wall_top_marker, kappa));
    }
    else
      Hermes::Mixins::Loggable::Static::warn("Are you sure that solid wall top and bottom markers should coincide?");

    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(0, 0, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(1, 0, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(2, 0, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(3, 0, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(0, 1, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(1, 1, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(2, 1, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(3, 1, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(0, 2, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(1, 2, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(2, 2, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(3, 2, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(0, 3, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(1, 3, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(2, 3, inlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(3, 3, inlet_marker, kappa));

    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(0, inlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(1, inlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(2, inlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(3, inlet_marker, kappa));

    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(0, 0, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(1, 0, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(2, 0, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(3, 0, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(0, 1, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(1, 1, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(2, 1, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(3, 1, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(0, 2, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(1, 2, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(2, 2, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(3, 2, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(0, 3, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(1, 3, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(2, 3, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet(3, 3, outlet_marker, kappa));

    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(0, outlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(1, outlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(2, outlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet(3, outlet_marker, kappa));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfvol.size();vector_form_i++) 
      vfvol.at(vector_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfsurf.size();vector_form_i++) 
      vfsurf.at(vector_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfDG.size();vector_form_i++) 
      vfDG.at(vector_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int matrix_form_i = 0;matrix_form_i < this->mfvol.size();matrix_form_i++) 
      mfvol.at(matrix_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int matrix_form_i = 0;matrix_form_i < this->mfsurf.size();matrix_form_i++) 
      mfsurf.at(matrix_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int matrix_form_i = 0;matrix_form_i < this->mfDG.size();matrix_form_i++) 
      mfDG.at(matrix_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
  };

  void set_time_step(double tau) 
  {
    this->tau = tau;
  }

  double get_tau() const 
  {
    return tau;
  }

  void set_stabilization(Solution<double>* prev_density, Solution<double>* prev_density_vel_x, Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, double nu_1, double nu_2) 
  {
    int mfvol_size = this->mfvol.size();
    int mfsurf_size = this->mfsurf.size();   

    add_matrix_form(new EulerEquationsFormStabilizationVol(0, nu_1));
    add_matrix_form(new EulerEquationsFormStabilizationVol(1, nu_1));
    add_matrix_form(new EulerEquationsFormStabilizationVol(2, nu_1));
    add_matrix_form(new EulerEquationsFormStabilizationVol(3, nu_1));

    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 3, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 3, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 3, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 3, nu_2));

    for(unsigned int matrix_form_i = mfvol_size;matrix_form_i < this->mfvol.size();matrix_form_i++) 
    {
      mfvol.at(matrix_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
    }

    for(unsigned int matrix_form_i = mfsurf_size;matrix_form_i < this->mfsurf.size();matrix_form_i++) 
    {
      mfsurf.at(matrix_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
    }
  }

  void set_discreteIndicator(bool* discreteIndicator)
  {
    this->discreteIndicator = discreteIndicator;
  }

protected:
  class EulerEquationsBilinearFormTime : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormTime(int i) : MatrixFormVol<double>(i, i) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormVol<double>* clone() { return new EulerEquationsBilinearFormTime(this->i); }
  };

  class EulerEquationsBilinearForm : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearForm(int i, int j) : MatrixFormVol<double>(i, j) {}
    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, 
      ExtData<double> *ext) const
    {
      double result = 0.;
      for (int point_i = 0; point_i < n;point_i++) 
      {
        switch(i)
        {
        case 0:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_0_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_0_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 1:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_0_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_0_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 2:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_0_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_0_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 3:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_0_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_0_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          }
          break;
        case 1:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_1_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_1_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 1:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_1_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_1_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 2:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_1_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_1_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 3:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_1_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_1_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          }
          break;
        case 2:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_2_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_2_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 1:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_2_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_2_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 2:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_2_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_2_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 3:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_2_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_2_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          }
          break;
        case 3:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_3_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], ext->fn[3]->val[point_i]) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_3_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], ext->fn[3]->val[point_i]) 
              * v->dy[point_i];
            break;
          case 1:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_3_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], ext->fn[3]->val[point_i]) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_3_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 2:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_3_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_3_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], ext->fn[3]->val[point_i]) 
              * v->dy[point_i];
            break;
          case 3:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_1_3_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf))->euler_fluxes->A_2_3_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          }
        }
      }

      return - result * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormVol<double>* clone()
    {
      EulerEquationsBilinearForm* form = new EulerEquationsBilinearForm(this->i, this->j);
      form->wf = this->wf;
      return form;
    }
  };

  class EulerEquationsMatrixFormSurfSemiImplicit : public MatrixFormDG<double>
  {
  public:
    EulerEquationsMatrixFormSurfSemiImplicit(int i, int j, double kappa) 
      : MatrixFormDG<double>(i, j), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) {}

    ~EulerEquationsMatrixFormSurfSemiImplicit() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double w[4];
      double result = 0.;
      if(!(u->sub_idx == static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->uTrn && 
          v->sub_idx == static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->vTrn && 
          e->id == static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->centralElemId && 
          static_cast<InterfaceGeom<double>*>(e)->get_neighbor_id() == static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->neighbElemId))
        for (int point_i = 0; point_i < n; point_i++) 
        {
          {
            w[0] = ext->fn[0]->get_val_central(point_i);
            w[1] = ext->fn[1]->get_val_central(point_i);
            w[2] = ext->fn[2]->get_val_central(point_i);
            w[3] = ext->fn[3]->get_val_central(point_i);

            double e_1_1[4] = {1, 0, 0, 0};
            double e_2_1[4] = {0, 1, 0, 0};
            double e_3_1[4] = {0, 0, 1, 0};
            double e_4_1[4] = {0, 0, 0, 1};

            num_flux->P_plus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i], w, e_1_1, e->nx[point_i], e->ny[point_i]);
            num_flux->P_plus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i] + 4, w, e_2_1, e->nx[point_i], e->ny[point_i]);
            num_flux->P_plus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i] + 8, w, e_3_1, e->nx[point_i], e->ny[point_i]);
            num_flux->P_plus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i] + 12, w, e_4_1, e->nx[point_i], e->ny[point_i]);

            w[0] = ext->fn[0]->get_val_neighbor(point_i);
            w[1] = ext->fn[1]->get_val_neighbor(point_i);
            w[2] = ext->fn[2]->get_val_neighbor(point_i);
            w[3] = ext->fn[3]->get_val_neighbor(point_i);

            double e_1_2[4] = {1, 0, 0, 0};
            double e_2_2[4] = {0, 1, 0, 0};
            double e_3_2[4] = {0, 0, 1, 0};
            double e_4_2[4] = {0, 0, 0, 1};

            num_flux->P_minus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i], w, e_1_2, e->nx[point_i], e->ny[point_i]);
            num_flux->P_minus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i] + 4, w, e_2_2, e->nx[point_i], e->ny[point_i]);
            num_flux->P_minus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i] + 8, w, e_3_2, e->nx[point_i], e->ny[point_i]);
            num_flux->P_minus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i] + 12, w, e_4_2, e->nx[point_i], e->ny[point_i]);

            static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->uTrn = u->sub_idx;
            static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->vTrn = v->sub_idx;
            static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->centralElemId = e->id;
            static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->neighbElemId = static_cast<InterfaceGeom<double>*>(e)->get_neighbor_id();
          }
        }
      for (int point_i = 0; point_i < n; point_i++) 
      {
        if(i == 0 && j == 0)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][0] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][0] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 0 && j == 1)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][4] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][4] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 0 && j == 2)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][8] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][8] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 0 && j == 3)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][12] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][12] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }

        if(i == 1 && j == 0)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][1] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][1] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 1 && j == 1)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][5] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][5] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 1 && j == 2)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][9] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][9] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 1 && j == 3)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][13] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][13] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }

        if(i == 2 && j == 0)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][2] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][2] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 2 && j == 1)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][6] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][6] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 2 && j == 2)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][10] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][10] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 2 && j == 3)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][14] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][14] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }

        if(i == 3 && j == 0)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][3] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][3] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 3 && j == 1)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][7] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][7] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 3 && j == 2)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][11] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][11] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 3 && j == 3)
        {
          result += wt[point_i] * (static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][15] * u->get_val_central(point_i) + static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_minus_cache[omp_get_thread_num()][point_i][15] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormDG<double>* clone() 
    { 
      EulerEquationsMatrixFormSurfSemiImplicit* form = new EulerEquationsMatrixFormSurfSemiImplicit(this->i, this->j, this->num_flux->kappa);
      form->wf = this->wf;
      return form;
    }

    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsMatrixFormSemiImplicitInletOutlet : public MatrixFormSurf<double>
  {
  public:
    EulerEquationsMatrixFormSemiImplicitInletOutlet(int i, int j, 
      std::string marker, double kappa) 
      : MatrixFormSurf<double>(i, j), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) { setArea(marker); }

    ~EulerEquationsMatrixFormSemiImplicitInletOutlet() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;

      if(!(u->sub_idx == static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->uTrn && 
          v->sub_idx == static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->vTrn && 
          e->id == static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->centralElemId && e->isurf == static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->iSurf))
        for (int point_i = 0; point_i < n; point_i++) 
        {
          double w_B[4], w_L[4], eigenvalues[4], alpha[4], q_ji_star[4], beta[4], q_ji[4], w_ji[4];

          // Inner state.
          w_L[0] = ext->fn[0]->val[point_i];
          w_L[1] = ext->fn[1]->val[point_i];
          w_L[2] = ext->fn[2]->val[point_i];
          w_L[3] = ext->fn[3]->val[point_i];

          // Transformation of the inner state to the local coordinates.
          num_flux->Q(num_flux->get_q(), w_L, e->nx[point_i], e->ny[point_i]);

          // Initialize the matrices.
          double T[4][4];
          double T_inv[4][4];
          for(unsigned int ai = 0; ai < 4; ai++) 
          {
            for(unsigned int aj = 0; aj < 4; aj++) 
            {
              T[ai][aj] = 0.0;
              T_inv[ai][aj] = 0.0;
            }
            alpha[ai] = 0;
            beta[ai] = 0;
            q_ji[ai] = 0;
            w_ji[ai] = 0;
            eigenvalues[ai] = 0;
          }

          // Calculate Lambda^-.
          num_flux->Lambda(eigenvalues);
          num_flux->T_1(T);
          num_flux->T_2(T);
          num_flux->T_3(T);
          num_flux->T_4(T);
          num_flux->T_inv_1(T_inv);
          num_flux->T_inv_2(T_inv);
          num_flux->T_inv_3(T_inv);
          num_flux->T_inv_4(T_inv);

          // "Prescribed" boundary state.
          w_B[0] = static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->rho_ext;
          w_B[1] = static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->rho_ext 
            * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->v1_ext;
          w_B[2] = static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->rho_ext 
            * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->v2_ext;
          w_B[3] = static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->energy_ext;

          num_flux->Q(q_ji_star, w_B, e->nx[point_i], e->ny[point_i]);

          for(unsigned int ai = 0; ai < 4; ai++)
            for(unsigned int aj = 0; aj < 4; aj++)
              alpha[ai] += T_inv[ai][aj] * num_flux->get_q()[aj];

          for(unsigned int bi = 0; bi < 4; bi++)
            for(unsigned int bj = 0; bj < 4; bj++)
              beta[bi] += T_inv[bi][bj] * q_ji_star[bj];

          for(unsigned int si = 0; si< 4; si++)
            for(unsigned int sj = 0; sj < 4; sj++)
              if(eigenvalues[sj] < 0)
                q_ji[si] += beta[sj] * T[si][sj];
              else
                q_ji[si] += alpha[sj] * T[si][sj];

          num_flux->Q_inv(w_ji, q_ji, e->nx[point_i], e->ny[point_i]);

          double P_minus[4];

          double w_temp[4];
          w_temp[0] = (w_ji[0] + w_L[0]) / 2;
          w_temp[1] = (w_ji[1] + w_L[1]) / 2;
          w_temp[2] = (w_ji[2] + w_L[2]) / 2;
          w_temp[3] = (w_ji[3] + w_L[3]) / 2;

          double e_1[4] = {1, 0, 0, 0};
          double e_2[4] = {0, 1, 0, 0};
          double e_3[4] = {0, 0, 1, 0};
          double e_4[4] = {0, 0, 0, 1};

          num_flux->P_plus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i], w_temp, e_1, e->nx[point_i], e->ny[point_i]);
          num_flux->P_plus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i] + 4, w_temp, e_2, e->nx[point_i], e->ny[point_i]);
          num_flux->P_plus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i] + 8, w_temp, e_3, e->nx[point_i], e->ny[point_i]);
          num_flux->P_plus(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i] + 12, w_temp, e_4, e->nx[point_i], e->ny[point_i]);

          static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->uTrn = u->sub_idx;
          static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->vTrn = v->sub_idx;
          static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->centralElemId = e->id;
          static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->iSurf = e->isurf;
        }

      for (int point_i = 0; point_i < n; point_i++) 
      {
        if(i == 0 && j == 0)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 0 && j == 1)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][4] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 0 && j == 2)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][8] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 0 && j == 3)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][12] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }

        if(i == 1 && j == 0)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 1 && j == 1)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][5] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 1 && j == 2)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][9] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 1 && j == 3)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][13] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }

        if(i == 2 && j == 0)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 2 && j == 1)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][6] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 2 && j == 2)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][10] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 2 && j == 3)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][14] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }

        if(i == 3 && j == 0)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 3 && j == 1)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][7] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 3 && j == 2)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][11] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 3 && j == 3)
        {
          result += wt[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->P_plus_cache[omp_get_thread_num()][point_i][15] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
      }

      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormSurf<double>* clone() 
    { 
      EulerEquationsMatrixFormSemiImplicitInletOutlet* form = new EulerEquationsMatrixFormSemiImplicitInletOutlet(this->i, this->j, this->areas[0], this->num_flux->kappa);
      form->wf = this->wf;
      return form;
    }

    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsVectorFormSemiImplicitInletOutlet : public VectorFormSurf<double>
  {
  public:
    EulerEquationsVectorFormSemiImplicitInletOutlet(int i, std::string marker, double kappa) 
      : VectorFormSurf<double>(i), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) {setArea(marker); }

    ~EulerEquationsVectorFormSemiImplicitInletOutlet() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;

      double w_B[4], w_L[4], eigenvalues[4], alpha[4], q_ji_star[4], beta[4], q_ji[4], w_ji[4];

      for (int point_i = 0; point_i < n; point_i++) 
      {
        // Inner state.
        w_L[0] = ext->fn[0]->val[point_i];
        w_L[1] = ext->fn[1]->val[point_i];
        w_L[2] = ext->fn[2]->val[point_i];
        w_L[3] = ext->fn[3]->val[point_i];

        // Transformation of the inner state to the local coordinates.
        num_flux->Q(num_flux->get_q(), w_L, e->nx[point_i], e->ny[point_i]);

        // Initialize the matrices.
        double T[4][4];
        double T_inv[4][4];
        for(unsigned int ai = 0; ai < 4; ai++) 
        {
          for(unsigned int aj = 0; aj < 4; aj++) 
          {
            T[ai][aj] = 0.0;
            T_inv[ai][aj] = 0.0;
          }
          alpha[ai] = 0;
          beta[ai] = 0;
          q_ji[ai] = 0;
          w_ji[ai] = 0;
          eigenvalues[ai] = 0;
        }

        // Calculate Lambda^-.
        num_flux->Lambda(eigenvalues);
        num_flux->T_1(T);
        num_flux->T_2(T);
        num_flux->T_3(T);
        num_flux->T_4(T);
        num_flux->T_inv_1(T_inv);
        num_flux->T_inv_2(T_inv);
        num_flux->T_inv_3(T_inv);
        num_flux->T_inv_4(T_inv);

        // "Prescribed" boundary state.
        w_B[0] = static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->rho_ext;
        w_B[1] = static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->rho_ext 
          * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->v1_ext;
        w_B[2] = static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->rho_ext 
          * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->v2_ext;
        w_B[3] = static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->energy_ext;

        num_flux->Q(q_ji_star, w_B, e->nx[point_i], e->ny[point_i]);

        for(unsigned int ai = 0; ai < 4; ai++)
          for(unsigned int aj = 0; aj < 4; aj++)
            alpha[ai] += T_inv[ai][aj] * num_flux->get_q()[aj];

        for(unsigned int bi = 0; bi < 4; bi++)
          for(unsigned int bj = 0; bj < 4; bj++)
            beta[bi] += T_inv[bi][bj] * q_ji_star[bj];

        for(unsigned int si = 0; si< 4; si++)
          for(unsigned int sj = 0; sj < 4; sj++)
            if(eigenvalues[sj] < 0)
              q_ji[si] += beta[sj] * T[si][sj];
            else
              q_ji[si] += alpha[sj] * T[si][sj];

        num_flux->Q_inv(w_ji, q_ji, e->nx[point_i], e->ny[point_i]);

        double P_minus[4];

        double w_temp[4];
        w_temp[0] = (w_ji[0] + w_L[0]) / 2;
        w_temp[1] = (w_ji[1] + w_L[1]) / 2;
        w_temp[2] = (w_ji[2] + w_L[2]) / 2;
        w_temp[3] = (w_ji[3] + w_L[3]) / 2;

        num_flux->P_minus(P_minus, w_temp, w_ji, e->nx[point_i], e->ny[point_i]);

        result += wt[point_i] * (P_minus[i]) * v->val[point_i];
      }

      return - result * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    VectorFormSurf<double>* clone() 
    { 
      EulerEquationsVectorFormSemiImplicitInletOutlet* form = new EulerEquationsVectorFormSemiImplicitInletOutlet(this->i, this->areas[0], this->num_flux->kappa);
      form->wf = this->wf;
      return form;
    }

    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormTime : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormTime(int i) 
      : VectorFormVol<double>(i) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
      ExtData<double> *ext) const 
    {
      return int_u_v<double, double>(n, wt, ext->fn[this->i], v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    VectorFormVol<double>* clone() { return new EulerEquationsLinearFormTime(this->i); }
  };

  class EulerEquationsMatrixFormSolidWall : public MatrixFormSurf<double>
  {
  public:
    EulerEquationsMatrixFormSolidWall(int i, int j, std::string marker, double kappa)
      : MatrixFormSurf<double>(i, j), kappa(kappa) {setArea(marker);}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;

      for (int point_i = 0; point_i < n; point_i++) 
      {
        double rho = ext->fn[0]->val[point_i];
        double v_1 = ext->fn[1]->val[point_i] / rho;
        double v_2 = ext->fn[2]->val[point_i] / rho;

        double P[4][4];
        for(unsigned int P_i = 0; P_i < 4; P_i++)
          for(unsigned int P_j = 0; P_j < 4; P_j++)
            P[P_i][P_j] = 0.0;

        P[1][0] = (kappa - 1) * (v_1 * v_1 + v_2 * v_2) * e->nx[point_i] / 2;
        P[1][1] = (kappa - 1) * (-v_1) * e->nx[point_i];
        P[1][2] = (kappa - 1) * (-v_2) * e->nx[point_i];
        P[1][3] = (kappa - 1) * e->nx[point_i];

        P[2][0] = (kappa - 1) * (v_1 * v_1 + v_2 * v_2) * e->ny[point_i] / 2;
        P[2][1] = (kappa - 1) * (-v_1) * e->ny[point_i];
        P[2][2] = (kappa - 1) * (-v_2) * e->ny[point_i];
        P[2][3] = (kappa - 1) * e->ny[point_i];

        if(i == 0 && j == 0)
        {
          result += wt[point_i] * P[0][0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 0 && j == 1)
        {
          result += wt[point_i] * P[0][1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 0 && j == 2)
        {
          result += wt[point_i] * P[0][2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 0 && j == 3)
        {
          result += wt[point_i] * P[0][3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }

        if(i == 1 && j == 0)
        {
          result += wt[point_i] * P[1][0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 1 && j == 1)
        {
          result += wt[point_i] * P[1][1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 1 && j == 2)
        {
          result += wt[point_i] * P[1][2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 1 && j == 3)
        {
          result += wt[point_i] * P[1][3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }

        if(i == 2 && j == 0)
        {
          result += wt[point_i] * P[2][0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 2 && j == 1)
        {
          result += wt[point_i] * P[2][1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 2 && j == 2)
        {
          result += wt[point_i] * P[2][2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 2 && j == 3)
        {
          result += wt[point_i] * P[2][3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }

        if(i == 3 && j == 0)
        {
          result += wt[point_i] * P[3][0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 3 && j == 1)
        {
          result += wt[point_i] * P[3][1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 3 && j == 2)
        {
          result += wt[point_i] * P[3][2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
        if(i == 3 && j == 3)
        {
          result += wt[point_i] * P[3][3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->get_tau();
          continue;
        }
      }

      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormSurf<double>* clone() 
    {
      EulerEquationsMatrixFormSolidWall* form = new EulerEquationsMatrixFormSolidWall(this->i, this->j, this->areas[0], this->kappa);
      form->wf = this->wf;
      return form;
    }

    // Members.
    double kappa;
  };

  class EulerEquationsFormStabilizationVol : public MatrixFormVol<double>
  {
  public:
    EulerEquationsFormStabilizationVol(int i, double nu_1) : MatrixFormVol<double>(i, i), nu_1(nu_1) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;
      if(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->discreteIndicator[e->id]) 
        result = int_grad_u_grad_v<double, double>(n, wt, u, v) * nu_1 * e->diam;
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
    }
    MatrixFormVol<double>* clone() 
    {
      EulerEquationsFormStabilizationVol* form = new EulerEquationsFormStabilizationVol(this->i, nu_1);
      form->wf = this->wf;
      return form;
    }
  private:
    double nu_1;
  };

  class EulerEquationsFormStabilizationSurf : public MatrixFormDG<double>
  {
  public:
    EulerEquationsFormStabilizationSurf(int i, int j, double nu_2) 
      : MatrixFormDG<double>(i, j), nu_2(nu_2) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;

      if(static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->discreteIndicator[e->id] && static_cast<EulerEquationsWeakFormSemiImplicit*>(wf)->discreteIndicator[e->get_neighbor_id()])
        for (int i = 0; i < n;i++)
          result += wt[i] * (u->get_val_central(i) - u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * nu_2;

      return result;
    }

    MatrixFormDG<double>* clone() 
    {
      EulerEquationsFormStabilizationSurf* form = new EulerEquationsFormStabilizationSurf(this->i, this->j, nu_2);
      form->wf = this->wf;
      return form;
    }

    double nu_2;
  };

  // Members.
  double rho_ext;
  double v1_ext;
  double v2_ext;
  double pressure_ext;
  double energy_ext;
  double tau;
  EulerFluxes* euler_fluxes;
  bool* discreteIndicator;

  // Cache for numerical flux.
  double*** P_plus_cache;
  double*** P_minus_cache;
  int centralElemId;
  int neighbElemId;
  uint64_t uTrn;
  uint64_t vTrn;
  int iSurf;
};

class EulerEquationsWeakFormSemiImplicitTwoInflows : public WeakForm<double>
{
public:
  // Constructor.
  EulerEquationsWeakFormSemiImplicitTwoInflows(double kappa, double rho_ext1, double v1_ext1, double v2_ext1, double pressure_ext1, double rho_ext2, double v1_ext2, double v2_ext2, double pressure_ext2, 
    std::string wall_marker, std::string inlet_marker1, std::string inlet_marker2, std::string outlet_marker,
    Solution<double>* prev_density, Solution<double>* prev_density_vel_x, Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, bool fvm_only = false, int num_of_equations = 4) :
  WeakForm<double>(num_of_equations), 
    rho_ext1(rho_ext1), v1_ext1(v1_ext1), v2_ext1(v2_ext1), pressure_ext1(pressure_ext1), 
    energy_ext1(QuantityCalculator::calc_energy(rho_ext1, rho_ext1* v1_ext1, rho_ext1 * v2_ext1, pressure_ext1, kappa)), 
    rho_ext2(rho_ext2), v1_ext2(v1_ext2), v2_ext2(v2_ext2), pressure_ext2(pressure_ext2), 
    energy_ext2(QuantityCalculator::calc_energy(rho_ext2, rho_ext2 * v1_ext2, rho_ext2 * v2_ext2, pressure_ext2, kappa)), 
    euler_fluxes(new EulerFluxes(kappa))
  {
    add_matrix_form(new EulerEquationsBilinearFormTime(0));
    add_matrix_form(new EulerEquationsBilinearFormTime(1));
    add_matrix_form(new EulerEquationsBilinearFormTime(2));
    add_matrix_form(new EulerEquationsBilinearFormTime(3));

    add_vector_form(new EulerEquationsLinearFormTime(0));
    add_vector_form(new EulerEquationsLinearFormTime(1));
    add_vector_form(new EulerEquationsLinearFormTime(2));
    add_vector_form(new EulerEquationsLinearFormTime(3));

    if(!fvm_only) 
    {
      add_matrix_form(new EulerEquationsBilinearForm(0, 0));
      add_matrix_form(new EulerEquationsBilinearForm(0, 1));
      add_matrix_form(new EulerEquationsBilinearForm(0, 2));
      add_matrix_form(new EulerEquationsBilinearForm(0, 3));
      add_matrix_form(new EulerEquationsBilinearForm(1, 0));
      add_matrix_form(new EulerEquationsBilinearForm(1, 1));
      add_matrix_form(new EulerEquationsBilinearForm(1, 2));
      add_matrix_form(new EulerEquationsBilinearForm(1, 3));
      add_matrix_form(new EulerEquationsBilinearForm(2, 0));
      add_matrix_form(new EulerEquationsBilinearForm(2, 1));
      add_matrix_form(new EulerEquationsBilinearForm(2, 2));
      add_matrix_form(new EulerEquationsBilinearForm(2, 3));
      add_matrix_form(new EulerEquationsBilinearForm(3, 0));
      add_matrix_form(new EulerEquationsBilinearForm(3, 1));
      add_matrix_form(new EulerEquationsBilinearForm(3, 2));
      add_matrix_form(new EulerEquationsBilinearForm(3, 3));
    }

    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(0, 0, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(1, 0, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(2, 0, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(3, 0, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(0, 1, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(1, 1, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(2, 1, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(3, 1, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(0, 2, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(1, 2, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(2, 2, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(3, 2, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(0, 3, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(1, 3, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(2, 3, kappa));
    add_matrix_form_DG(new EulerEquationsMatrixFormSurfSemiImplicit(3, 3, kappa));

    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 0, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 0, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 0, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 0, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 1, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 1, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 1, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 1, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 2, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 2, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 2, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 2, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(0, 3, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(1, 3, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(2, 3, wall_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSolidWall(3, 3, wall_marker, kappa));

    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(0, 0, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(1, 0, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(2, 0, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(3, 0, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(0, 1, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(1, 1, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(2, 1, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(3, 1, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(0, 2, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(1, 2, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(2, 2, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(3, 2, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(0, 3, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(1, 3, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(2, 3, inlet_marker1, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(3, 3, inlet_marker1, kappa));

    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet1(0, inlet_marker1, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet1(1, inlet_marker1, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet1(2, inlet_marker1, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet1(3, inlet_marker1, kappa));

    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(0, 0, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(1, 0, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(2, 0, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(3, 0, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(0, 1, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(1, 1, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(2, 1, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(3, 1, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(0, 2, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(1, 2, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(2, 2, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(3, 2, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(0, 3, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(1, 3, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(2, 3, inlet_marker2, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet2(3, 3, inlet_marker2, kappa));

    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet2(0, inlet_marker2, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet2(1, inlet_marker2, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet2(2, inlet_marker2, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet2(3, inlet_marker2, kappa));


    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(0, 0, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(1, 0, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(2, 0, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(3, 0, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(0, 1, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(1, 1, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(2, 1, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(3, 1, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(0, 2, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(1, 2, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(2, 2, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(3, 2, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(0, 3, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(1, 3, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(2, 3, outlet_marker, kappa));
    add_matrix_form_surf(new EulerEquationsMatrixFormSemiImplicitInletOutlet1(3, 3, outlet_marker, kappa));

    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet1(0, outlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet1(1, outlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet1(2, outlet_marker, kappa));
    add_vector_form_surf(new EulerEquationsVectorFormSemiImplicitInletOutlet1(3, outlet_marker, kappa));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfvol.size();vector_form_i++) 
      vfvol.at(vector_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfsurf.size();vector_form_i++) 
      vfsurf.at(vector_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int vector_form_i = 0;vector_form_i < this->vfDG.size();vector_form_i++) 
      vfDG.at(vector_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int matrix_form_i = 0;matrix_form_i < this->mfvol.size();matrix_form_i++) 
      mfvol.at(matrix_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int matrix_form_i = 0;matrix_form_i < this->mfsurf.size();matrix_form_i++) 
      mfsurf.at(matrix_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int matrix_form_i = 0;matrix_form_i < this->mfDG.size();matrix_form_i++) 
      mfDG.at(matrix_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
  };

  void set_time_step(double tau) 
  {
    this->tau = tau;
  }

  double get_tau() const 
  {
    return tau;
  }

  void set_stabilization(Solution<double>* prev_density, Solution<double>* prev_density_vel_x, Solution<double>* prev_density_vel_y, Solution<double>* prev_energy, double nu_1, double nu_2) 
  {
    int mfvol_size = this->mfvol.size();
    int mfsurf_size = this->mfsurf.size();   

    add_matrix_form(new EulerEquationsFormStabilizationVol(0, nu_1));
    add_matrix_form(new EulerEquationsFormStabilizationVol(1, nu_1));
    add_matrix_form(new EulerEquationsFormStabilizationVol(2, nu_1));
    add_matrix_form(new EulerEquationsFormStabilizationVol(3, nu_1));

    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(0, 3, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(1, 3, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(2, 3, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 0, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 1, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 2, nu_2));
    add_matrix_form_DG(new EulerEquationsFormStabilizationSurf(3, 3, nu_2));

    for(unsigned int matrix_form_i = mfvol_size;matrix_form_i < this->mfvol.size();matrix_form_i++) 
    {
      mfvol.at(matrix_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
    }

    for(unsigned int matrix_form_i = mfsurf_size;matrix_form_i < this->mfsurf.size();matrix_form_i++) 
    {
      mfsurf.at(matrix_form_i)->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
    }
  }

  void set_discreteIndicator(bool* discreteIndicator)
  {
    this->discreteIndicator = discreteIndicator;
  }

protected:
  class EulerEquationsBilinearFormTime : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearFormTime(int i) : MatrixFormVol<double>(i, i) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    MatrixFormVol<double>* clone() { return new EulerEquationsBilinearFormTime(this->i); }
  };

  class EulerEquationsBilinearForm : public MatrixFormVol<double>
  {
  public:
    EulerEquationsBilinearForm(int i, int j) : MatrixFormVol<double>(i, j) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, 
      ExtData<double> *ext) const
    {
      double result = 0.;
      for (int point_i = 0; point_i < n;point_i++) 
      {
        switch(i)
        {
        case 0:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_0_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_0_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 1:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_0_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_0_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 2:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_0_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_0_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 3:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_0_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_0_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          }
          break;
        case 1:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_1_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_1_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 1:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_1_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_1_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 2:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_1_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_1_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 3:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_1_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_1_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          }
          break;
        case 2:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_2_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_2_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 1:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_2_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_2_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 2:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_2_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_2_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 3:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_2_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_2_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          }
          break;
        case 3:
          switch(j)
          {
          case 0:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_3_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], ext->fn[3]->val[point_i]) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_3_0<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], ext->fn[3]->val[point_i]) 
              * v->dy[point_i];
            break;
          case 1:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_3_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], ext->fn[3]->val[point_i]) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_3_1<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          case 2:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_3_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_3_2<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], ext->fn[3]->val[point_i]) 
              * v->dy[point_i];
            break;
          case 3:
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_1_3_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dx[point_i];
            result += wt[point_i] * u->val[point_i]
            * (static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf))->euler_fluxes->A_2_3_3<double>(ext->fn[0]->val[point_i], ext->fn[1]->val[point_i], ext->fn[2]->val[point_i], 0) 
              * v->dy[point_i];
            break;
          }
        }
      }

      return - result * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }
    MatrixFormVol<double>* clone()
    {
      EulerEquationsBilinearForm* form = new EulerEquationsBilinearForm(this->i, this->j);
      form->wf = this->wf;
      return form;
    }
  };

  class EulerEquationsMatrixFormSurfSemiImplicit : public MatrixFormDG<double>
  {
  public:
    EulerEquationsMatrixFormSurfSemiImplicit(int i, int j, double kappa) 
      : MatrixFormDG<double>(i, j), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) {}

    virtual ~EulerEquationsMatrixFormSurfSemiImplicit() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double w[4];
      double result = 0.;
      for (int point_i = 0; point_i < n; point_i++) 
      {
        w[0] = ext->fn[0]->get_val_central(point_i);
        w[1] = ext->fn[1]->get_val_central(point_i);
        w[2] = ext->fn[2]->get_val_central(point_i);
        w[3] = ext->fn[3]->get_val_central(point_i);

        double e_1_1[4] = {1, 0, 0, 0};
        double e_2_1[4] = {0, 1, 0, 0};
        double e_3_1[4] = {0, 0, 1, 0};
        double e_4_1[4] = {0, 0, 0, 1};

        double P_plus_1[4];
        double P_plus_2[4];
        double P_plus_3[4];
        double P_plus_4[4];

        num_flux->P_plus(P_plus_1, w, e_1_1, e->nx[point_i], e->ny[point_i]);
        num_flux->P_plus(P_plus_2, w, e_2_1, e->nx[point_i], e->ny[point_i]);
        num_flux->P_plus(P_plus_3, w, e_3_1, e->nx[point_i], e->ny[point_i]);
        num_flux->P_plus(P_plus_4, w, e_4_1, e->nx[point_i], e->ny[point_i]);

        w[0] = ext->fn[0]->get_val_neighbor(point_i);
        w[1] = ext->fn[1]->get_val_neighbor(point_i);
        w[2] = ext->fn[2]->get_val_neighbor(point_i);
        w[3] = ext->fn[3]->get_val_neighbor(point_i);

        double e_1_2[4] = {1, 0, 0, 0};
        double e_2_2[4] = {0, 1, 0, 0};
        double e_3_2[4] = {0, 0, 1, 0};
        double e_4_2[4] = {0, 0, 0, 1};

        double P_minus_1[4];
        double P_minus_2[4];
        double P_minus_3[4];
        double P_minus_4[4];

        num_flux->P_minus(P_minus_1, w, e_1_2, e->nx[point_i], e->ny[point_i]);
        num_flux->P_minus(P_minus_2, w, e_2_2, e->nx[point_i], e->ny[point_i]);
        num_flux->P_minus(P_minus_3, w, e_3_2, e->nx[point_i], e->ny[point_i]);
        num_flux->P_minus(P_minus_4, w, e_4_2, e->nx[point_i], e->ny[point_i]);

        if(i == 0 && j == 0)
          result += wt[point_i] * (P_plus_1[0] * u->get_val_central(point_i) + P_minus_1[0] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 1)
          result += wt[point_i] * (P_plus_2[0] * u->get_val_central(point_i) + P_minus_2[0] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 2)
          result += wt[point_i] * (P_plus_3[0] * u->get_val_central(point_i) + P_minus_3[0] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 3)
          result += wt[point_i] * (P_plus_4[0] * u->get_val_central(point_i) + P_minus_4[0] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 1 && j == 0)
          result += wt[point_i] * (P_plus_1[1] * u->get_val_central(point_i) + P_minus_1[1] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 1)
          result += wt[point_i] * (P_plus_2[1] * u->get_val_central(point_i) + P_minus_2[1] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 2)
          result += wt[point_i] * (P_plus_3[1] * u->get_val_central(point_i) + P_minus_3[1] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 3)
          result += wt[point_i] * (P_plus_4[1] * u->get_val_central(point_i) + P_minus_4[1] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 2 && j == 0)
          result += wt[point_i] * (P_plus_1[2] * u->get_val_central(point_i) + P_minus_1[2] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 1)
          result += wt[point_i] * (P_plus_2[2] * u->get_val_central(point_i) + P_minus_2[2] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 2)
          result += wt[point_i] * (P_plus_3[2] * u->get_val_central(point_i) + P_minus_3[2] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 3)
          result += wt[point_i] * (P_plus_4[2] * u->get_val_central(point_i) + P_minus_4[2] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 3 && j == 0)
          result += wt[point_i] * (P_plus_1[3] * u->get_val_central(point_i) + P_minus_1[3] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 1)
          result += wt[point_i] * (P_plus_2[3] * u->get_val_central(point_i) + P_minus_2[3] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 2)
          result += wt[point_i] * (P_plus_3[3] * u->get_val_central(point_i) + P_minus_3[3] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 3)
          result += wt[point_i] * (P_plus_4[3] * u->get_val_central(point_i) + P_minus_4[3] * u->get_val_neighbor(point_i)) * (v->get_val_central(point_i) - v->get_val_neighbor(point_i)) * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
      }

      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormDG<double>* clone()
    {
      EulerEquationsMatrixFormSurfSemiImplicit* form = new EulerEquationsMatrixFormSurfSemiImplicit(this->i, this->j, this->num_flux->kappa);
      form->wf = this->wf;
      return form;
    }
    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsMatrixFormSemiImplicitInletOutlet1 : public MatrixFormSurf<double>
  {
  public:
    EulerEquationsMatrixFormSemiImplicitInletOutlet1(int i, int j, 
      std::string marker, double kappa) 
      : MatrixFormSurf<double>(i, j), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) { setArea(marker);}

    virtual ~EulerEquationsMatrixFormSemiImplicitInletOutlet1() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;

      double w_B[4], w_L[4], eigenvalues[4], alpha[4], q_ji_star[4], beta[4], q_ji[4], w_ji[4];

      for (int point_i = 0; point_i < n; point_i++) 
      {
        // Inner state.
        w_L[0] = ext->fn[0]->val[point_i];
        w_L[1] = ext->fn[1]->val[point_i];
        w_L[2] = ext->fn[2]->val[point_i];
        w_L[3] = ext->fn[3]->val[point_i];

        // Transformation of the inner state to the local coordinates.
        num_flux->Q(num_flux->get_q(), w_L, e->nx[point_i], e->ny[point_i]);

        // Initialize the matrices.
        double T[4][4];
        double T_inv[4][4];
        for(unsigned int ai = 0; ai < 4; ai++) 
        {
          for(unsigned int aj = 0; aj < 4; aj++) 
          {
            T[ai][aj] = 0.0;
            T_inv[ai][aj] = 0.0;
          }
          alpha[ai] = 0;
          beta[ai] = 0;
          q_ji[ai] = 0;
          w_ji[ai] = 0;
          eigenvalues[ai] = 0;
        }

        // Calculate Lambda^-.
        num_flux->Lambda(eigenvalues);
        num_flux->T_1(T);
        num_flux->T_2(T);
        num_flux->T_3(T);
        num_flux->T_4(T);
        num_flux->T_inv_1(T_inv);
        num_flux->T_inv_2(T_inv);
        num_flux->T_inv_3(T_inv);
        num_flux->T_inv_4(T_inv);

        // "Prescribed" boundary state.
        w_B[0] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext1;
        w_B[1] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext1 
          * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->v1_ext1;
        w_B[2] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext1
          * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->v2_ext1;
        w_B[3] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->energy_ext1;

        num_flux->Q(q_ji_star, w_B, e->nx[point_i], e->ny[point_i]);

        for(unsigned int ai = 0; ai < 4; ai++)
          for(unsigned int aj = 0; aj < 4; aj++)
            alpha[ai] += T_inv[ai][aj] * num_flux->get_q()[aj];

        for(unsigned int bi = 0; bi < 4; bi++)
          for(unsigned int bj = 0; bj < 4; bj++)
            beta[bi] += T_inv[bi][bj] * q_ji_star[bj];

        for(unsigned int si = 0; si< 4; si++)
          for(unsigned int sj = 0; sj < 4; sj++)
            if(eigenvalues[sj] < 0)
              q_ji[si] += beta[sj] * T[si][sj];
            else
              q_ji[si] += alpha[sj] * T[si][sj];

        num_flux->Q_inv(w_ji, q_ji, e->nx[point_i], e->ny[point_i]);

        double P_minus[4];

        double w_temp[4];
        w_temp[0] = (w_ji[0] + w_L[0]) / 2;
        w_temp[1] = (w_ji[1] + w_L[1]) / 2;
        w_temp[2] = (w_ji[2] + w_L[2]) / 2;
        w_temp[3] = (w_ji[3] + w_L[3]) / 2;

        double e_1[4] = {1, 0, 0, 0};
        double e_2[4] = {0, 1, 0, 0};
        double e_3[4] = {0, 0, 1, 0};
        double e_4[4] = {0, 0, 0, 1};

        double P_plus_1[4];
        double P_plus_2[4];
        double P_plus_3[4];
        double P_plus_4[4];

        num_flux->P_plus(P_plus_1, w_temp, e_1, e->nx[point_i], e->ny[point_i]);
        num_flux->P_plus(P_plus_2, w_temp, e_2, e->nx[point_i], e->ny[point_i]);
        num_flux->P_plus(P_plus_3, w_temp, e_3, e->nx[point_i], e->ny[point_i]);
        num_flux->P_plus(P_plus_4, w_temp, e_4, e->nx[point_i], e->ny[point_i]);

        if(i == 0 && j == 0)
          result += wt[point_i] * P_plus_1[0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 1)
          result += wt[point_i] * P_plus_2[0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 2)
          result += wt[point_i] * P_plus_3[0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 3)
          result += wt[point_i] * P_plus_4[0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 1 && j == 0)
          result += wt[point_i] * P_plus_1[1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 1)
          result += wt[point_i] * P_plus_2[1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 2)
          result += wt[point_i] * P_plus_3[1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 3)
          result += wt[point_i] * P_plus_4[1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 2 && j == 0)
          result += wt[point_i] * P_plus_1[2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 1)
          result += wt[point_i] * P_plus_2[2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 2)
          result += wt[point_i] * P_plus_3[2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 3)
          result += wt[point_i] * P_plus_4[2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 3 && j == 0)
          result += wt[point_i] * P_plus_1[3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 1)
          result += wt[point_i] * P_plus_2[3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 2)
          result += wt[point_i] * P_plus_3[3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 3)
          result += wt[point_i] * P_plus_4[3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
      }

      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormSurf<double>* clone()
    {
      EulerEquationsMatrixFormSemiImplicitInletOutlet1* form = new EulerEquationsMatrixFormSemiImplicitInletOutlet1(this->i, this->j, this->areas[0], this->num_flux->kappa);
      form->wf = this->wf;
      return form;
    }
    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsVectorFormSemiImplicitInletOutlet1 : public VectorFormSurf<double>
  {
  public:
    EulerEquationsVectorFormSemiImplicitInletOutlet1(int i, std::string marker, double kappa) 
      : VectorFormSurf<double>(i), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) { setArea(marker);}

    ~EulerEquationsVectorFormSemiImplicitInletOutlet1() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;

      double w_B[4], w_L[4], eigenvalues[4], alpha[4], q_ji_star[4], beta[4], q_ji[4], w_ji[4];

      for (int point_i = 0; point_i < n; point_i++) 
      {
        // Inner state.
        w_L[0] = ext->fn[0]->val[point_i];
        w_L[1] = ext->fn[1]->val[point_i];
        w_L[2] = ext->fn[2]->val[point_i];
        w_L[3] = ext->fn[3]->val[point_i];

        // Transformation of the inner state to the local coordinates.
        num_flux->Q(num_flux->get_q(), w_L, e->nx[point_i], e->ny[point_i]);

        // Initialize the matrices.
        double T[4][4];
        double T_inv[4][4];
        for(unsigned int ai = 0; ai < 4; ai++) 
        {
          for(unsigned int aj = 0; aj < 4; aj++) 
          {
            T[ai][aj] = 0.0;
            T_inv[ai][aj] = 0.0;
          }
          alpha[ai] = 0;
          beta[ai] = 0;
          q_ji[ai] = 0;
          w_ji[ai] = 0;
          eigenvalues[ai] = 0;
        }

        // Calculate Lambda^-.
        num_flux->Lambda(eigenvalues);
        num_flux->T_1(T);
        num_flux->T_2(T);
        num_flux->T_3(T);
        num_flux->T_4(T);
        num_flux->T_inv_1(T_inv);
        num_flux->T_inv_2(T_inv);
        num_flux->T_inv_3(T_inv);
        num_flux->T_inv_4(T_inv);

        // "Prescribed" boundary state.
        w_B[0] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext1;
        w_B[1] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext1
          * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->v1_ext1;
        w_B[2] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext1 
          * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->v2_ext1;
        w_B[3] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->energy_ext1;

        num_flux->Q(q_ji_star, w_B, e->nx[point_i], e->ny[point_i]);

        for(unsigned int ai = 0; ai < 4; ai++)
          for(unsigned int aj = 0; aj < 4; aj++)
            alpha[ai] += T_inv[ai][aj] * num_flux->get_q()[aj];

        for(unsigned int bi = 0; bi < 4; bi++)
          for(unsigned int bj = 0; bj < 4; bj++)
            beta[bi] += T_inv[bi][bj] * q_ji_star[bj];

        for(unsigned int si = 0; si< 4; si++)
          for(unsigned int sj = 0; sj < 4; sj++)
            if(eigenvalues[sj] < 0)
              q_ji[si] += beta[sj] * T[si][sj];
            else
              q_ji[si] += alpha[sj] * T[si][sj];

        num_flux->Q_inv(w_ji, q_ji, e->nx[point_i], e->ny[point_i]);

        double P_minus[4];

        double w_temp[4];
        w_temp[0] = (w_ji[0] + w_L[0]) / 2;
        w_temp[1] = (w_ji[1] + w_L[1]) / 2;
        w_temp[2] = (w_ji[2] + w_L[2]) / 2;
        w_temp[3] = (w_ji[3] + w_L[3]) / 2;

        num_flux->P_minus(P_minus, w_temp, w_ji, e->nx[point_i], e->ny[point_i]);

        result += wt[point_i] * (P_minus[i]) * v->val[point_i];
      }

      return - result * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    VectorFormSurf<double>* clone()
    {
      EulerEquationsVectorFormSemiImplicitInletOutlet1* form = new EulerEquationsVectorFormSemiImplicitInletOutlet1(this->i, this->areas[0], this->num_flux->kappa);
      form->wf = this->wf;
      return form;
    }
    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsMatrixFormSemiImplicitInletOutlet2 : public MatrixFormSurf<double>
  {
  public:
    EulerEquationsMatrixFormSemiImplicitInletOutlet2(int i, int j, 
      std::string marker, double kappa) 
      : MatrixFormSurf<double>(i, j), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) { setArea(marker);}

    virtual ~EulerEquationsMatrixFormSemiImplicitInletOutlet2() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;

      double w_B[4], w_L[4], eigenvalues[4], alpha[4], q_ji_star[4], beta[4], q_ji[4], w_ji[4];

      for (int point_i = 0; point_i < n; point_i++) 
      {
        // Inner state.
        w_L[0] = ext->fn[0]->val[point_i];
        w_L[1] = ext->fn[1]->val[point_i];
        w_L[2] = ext->fn[2]->val[point_i];
        w_L[3] = ext->fn[3]->val[point_i];

        // Transformation of the inner state to the local coordinates.
        num_flux->Q(num_flux->get_q(), w_L, e->nx[point_i], e->ny[point_i]);

        // Initialize the matrices.
        double T[4][4];
        double T_inv[4][4];
        for(unsigned int ai = 0; ai < 4; ai++) 
        {
          for(unsigned int aj = 0; aj < 4; aj++) 
          {
            T[ai][aj] = 0.0;
            T_inv[ai][aj] = 0.0;
          }
          alpha[ai] = 0;
          beta[ai] = 0;
          q_ji[ai] = 0;
          w_ji[ai] = 0;
          eigenvalues[ai] = 0;
        }

        // Calculate Lambda^-.
        num_flux->Lambda(eigenvalues);
        num_flux->T_1(T);
        num_flux->T_2(T);
        num_flux->T_3(T);
        num_flux->T_4(T);
        num_flux->T_inv_1(T_inv);
        num_flux->T_inv_2(T_inv);
        num_flux->T_inv_3(T_inv);
        num_flux->T_inv_4(T_inv);

        // "Prescribed" boundary state.
        w_B[0] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext2;
        w_B[1] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext2 
          * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->v1_ext2;
        w_B[2] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext2 
          * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->v2_ext2;
        w_B[3] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->energy_ext2;

        num_flux->Q(q_ji_star, w_B, e->nx[point_i], e->ny[point_i]);

        for(unsigned int ai = 0; ai < 4; ai++)
          for(unsigned int aj = 0; aj < 4; aj++)
            alpha[ai] += T_inv[ai][aj] * num_flux->get_q()[aj];

        for(unsigned int bi = 0; bi < 4; bi++)
          for(unsigned int bj = 0; bj < 4; bj++)
            beta[bi] += T_inv[bi][bj] * q_ji_star[bj];

        for(unsigned int si = 0; si< 4; si++)
          for(unsigned int sj = 0; sj < 4; sj++)
            if(eigenvalues[sj] < 0)
              q_ji[si] += beta[sj] * T[si][sj];
            else
              q_ji[si] += alpha[sj] * T[si][sj];

        num_flux->Q_inv(w_ji, q_ji, e->nx[point_i], e->ny[point_i]);

        double P_minus[4];

        double w_temp[4];
        w_temp[0] = (w_ji[0] + w_L[0]) / 2;
        w_temp[1] = (w_ji[1] + w_L[1]) / 2;
        w_temp[2] = (w_ji[2] + w_L[2]) / 2;
        w_temp[3] = (w_ji[3] + w_L[3]) / 2;

        double e_1[4] = {1, 0, 0, 0};
        double e_2[4] = {0, 1, 0, 0};
        double e_3[4] = {0, 0, 1, 0};
        double e_4[4] = {0, 0, 0, 1};

        double P_plus_1[4];
        double P_plus_2[4];
        double P_plus_3[4];
        double P_plus_4[4];

        num_flux->P_plus(P_plus_1, w_temp, e_1, e->nx[point_i], e->ny[point_i]);
        num_flux->P_plus(P_plus_2, w_temp, e_2, e->nx[point_i], e->ny[point_i]);
        num_flux->P_plus(P_plus_3, w_temp, e_3, e->nx[point_i], e->ny[point_i]);
        num_flux->P_plus(P_plus_4, w_temp, e_4, e->nx[point_i], e->ny[point_i]);

        if(i == 0 && j == 0)
          result += wt[point_i] * P_plus_1[0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 1)
          result += wt[point_i] * P_plus_2[0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 2)
          result += wt[point_i] * P_plus_3[0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 3)
          result += wt[point_i] * P_plus_4[0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 1 && j == 0)
          result += wt[point_i] * P_plus_1[1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 1)
          result += wt[point_i] * P_plus_2[1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 2)
          result += wt[point_i] * P_plus_3[1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 3)
          result += wt[point_i] * P_plus_4[1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 2 && j == 0)
          result += wt[point_i] * P_plus_1[2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 1)
          result += wt[point_i] * P_plus_2[2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 2)
          result += wt[point_i] * P_plus_3[2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 3)
          result += wt[point_i] * P_plus_4[2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 3 && j == 0)
          result += wt[point_i] * P_plus_1[3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 1)
          result += wt[point_i] * P_plus_2[3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 2)
          result += wt[point_i] * P_plus_3[3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 3)
          result += wt[point_i] * P_plus_4[3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
      }

      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormSurf<double>* clone()
    {
      EulerEquationsMatrixFormSemiImplicitInletOutlet2* form = new EulerEquationsMatrixFormSemiImplicitInletOutlet2(this->i, this->j, this->areas[0], this->num_flux->kappa);
      form->wf = this->wf;
      return form;
    }
    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsVectorFormSemiImplicitInletOutlet2 : public VectorFormSurf<double>
  {
  public:
    EulerEquationsVectorFormSemiImplicitInletOutlet2(int i, std::string marker, double kappa) 
      : VectorFormSurf<double>(i), 
      num_flux(new StegerWarmingNumericalFlux(kappa)) { setArea(marker);}

    ~EulerEquationsVectorFormSemiImplicitInletOutlet2() 
    {
      delete num_flux;
    }

    double value(int n, double *wt, Func<double> *u_ext[],
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;

      double w_B[4], w_L[4], eigenvalues[4], alpha[4], q_ji_star[4], beta[4], q_ji[4], w_ji[4];

      for (int point_i = 0; point_i < n; point_i++) 
      {
        // Inner state.
        w_L[0] = ext->fn[0]->val[point_i];
        w_L[1] = ext->fn[1]->val[point_i];
        w_L[2] = ext->fn[2]->val[point_i];
        w_L[3] = ext->fn[3]->val[point_i];

        // Transformation of the inner state to the local coordinates.
        num_flux->Q(num_flux->get_q(), w_L, e->nx[point_i], e->ny[point_i]);

        // Initialize the matrices.
        double T[4][4];
        double T_inv[4][4];
        for(unsigned int ai = 0; ai < 4; ai++) 
        {
          for(unsigned int aj = 0; aj < 4; aj++) 
          {
            T[ai][aj] = 0.0;
            T_inv[ai][aj] = 0.0;
          }
          alpha[ai] = 0;
          beta[ai] = 0;
          q_ji[ai] = 0;
          w_ji[ai] = 0;
          eigenvalues[ai] = 0;
        }

        // Calculate Lambda^-.
        num_flux->Lambda(eigenvalues);
        num_flux->T_1(T);
        num_flux->T_2(T);
        num_flux->T_3(T);
        num_flux->T_4(T);
        num_flux->T_inv_1(T_inv);
        num_flux->T_inv_2(T_inv);
        num_flux->T_inv_3(T_inv);
        num_flux->T_inv_4(T_inv);

        // "Prescribed" boundary state.
        w_B[0] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext2;
        w_B[1] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext2
          * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->v1_ext2;
        w_B[2] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->rho_ext2 
          * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->v2_ext2;
        w_B[3] = static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->energy_ext2;

        num_flux->Q(q_ji_star, w_B, e->nx[point_i], e->ny[point_i]);

        for(unsigned int ai = 0; ai < 4; ai++)
          for(unsigned int aj = 0; aj < 4; aj++)
            alpha[ai] += T_inv[ai][aj] * num_flux->get_q()[aj];

        for(unsigned int bi = 0; bi < 4; bi++)
          for(unsigned int bj = 0; bj < 4; bj++)
            beta[bi] += T_inv[bi][bj] * q_ji_star[bj];

        for(unsigned int si = 0; si< 4; si++)
          for(unsigned int sj = 0; sj < 4; sj++)
            if(eigenvalues[sj] < 0)
              q_ji[si] += beta[sj] * T[si][sj];
            else
              q_ji[si] += alpha[sj] * T[si][sj];

        num_flux->Q_inv(w_ji, q_ji, e->nx[point_i], e->ny[point_i]);

        double P_minus[4];

        double w_temp[4];
        w_temp[0] = (w_ji[0] + w_L[0]) / 2;
        w_temp[1] = (w_ji[1] + w_L[1]) / 2;
        w_temp[2] = (w_ji[2] + w_L[2]) / 2;
        w_temp[3] = (w_ji[3] + w_L[3]) / 2;

        num_flux->P_minus(P_minus, w_temp, w_ji, e->nx[point_i], e->ny[point_i]);

        result += wt[point_i] * (P_minus[i]) * v->val[point_i];
      }

      return - result * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    VectorFormSurf<double>* clone()
    {
      EulerEquationsVectorFormSemiImplicitInletOutlet2* form = new EulerEquationsVectorFormSemiImplicitInletOutlet2(this->i, this->areas[0], this->num_flux->kappa);
      form->wf = this->wf;
      return form;
    }
    StegerWarmingNumericalFlux* num_flux;
  };

  class EulerEquationsLinearFormTime : public VectorFormVol<double>
  {
  public:
    EulerEquationsLinearFormTime(int i) 
      : VectorFormVol<double>(i) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
      ExtData<double> *ext) const 
    {
      return int_u_v<double, double>(n, wt, ext->fn[this->i], v);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }
    VectorFormVol<double>* clone() { return new EulerEquationsLinearFormTime(this->i); }
  };

  class EulerEquationsMatrixFormSolidWall : public MatrixFormSurf<double>
  {
  public:
    EulerEquationsMatrixFormSolidWall(int i, int j, std::string marker, double kappa)
      : MatrixFormSurf<double>(i, j), kappa(kappa) {setArea(marker);}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;

      for (int point_i = 0; point_i < n; point_i++) 
      {
        double rho = ext->fn[0]->val[point_i];
        double v_1 = ext->fn[1]->val[point_i] / rho;
        double v_2 = ext->fn[2]->val[point_i] / rho;

        double P[4][4];
        for(unsigned int P_i = 0; P_i < 4; P_i++)
          for(unsigned int P_j = 0; P_j < 4; P_j++)
            P[P_i][P_j] = 0.0;

        P[1][0] = (kappa - 1) * (v_1 * v_1 + v_2 * v_2) * e->nx[point_i] / 2;
        P[1][1] = (kappa - 1) * (-v_1) * e->nx[point_i];
        P[1][2] = (kappa - 1) * (-v_2) * e->nx[point_i];
        P[1][3] = (kappa - 1) * e->nx[point_i];

        P[2][0] = (kappa - 1) * (v_1 * v_1 + v_2 * v_2) * e->ny[point_i] / 2;
        P[2][1] = (kappa - 1) * (-v_1) * e->ny[point_i];
        P[2][2] = (kappa - 1) * (-v_2) * e->ny[point_i];
        P[2][3] = (kappa - 1) * e->ny[point_i];

        if(i == 0 && j == 0)
          result += wt[point_i] * P[0][0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 1)
          result += wt[point_i] * P[0][1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 2)
          result += wt[point_i] * P[0][2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 0 && j == 3)
          result += wt[point_i] * P[0][3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 1 && j == 0)
          result += wt[point_i] * P[1][0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 1)
          result += wt[point_i] * P[1][1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 2)
          result += wt[point_i] * P[1][2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 1 && j == 3)
          result += wt[point_i] * P[1][3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 2 && j == 0)
          result += wt[point_i] * P[2][0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 1)
          result += wt[point_i] * P[2][1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 2)
          result += wt[point_i] * P[2][2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 2 && j == 3)
          result += wt[point_i] * P[2][3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();

        if(i == 3 && j == 0)
          result += wt[point_i] * P[3][0] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 1)
          result += wt[point_i] * P[3][1] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 2)
          result += wt[point_i] * P[3][2] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
        if(i == 3 && j == 3)
          result += wt[point_i] * P[3][3] * u->val[point_i] * v->val[point_i] * static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->get_tau();
      }

      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormSurf<double>* clone()
    {
      EulerEquationsMatrixFormSolidWall* form = new EulerEquationsMatrixFormSolidWall(this->i, this->j, this->areas[0], this->kappa);
      form->wf = this->wf;
      return form;
    }
    // Members.
    double kappa;
  };

  class EulerEquationsFormStabilizationVol : public MatrixFormVol<double>
  {
  public:
    EulerEquationsFormStabilizationVol(int i, double nu_1) : MatrixFormVol<double>(i, i), nu_1(nu_1) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      double result_i = 0.;
      if(static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->discreteIndicator[e->id]) 
        return int_grad_u_grad_v<double, double>(n, wt, u, v) * nu_1 * e->diam;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormVol<double>* clone() { return new EulerEquationsFormStabilizationVol(this->i, this->nu_1); }
  private:
    double nu_1;
  };

  class EulerEquationsFormStabilizationSurf : public MatrixFormDG<double>
  {
  public:
    EulerEquationsFormStabilizationSurf(int i, int j, double nu_2) 
      : MatrixFormDG<double>(i, j), nu_2(nu_2) {}

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
    {
      double result = 0.;

      if(static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->discreteIndicator[e->id] && static_cast<EulerEquationsWeakFormSemiImplicitTwoInflows*>(wf)->discreteIndicator[e->get_neighbor_id()])
        for (int i = 0;i < n;i++)
          result += wt[i] * (u->get_val_central(i) - u->get_val_neighbor(i)) * (v->get_val_central(i) - v->get_val_neighbor(i)) * nu_2;

      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
      Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormDG<double>* clone() { return new EulerEquationsFormStabilizationSurf(this->i, this->j, this->nu_2); }

    double nu_2;
  };

  // Members.
  double rho_ext1;
  double v1_ext1;
  double v2_ext1;
  double pressure_ext1;
  double energy_ext1;

  double rho_ext2;
  double v1_ext2;
  double v2_ext2;
  double pressure_ext2;
  double energy_ext2;

  double tau;
  EulerFluxes* euler_fluxes;
  bool* discreteIndicator;
};

class EulerEquationsWeakFormExplicitCoupled : public EulerEquationsWeakFormExplicit
{
public:
  // Constructor.
  EulerEquationsWeakFormExplicitCoupled(double kappa,
    double rho_ext, double v1_ext, double v2_ext,
    double pressure_ext, std::string solid_wall_marker_bottom, std::string solid_wall_marker_top,
    std::string inlet_marker, std::string outlet_marker,
    Hermes::vector<std::string> natural_bc_concentration_markers,
    Solution<double>* prev_density, Solution<double>* prev_density_vel_x,
    Solution<double>* prev_density_vel_y, Solution<double>* prev_energy,
    Solution<double>* prev_concentration, double epsilon, bool fvm_only = false)
    : EulerEquationsWeakFormExplicit(kappa, rho_ext, v1_ext, v2_ext, pressure_ext,
    solid_wall_marker_bottom, solid_wall_marker_top, inlet_marker,
    outlet_marker, prev_density, prev_density_vel_x,
    prev_density_vel_y, prev_energy, fvm_only, 5) 
  {
    add_matrix_form(new EulerEquationsWeakFormExplicit::EulerEquationsBilinearFormTime(4));

    add_vector_form(new VectorFormConcentrationAdvectionDiffusion(4, epsilon));
    vfvol.back()->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy, prev_concentration));

    for(unsigned int i = 0;i < natural_bc_concentration_markers.size();i++) 
    {
      add_vector_form_surf(new VectorFormConcentrationNatural(4, natural_bc_concentration_markers[i]));
      vfsurf.back()->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy, prev_concentration));
    }

    EulerEquationsWeakFormExplicit::EulerEquationsLinearFormTime* vector_form_time = new EulerEquationsWeakFormExplicit::EulerEquationsLinearFormTime(4);
    vector_form_time->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy, prev_concentration));
    add_vector_form(vector_form_time);
  };

  // Destructor.
  ~EulerEquationsWeakFormExplicitCoupled() {};
protected:

  class VectorFormConcentrationAdvectionDiffusion : public VectorFormVol<double>
  {
  public:
    VectorFormConcentrationAdvectionDiffusion(int i, double epsilon) 
      : VectorFormVol<double>(i), epsilon(epsilon) {}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Scalar result = Scalar(0);
      Real h_e = e->diam;
      Real s_c = Real(0.9);

      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];
      Func<Real>* concentration_prev = ext->fn[4];

      for (int point_i=0; point_i < n; point_i++) 
      {
        Scalar v_1 = density_vel_x_prev->val[point_i] / density_prev->val[point_i];
        Scalar v_2 = density_vel_y_prev->val[point_i] / density_prev->val[point_i];

        result += wt[point_i] * (epsilon * (concentration_prev->dx[point_i]*v->dx[point_i] + concentration_prev->dy[point_i]*v->dy[point_i])
          - (v_1 * concentration_prev->val[point_i] * v->dx[point_i] + v_2 * concentration_prev->val[point_i] * v->dy[point_i]));

        result += 100 * wt[point_i] * ((v_1 * concentration_prev->dx[point_i] + v_2 * concentration_prev->dy[point_i]) - (epsilon * (concentration_prev->dx[point_i]*concentration_prev->dx[point_i] + concentration_prev->dy[point_i]*concentration_prev->dy[point_i])))
          * (v_1 * v->dx[point_i] + v_2 * v->dy[point_i]) * h_e / (2 * std::sqrt(v_1*v_1 + v_2*v_2));

        /*

        Real R_squared = Hermes::pow(v_1 * u->dx[i] + v_2 * u->dy[i], 2.);
        Real R = Hermes::sqrt(R_squared); //This just does fabs(b1 * u->dx[i] + b2 * u->dy[i]); but it can be parsed
        result += wt[i] * s_c * 0.5 * h_e * R * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) / (Hermes::sqrt(Hermes::pow(u->dx[i], 2) + Hermes::pow(u->dy[i], 2)) + 1.e-8);

        Scalar b_norm = Hermes::sqrt(v_1 * v_1 + v_2 * v_2);
        Real tau = 1. / Hermes::sqrt( 9 * Hermes::pow(4 * epsilon / Hermes::pow(h_e, 2), 2) + Hermes::pow(2 * b_norm / h_e, 2));
        result += wt[i] * tau * (-v_1 * v->dx[i] - v_2 * v->dy[i] + epsilon * v->laplace[i]) * (-v_1 * u->dx[i] - v_2 * u->dy[i] + epsilon * u->laplace[i]);
        */
      }
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    VectorFormVol<double>* clone() 
    {
      VectorFormConcentrationAdvectionDiffusion* form = new VectorFormConcentrationAdvectionDiffusion(this->i, this->epsilon); 
      form->wf = this->wf;
      return form;
    }

    // Member.
    double epsilon;
  };

  class VectorFormConcentrationNatural : public VectorFormSurf<double>
  {
  public:
    VectorFormConcentrationNatural(int i, std::string marker) 
      : VectorFormSurf<double>(i) {setArea(marker);}

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];
      Func<Real>* concentration_prev = ext->fn[4];

      Scalar result = Scalar(0);
      for (int i = 0;i < n;i++)
        result += wt[i] * v->val[i] * concentration_prev->val[i]
      * (density_vel_x_prev->val[i] * e->nx[i] + density_vel_y_prev->val[i] * e->ny[i])
        / density_prev->val[i];
      // (OR: for inlet/outlet) result += wt[i] * v->val[i] * concentration_prev->val[i] 
      //      * (V1_EXT * e->nx[i] + V2_EXT * e->ny[i]);
      return result * static_cast<EulerEquationsWeakFormExplicit*>(wf)->get_tau();

    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    VectorFormSurf<double>* clone() 
    {
      VectorFormConcentrationNatural* form = new VectorFormConcentrationNatural(this->i, this->areas[0]);
      form->wf = this->wf;
      return form;
    }
  };
};

class EulerEquationsWeakFormSemiImplicitCoupled : public EulerEquationsWeakFormSemiImplicit
{
public:
  // Constructor.
  EulerEquationsWeakFormSemiImplicitCoupled(double kappa,
    double rho_ext, double v1_ext, double v2_ext,
    double pressure_ext, std::string solid_wall_marker_bottom, std::string solid_wall_marker_top,
    std::string inlet_marker, std::string outlet_marker,
    Hermes::vector<std::string> natural_bc_concentration_markers,
    Solution<double>* prev_density, Solution<double>* prev_density_vel_x,
    Solution<double>* prev_density_vel_y, Solution<double>* prev_energy,
    Solution<double>* prev_concentration, double epsilon, bool fvm_only = false)
    : EulerEquationsWeakFormSemiImplicit(kappa, rho_ext, v1_ext, v2_ext, pressure_ext,
    solid_wall_marker_bottom, solid_wall_marker_top, inlet_marker,
    outlet_marker, prev_density, prev_density_vel_x,
    prev_density_vel_y, prev_energy, fvm_only, 5) 
  {
    add_matrix_form(new EulerEquationsWeakFormSemiImplicit::EulerEquationsBilinearFormTime(4));

    add_matrix_form(new MatrixFormConcentrationAdvectionDiffusion(4, 4, epsilon));
    mfvol.back()->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));

    for(unsigned int i = 0;i < natural_bc_concentration_markers.size();i++) 
    {
      add_matrix_form_surf(new MatrixFormConcentrationNatural(4, 4, natural_bc_concentration_markers[i]));
      mfsurf.back()->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy));
    }

    EulerEquationsWeakFormSemiImplicit::EulerEquationsLinearFormTime* vector_form_time = new EulerEquationsWeakFormSemiImplicit::EulerEquationsLinearFormTime(4);
    vector_form_time->setExt(Hermes::vector<MeshFunction<double>*>(prev_density, prev_density_vel_x, prev_density_vel_y, prev_energy, prev_concentration));
    add_vector_form(vector_form_time);
  };

  // Destructor.
  ~EulerEquationsWeakFormSemiImplicitCoupled() {};
protected:

  class MatrixFormConcentrationAdvectionDiffusion : public MatrixFormVol<double>
  {
  public:
    MatrixFormConcentrationAdvectionDiffusion(int i, int j, double epsilon) 
      : MatrixFormVol<double>(i, j), epsilon(epsilon) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Scalar result = Scalar(0);
      Real h_e = e->diam;
      Real s_c = Real(0.9);

      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];

      for (int point_i=0; point_i < n; point_i++) 
      {
        Scalar v_1 = density_vel_x_prev->val[point_i] / density_prev->val[point_i];
        Scalar v_2 = density_vel_y_prev->val[point_i] / density_prev->val[point_i];

        result += wt[point_i] * (epsilon * (u->dx[point_i]*v->dx[point_i] + u->dy[point_i]*v->dy[point_i])
          - (v_1 * u->dx[point_i] * v->val[point_i] + v_2 * u->dy[point_i] * v->val[point_i]));

        //result += 1. * wt[point_i] * ((v_1 * u->dx[point_i] + v_2 * u->dy[point_i]) - (epsilon * (u->dx[point_i]*u->dx[point_i] + u->dy[point_i]*u->dy[point_i])))
        //  * (v_1 * v->dx[point_i] + v_2 * v->dy[point_i]) * h_e / (2 * std::sqrt(v_1*v_1 + v_2*v_2));

        /*

        Real R_squared = Hermes::pow(v_1 * u->dx[i] + v_2 * u->dy[i], 2.);
        Real R = Hermes::sqrt(R_squared); //This just does fabs(b1 * u->dx[i] + b2 * u->dy[i]); but it can be parsed
        result += wt[i] * s_c * 0.5 * h_e * R * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) / (Hermes::sqrt(Hermes::pow(u->dx[i], 2) + Hermes::pow(u->dy[i], 2)) + 1.e-8);

        Scalar b_norm = Hermes::sqrt(v_1 * v_1 + v_2 * v_2);
        Real tau = 1. / Hermes::sqrt( 9 * Hermes::pow(4 * epsilon / Hermes::pow(h_e, 2), 2) + Hermes::pow(2 * b_norm / h_e, 2));
        result += wt[i] * tau * (-v_1 * v->dx[i] - v_2 * v->dy[i] + epsilon * v->laplace[i]) * (-v_1 * u->dx[i] - v_2 * u->dy[i] + epsilon * u->laplace[i]);
        */
      }
      return result * static_cast<EulerEquationsWeakFormSemiImplicitCoupled*>(wf)->get_tau();
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormVol<double>* clone() 
    { 
      MatrixFormConcentrationAdvectionDiffusion* form = new MatrixFormConcentrationAdvectionDiffusion(this->i, this->j, this->epsilon);
      form->wf = this->wf;
      return form;
    }

    // Member.
    double epsilon;
  };

  class MatrixFormConcentrationNatural : public MatrixFormSurf<double>
  {
  public:
    MatrixFormConcentrationNatural(int i, int j, std::string marker) 
      : MatrixFormSurf<double>(i, j) {setArea(marker);}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
      Geom<Real> *e, ExtData<Scalar> *ext) const 
    {
      Func<Real>* density_prev = ext->fn[0];
      Func<Real>* density_vel_x_prev = ext->fn[1];
      Func<Real>* density_vel_y_prev = ext->fn[2];

      return Scalar(0.0);
      Scalar result = Scalar(0);
      for (int point_i = 0; point_i < n; point_i++)
        result += wt[i] * v->val[i] * u->val[i] * (density_vel_x_prev->val[i] * e->nx[i] + density_vel_y_prev->val[i] * e->ny[i])
        / density_prev->val[i];
      return result * static_cast<EulerEquationsWeakFormSemiImplicitCoupled*>(wf)->get_tau();
    }

    double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
      Geom<double> *e, ExtData<double> *ext) const 
    {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
      ExtData<Ord> *ext) const 
    {
      return Ord(24);
    }

    MatrixFormSurf<double>* clone()
    { 
      MatrixFormConcentrationNatural* form = new MatrixFormConcentrationNatural(this->i, this->j, this->areas[0]);
      form->wf = this->wf;
      return form;
    }
  };
};