#include "definitions.h"

// The pressure head is raised by H_OFFSET
// so that the initial condition can be taken
// as the zero vector. Note: the resulting
// pressure head will also be greater than the
// true one by this offset.
double H_OFFSET = 1000;

/* Custom non-constant Dirichlet condition */

EssentialBoundaryCondition<double>::EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const
{
  return EssentialBoundaryCondition<double>::BC_FUNCTION;
}

double CustomEssentialBCNonConst::value(double x, double y, double n_x, double n_y,
  double t_x, double t_y) const
{
  return x*(100. - x) / 2.5 * y / 100 - 1000. + H_OFFSET;
}

/* Custom weak forms */

CustomWeakFormRichardsIE::CustomWeakFormRichardsIE(double time_step, MeshFunctionSharedPtr<double>  h_time_prev, ConstitutiveRelations* constitutive) : WeakForm<double>(1), constitutive(constitutive)
{
  // Jacobian volumetric part.
  CustomJacobianFormVol* jac_form_vol = new CustomJacobianFormVol(0, 0, time_step);
  jac_form_vol->set_ext(h_time_prev);
  add_matrix_form(jac_form_vol);

  // Residual - volumetric.
  CustomResidualFormVol* res_form_vol = new CustomResidualFormVol(0, time_step);
  res_form_vol->set_ext(h_time_prev);
  add_vector_form(res_form_vol);
}

double CustomWeakFormRichardsIE::CustomJacobianFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
  Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext[0];
  for (int i = 0; i < n; i++)
  {
    double h_val_i = h_prev_newton->val[i] - H_OFFSET;
    result += wt[i] * (static_cast<CustomWeakFormRichardsIE*>(wf)->constitutive->dCdh(h_val_i) * u->val[i] * (h_prev_newton->val[i] - h_prev_time->val[i])
      * v->val[i] + static_cast<CustomWeakFormRichardsIE*>(wf)->constitutive->C(h_val_i) * u->val[i] * v->val[i]
      + static_cast<CustomWeakFormRichardsIE*>(wf)->constitutive->dKdh(h_val_i) * u->val[i] * (h_prev_newton->dx[i] * v->dx[i]
      + h_prev_newton->dy[i] * v->dy[i]) * time_step
      + static_cast<CustomWeakFormRichardsIE*>(wf)->constitutive->K(h_val_i) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) * time_step
      - static_cast<CustomWeakFormRichardsIE*>(wf)->constitutive->ddKdhh(h_val_i) * u->val[i] * h_prev_newton->dy[i] * v->val[i] * time_step
      - static_cast<CustomWeakFormRichardsIE*>(wf)->constitutive->dKdh(h_val_i) * u->dy[i] * v->val[i] * time_step
      );
  }
  return result;
}

Ord CustomWeakFormRichardsIE::CustomJacobianFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
  Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
{
  return Ord(10);
}

MatrixFormVol<double>* CustomWeakFormRichardsIE::CustomJacobianFormVol::clone() const
{
  return new CustomJacobianFormVol(*this);
}

double CustomWeakFormRichardsIE::CustomResidualFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
  Func<double>* *ext) const
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext[0];
  for (int i = 0; i < n; i++)
  {
    double h_val_i = h_prev_newton->val[i] - H_OFFSET;
    result += wt[i] * (static_cast<CustomWeakFormRichardsIE*>(wf)->constitutive->C(h_val_i) * (h_val_i - (h_prev_time->val[i] - H_OFFSET)) * v->val[i]
      + static_cast<CustomWeakFormRichardsIE*>(wf)->constitutive->K(h_val_i) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i]) * time_step
      - static_cast<CustomWeakFormRichardsIE*>(wf)->constitutive->dKdh(h_val_i) * h_prev_newton->dy[i] * v->val[i] * time_step
      );
  }
  return result;
}

Ord CustomWeakFormRichardsIE::CustomResidualFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const
{
  return Ord(10);
}

VectorFormVol<double>* CustomWeakFormRichardsIE::CustomResidualFormVol::clone() const
{
  return new CustomResidualFormVol(*this);
}