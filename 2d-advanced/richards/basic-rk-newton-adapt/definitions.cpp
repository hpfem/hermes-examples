#include "definitions.h"

// The pressure head is raised by H_OFFSET
// so that the initial condition can be taken
// as the zero vector. Note: the resulting
// pressure head will also be greater than the
// true one by this offset.
double H_OFFSET = 1000;

/* Custom non-constant Dirichlet condition */

EssentialBCValueType CustomEssentialBCNonConst::get_value_type() const
{
  return BC_FUNCTION;
}

double CustomEssentialBCNonConst::value(double x, double y) const
{
  return x*(100. - x) / 2.5 * y / 100 - 1000. + H_OFFSET;
}

/* Custom weak forms */

CustomWeakFormRichardsRK::CustomWeakFormRichardsRK(ConstitutiveRelations* constitutive) : WeakForm<double>(1), constitutive(constitutive)
{
  // Jacobian volumetric part.
  CustomJacobianFormVol* jac_form_vol = new CustomJacobianFormVol(0, 0, constitutive);
  add_matrix_form(jac_form_vol);

  // Residual - volumetric.
  CustomResidualFormVol* res_form_vol = new CustomResidualFormVol(0, constitutive);
  add_vector_form(res_form_vol);
}

double CustomWeakFormRichardsRK::CustomJacobianFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
  Func<double> *v, GeomVol<double> *e, Func<double>* *ext) const
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
  {
    double h_val_i = h_prev_newton->val[i] - H_OFFSET;

    double C2 = constitutive->C(h_val_i) * constitutive->C(h_val_i);
    double a1_1 = (constitutive->dKdh(h_val_i) * constitutive->C(h_val_i) - constitutive->dCdh(h_val_i) * constitutive->K(h_val_i)) / C2;
    double a1_2 = constitutive->K(h_val_i) / constitutive->C(h_val_i);

    double a2_1 = ((constitutive->dKdh(h_val_i) * constitutive->dCdh(h_val_i) + constitutive->K(h_val_i) * constitutive->ddCdhh(h_val_i)) * C2
      - 2 * constitutive->K(h_val_i) * constitutive->C(h_val_i) * constitutive->dCdh(h_val_i) * constitutive->dCdh(h_val_i)) / (C2 * C2);
    double a2_2 = 2 * constitutive->K(h_val_i) * constitutive->dCdh(h_val_i) / C2;

    double a3_1 = (constitutive->ddKdhh(h_val_i) * constitutive->C(h_val_i) - constitutive->dKdh(h_val_i) * constitutive->dCdh(h_val_i)) / C2;
    double a3_2 = constitutive->dKdh(h_val_i) / constitutive->C(h_val_i);

    result += wt[i] * (-a1_1 * u->val[i] * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
      - a1_2 * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
      + a2_1 * u->val[i] * v->val[i] * (h_prev_newton->dx[i] * h_prev_newton->dx[i]
      + h_prev_newton->dy[i] * h_prev_newton->dy[i])
      + a2_2 * v->val[i] * (u->dx[i] * h_prev_newton->dx[i]
      + u->dy[i] * h_prev_newton->dy[i])
      + a3_1 * u->val[i] * v->val[i] * h_prev_newton->dy[i]
      + a3_2 * v->val[i] * u->dy[i]
      );
  }
  return result;
}

Ord CustomWeakFormRichardsRK::CustomJacobianFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
  Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const
{
  return Ord(10);
}

MatrixFormVol<double>* CustomWeakFormRichardsRK::CustomJacobianFormVol::clone() const
{
  return new CustomJacobianFormVol(*this);
}

double CustomWeakFormRichardsRK::CustomResidualFormVol::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomVol<double> *e,
  Func<double>* *ext) const
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
  {
    double h_val_i = h_prev_newton->val[i] - H_OFFSET;
    double r1 = (constitutive->K(h_val_i) / constitutive->C(h_val_i));
    double r2 = constitutive->K(h_val_i) * constitutive->dCdh(h_val_i) / (constitutive->C(h_val_i) * constitutive->C(h_val_i));
    double r3 = constitutive->dKdh(h_val_i) / constitutive->C(h_val_i);

    result += wt[i] * (-r1 * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
      + r2 * v->val[i] * (h_prev_newton->dx[i] * h_prev_newton->dx[i]
      + h_prev_newton->dy[i] * h_prev_newton->dy[i])
      + r3 * v->val[i] * h_prev_newton->dy[i]
      );
  }
  return result;
}

Ord CustomWeakFormRichardsRK::CustomResidualFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord>* *ext) const
{
  return Ord(10);
}

VectorFormVol<double>* CustomWeakFormRichardsRK::CustomResidualFormVol::clone() const
{
  return new CustomResidualFormVol(*this);
}