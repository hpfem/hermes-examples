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
  return x*(100. - x)/2.5 * y/100 - 1000. + H_OFFSET;
}

/* Custom weak forms */

CustomWeakFormRichardsIEPicard::CustomWeakFormRichardsIEPicard(double time_step, MeshFunctionSharedPtr<double>  h_time_prev, MeshFunctionSharedPtr<double>  h_iter_prev, ConstitutiveRelations* constitutive) : WeakForm<double>(1), constitutive(constitutive) 
{
  // Jacobian.
  CustomJacobian* matrix_form = new CustomJacobian(0, 0, time_step);
  matrix_form->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(h_time_prev, h_iter_prev));
  add_matrix_form(matrix_form);

  // Residual.
  CustomResidual* vector_form = new CustomResidual(0, time_step);
  vector_form->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(h_time_prev, h_iter_prev));
  add_vector_form(vector_form);
}

double CustomWeakFormRichardsIEPicard::CustomJacobian::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                                                             Func<double> *v, Geom<double> *e, Func<double>* *ext) const 
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext[0];
  Func<double>* h_prev_picard = ext[1];
  for (int i = 0; i < n; i++)
  {
    double h_prev_newton_i = h_prev_newton->val[i] - H_OFFSET;
    double h_prev_picard_i = h_prev_picard->val[i] - H_OFFSET;
    double h_prev_time_i = h_prev_time->val[i] - H_OFFSET;
    result += wt[i] * (   static_cast<CustomWeakFormRichardsIEPicard*>(wf)->constitutive->C(h_prev_picard_i) * u->val[i] * v->val[i]
                        + static_cast<CustomWeakFormRichardsIEPicard*>(wf)->constitutive->K(h_prev_picard_i) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) * time_step
                        - static_cast<CustomWeakFormRichardsIEPicard*>(wf)->constitutive->dKdh(h_prev_picard_i) * u->dy[i] * v->val[i] * time_step
                      );
  }
  return result;
}

Ord CustomWeakFormRichardsIEPicard::CustomJacobian::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                                                        Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return Ord(10);
}

MatrixFormVol<double>* CustomWeakFormRichardsIEPicard::CustomJacobian::clone() const 
{
  return new CustomJacobian(*this);
}

double CustomWeakFormRichardsIEPicard::CustomResidual::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
                                                             Func<double>* *ext) const 
{
  double result = 0;
  Func<double>* h_prev_newton = u_ext[0];
  Func<double>* h_prev_time = ext[0];
  Func<double>* h_prev_picard = ext[1];
  for (int i = 0; i < n; i++)
  {
    double h_prev_newton_i = h_prev_newton->val[i] - H_OFFSET;
    double h_prev_picard_i = h_prev_picard->val[i] - H_OFFSET;
    double h_prev_time_i = h_prev_time->val[i] - H_OFFSET;
    result += wt[i] * (   static_cast<CustomWeakFormRichardsIEPicard*>(wf)->constitutive->C(h_prev_picard_i) * (h_prev_newton_i - h_prev_time_i) * v->val[i]
                        + static_cast<CustomWeakFormRichardsIEPicard*>(wf)->constitutive->K(h_prev_picard_i) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i]) * time_step
                        - static_cast<CustomWeakFormRichardsIEPicard*>(wf)->constitutive->dKdh(h_prev_picard_i) * h_prev_newton->dy[i] * v->val[i] * time_step
                       );
  }
  return result;
}

Ord CustomWeakFormRichardsIEPicard::CustomResidual::ord(int n, double *wt, Func<Ord> *u_ext[], 
                                                               Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const 
{
  return Ord(10);
}

VectorFormVol<double>* CustomWeakFormRichardsIEPicard::CustomResidual::clone() const 
{
  return new CustomResidual(*this);
}


