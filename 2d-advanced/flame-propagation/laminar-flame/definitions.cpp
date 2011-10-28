#include "definitions.h"

CustomWeakForm::CustomWeakForm(double Le, double alpha, double beta, double kappa, double x1, Hermes::Hermes2D::Filter<double>* omega, Hermes::Hermes2D::Filter<double>* omega_dt, Hermes::Hermes2D::Filter<double>* omega_dc) : WeakForm<double>(2), Le(Le), alpha(alpha), beta(beta), kappa(kappa), x1(x1)
{
  // Stationary Jacobian.
  MatrixFormVol<double>* mfv = new JacobianFormVol_0_0();
  mfv->ext.push_back(omega_dt);
  add_matrix_form(mfv);
  MatrixFormSurf<double>* mfs = new JacobianFormSurf_0_0("Neumann", kappa);
  add_matrix_form_surf(mfs);
  mfv = new JacobianFormVol_0_1();
  mfv->ext.push_back(omega_dc);
  add_matrix_form(mfv);
  mfv = new JacobianFormVol_1_0();
  mfv->ext.push_back(omega_dt);
  add_matrix_form(mfv);
  mfv = new JacobianFormVol_1_1(Le);
  mfv->ext.push_back(omega_dc);
  add_matrix_form(mfv);
  
  // Stationary residual.
  VectorFormVol<double>* vfv = new ResidualFormVol_0();
  vfv->ext.push_back(omega);
  add_vector_form(vfv);
  VectorFormSurf<double>* vfs = new ResidualFormSurf_0("Neumann", kappa);
  add_vector_form_surf(vfs);
  vfv = new ResidualFormVol_1(Le);
  vfv->ext.push_back(omega);
  add_vector_form(vfv);
}

double CustomWeakForm::JacobianFormVol_0_0::value(int n, double *wt, Func<double> *u_ext[], 
                               Func<double> *vj, Func<double> *vi, Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  Func<double>* domegadt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (  vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]
                      - domegadt->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

Ord CustomWeakForm::JacobianFormVol_0_0::ord(int n, double *wt, Func<Ord> *u_ext[], 
                            Func<Ord> *vj, Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  Ord result = Ord(0);
  Func<Ord>* domegadt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (  vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]
                      - domegadt->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

MatrixFormVol<double>* CustomWeakForm::JacobianFormVol_0_0::clone() 
{
  return new JacobianFormVol_0_0(*this);
}

double CustomWeakForm::JacobianFormSurf_0_0::value(int n, double *wt, Func<double> *u_ext[], 
                               Func<double> *vj, Func<double> *vi, Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * vj->val[i] * vi->val[i]);
  return result;
}

Ord CustomWeakForm::JacobianFormSurf_0_0::ord(int n, double *wt, Func<Ord> *u_ext[], 
                            Func<Ord> *vj, Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * vj->val[i] * vi->val[i]);
  return result;
}

MatrixFormSurf<double>* CustomWeakForm::JacobianFormSurf_0_0::clone() 
{
  return new JacobianFormSurf_0_0(*this);
}

double CustomWeakForm::JacobianFormVol_0_1::value(int n, double *wt, Func<double> *u_ext[], 
                               Func<double> *vj, Func<double> *vi, Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  Func<double>* domegady = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (- domegady->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

Ord CustomWeakForm::JacobianFormVol_0_1::ord(int n, double *wt, Func<Ord> *u_ext[], 
                            Func<Ord> *vj, Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  Ord result = Ord(0);
  Func<Ord>* domegady = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (- domegady->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

MatrixFormVol<double>* CustomWeakForm::JacobianFormVol_0_1::clone() 
{
  return new JacobianFormVol_0_1(*this);
}

double CustomWeakForm::JacobianFormVol_1_0::value(int n, double *wt, Func<double> *u_ext[], 
                               Func<double> *vj, Func<double> *vi, Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  Func<double>* domegadt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( domegadt->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

Ord CustomWeakForm::JacobianFormVol_1_0::ord(int n, double *wt, Func<Ord> *u_ext[], 
                            Func<Ord> *vj, Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  Ord result = Ord(0);
  Func<Ord>* domegadt = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( domegadt->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

MatrixFormVol<double>* CustomWeakForm::JacobianFormVol_1_0::clone() 
{
  return new JacobianFormVol_1_0(*this);
}

double CustomWeakForm::JacobianFormVol_1_1::value(int n, double *wt, Func<double> *u_ext[], 
                               Func<double> *vj, Func<double> *vi, Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  Func<double>* domegady = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]) / Le
                      + domegady->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

Ord CustomWeakForm::JacobianFormVol_1_1::ord(int n, double *wt, Func<Ord> *u_ext[], 
                            Func<Ord> *vj, Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  Ord result = Ord(0);
  Func<Ord>* domegady = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i] / Le
                      + domegady->val[i] * vj->val[i] * vi->val[i] );
  return result;
}

MatrixFormVol<double>* CustomWeakForm::JacobianFormVol_1_1::clone() 
{
  return new JacobianFormVol_1_1(*this);
}

double CustomWeakForm::ResidualFormVol_0::value(int n, double *wt, Func<double> *u_ext[], 
                               Func<double> *vi, Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  Func<double>* titer = u_ext[0];
  Func<double>* omega = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (titer->dx[i] * vi->dx[i] + titer->dy[i] * vi->dy[i]) -
                        omega->val[i] * vi->val[i]);
  return result;
}

Ord CustomWeakForm::ResidualFormVol_0::ord(int n, double *wt, Func<Ord> *u_ext[], 
                            Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  Ord result = Ord(0);
  Func<Ord>* titer = u_ext[0];
  Func<Ord>* omega = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (titer->dx[i] * vi->dx[i] + titer->dy[i] * vi->dy[i]) -
                        omega->val[i] * vi->val[i]);
  return result;
}

VectorFormVol<double>* CustomWeakForm::ResidualFormVol_0::clone() 
{
  return new ResidualFormVol_0(*this);
}

double CustomWeakForm::ResidualFormSurf_0::value(int n, double *wt, Func<double> *u_ext[], 
                               Func<double> *vi, Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  Func<double>* t_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * t_prev_newton->val[i] * vi->val[i]);
  return result;
}

Ord CustomWeakForm::ResidualFormSurf_0::ord(int n, double *wt, Func<Ord> *u_ext[], 
                            Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  Ord result = Ord(0);
  Func<Ord>* t_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (kappa * t_prev_newton->val[i] * vi->val[i]);
  return result;
}

VectorFormSurf<double>* CustomWeakForm::ResidualFormSurf_0::clone() 
{
  return new ResidualFormSurf_0(*this);
}

double CustomWeakForm::ResidualFormVol_1::value(int n, double *wt, Func<double> *u_ext[], 
                               Func<double> *vi, Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  Func<double>* c_prev_newton = u_ext[1];
  Func<double>* omega = ext->fn[0];
  for (int i = 0; i < n; i++)
		result += wt[i] * ( (c_prev_newton->dx[i] * vi->dx[i] + c_prev_newton->dy[i] * vi->dy[i]) / Le +
                        omega->val[i] * vi->val[i]);
  return result;
}

Ord CustomWeakForm::ResidualFormVol_1::ord(int n, double *wt, Func<Ord> *u_ext[], 
                            Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  Ord result = Ord(0);
  Func<Ord>* c_prev_newton = u_ext[1];
  Func<Ord>* omega = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ( (c_prev_newton->dx[i] * vi->dx[i] + c_prev_newton->dy[i] * vi->dy[i]) / Le +
                        omega->val[i] * vi->val[i]);
  return result;
}

VectorFormVol<double>* CustomWeakForm::ResidualFormVol_1::clone() 
{
  return new ResidualFormVol_1(*this);
}

// Preconditioner weak forms.
double CustomWeakForm::PreconditionerForm_0::value(int n, double *wt, Func<double>* u_ext[], Func<double> *vj, 
                                     Func<double> *vi, Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]);
  return result;
}

Ord CustomWeakForm::PreconditionerForm_0::ord(int n, double *wt, Func<Ord>* u_ext[], Func<Ord> *vj, 
                                  Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  return vj->val[0] * vi->val[0] +  vj->dx[0] * vi->dx[0] + vj->dy[0] * vi->dy[0];
}

MatrixFormVol<double>* CustomWeakForm::PreconditionerForm_0::clone() 
{
  return new PreconditionerForm_0(*this);
}

double CustomWeakForm::PreconditionerForm_1::value(int n, double *wt, Func<double>* u_ext[], Func<double> *vj, 
                                     Func<double> *vi, Geom<double> *e, ExtData<double> *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((vj->dx[i] * vi->dx[i] + vj->dy[i] * vi->dy[i]) / Le );
  return result;
}

Ord CustomWeakForm::PreconditionerForm_1::ord(int n, double *wt, Func<Ord>* u_ext[], Func<Ord> *vj, 
                                  Func<Ord> *vi, Geom<Ord> *e, ExtData<Ord> *ext) const
{
  return vj->val[0] * vi->val[0] +  vj->dx[0] * vi->dx[0] + vj->dy[0] * vi->dy[0];
}

MatrixFormVol<double>* CustomWeakForm::PreconditionerForm_1::clone() 
{
  return new PreconditionerForm_1(*this);
}

void CustomFilter::filter_fn(int n, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                             double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = std::max(values.at(0)[i],0.0) - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    double t5 = (beta / (t3 * t3)) * values.at(1)[i];
    out[i] = t4 * values.at(1)[i];
    outdx[i] = t4 * (dx.at(1)[i] + dx.at(0)[i] * t5);
    outdy[i] = t4 * (dy.at(1)[i] + dy.at(0)[i] * t5);
  }
}

void CustomFilterDt::filter_fn(int n, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                               double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = std::max(values.at(0)[i],0.0) - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    double t5 = (beta / (t3 * t3));
    out[i] = t4 * t5 * values.at(1)[i];
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}

void CustomFilterDc::filter_fn(int n, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                               double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double t1 = std::max(values.at(0)[i],0.0) - 1.0;
    double t2 = t1 * beta;
    double t3 = 1.0 + t1 * alpha;
    double t4 = sqr(beta) / (2.0*Le) * exp(t2 / t3);
    out[i] = t4;
    outdx[i] = 0.0;
    outdy[i] = 0.0; // not important
  }
}
