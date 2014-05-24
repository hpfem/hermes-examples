#include "definitions.h"

CustomWeakFormLinearElasticity::CustomWeakFormLinearElasticity(double E, double nu, Hermes2DFunction<double>* Eo11d, Hermes2DFunction<double>* Eo12d, Hermes2DFunction<double>* Eo22d, Hermes2DFunction<double>* Eo33d, double Eo11, double Eo12, double Eo22, double Eo33, double rho_g,
                                 std::string surface_force_bdy, double f0, double f1) : WeakForm<double>(3)
{
  double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
  double mu = E / (2 * (1 + nu));
  double kappa = 2.0 * (2.0 * E - 3.0 * nu) / 3 / (1 + nu) / (1 - 2 * nu);

  // SINGLE-COMPONENT FORMS. USEFUL FOR MULTIMESH, DO NOT REMOVE.
  // Jacobian.
  add_matrix_form(new CustomJacobianElast00(0, 0, lambda, mu, kappa, Eo11, Eo12, Eo22, Eo33));
  add_matrix_form(new CustomJacobianElast01(0, 1, lambda, mu, kappa, Eo11, Eo12, Eo22, Eo33));
  add_matrix_form(new CustomJacobianElast02(0, 2, lambda, mu, kappa, Eo11, Eo12, Eo22, Eo33));
  add_matrix_form(new CustomJacobianElast10(1, 0, lambda, mu, kappa, Eo11, Eo12, Eo22, Eo33));
  add_matrix_form(new CustomJacobianElast11(1, 1, lambda, mu, kappa, Eo11, Eo12, Eo22, Eo33));
  add_matrix_form(new CustomJacobianElast12(1, 2, lambda, mu, kappa, Eo11, Eo12, Eo22, Eo33));
  add_matrix_form(new CustomJacobianElast20(2, 0, lambda, mu, kappa, Eo11, Eo12, Eo22, Eo33));
  add_matrix_form(new CustomJacobianElast21(2, 1, lambda, mu, kappa, Eo11, Eo12, Eo22, Eo33));
  add_matrix_form(new CustomJacobianElast22(2, 2, lambda, mu, kappa, Eo11, Eo12, Eo22, Eo33));

  //Residuals
  add_vector_form(new CustomVectorRes0(0, lambda, mu, kappa, Eo11d, Eo12d, Eo22d, Eo33d, Eo11, Eo12, Eo22, Eo33));
  add_vector_form(new CustomVectorRes1(1, lambda, mu, kappa, Eo11d, Eo12d, Eo22d, Eo33d, Eo11, Eo12, Eo22, Eo33));
  add_vector_form(new CustomVectorRes2(2, lambda, mu, kappa, Eo11d, Eo12d, Eo22d, Eo33d, Eo11, Eo12, Eo22, Eo33));
}

// Jacobian lin elast 0-0
double CustomWeakFormLinearElasticity::CustomJacobianElast00::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * ((lambda + 2.0 * mu) * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
    //result += wt[i] * ((kappa + 4.0 * mu / 3.0) * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
    result += - wt[i] * (4.0*mu/3.0 * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomJacobianElast00::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * ((lambda + 2.0 * mu) * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
	//result += wt[i] * ((kappa + 4.0 * mu / 3.0) * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
	result += - wt[i] * (4.0*mu/3.0 * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
  }
  return result;
}
MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast00::clone() const
{
        return new CustomJacobianElast00(*this);
}

// Jacobian lin elast 0-1
double CustomWeakFormLinearElasticity::CustomJacobianElast01::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += - wt[i] * (-2.0*mu/3.0 * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomJacobianElast01::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
    result += - wt[i] * (-2.0*mu/3.0 * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
  }
  return result;
}
MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast01::clone() const
{
        return new CustomJacobianElast01(*this);
}

// Jacobian lin elast 0-2
double CustomWeakFormLinearElasticity::CustomJacobianElast02::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += wt[i] * (u->val[i] * v->dx[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomJacobianElast02::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += wt[i] * (u->val[i] * v->dx[i]);
  }
  return result;
}
MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast02::clone() const
{
        return new CustomJacobianElast02(*this);
}

// Jacobian lin elast 1-0
double CustomWeakFormLinearElasticity::CustomJacobianElast10::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * (mu * u->dy[i] * v->dx[i] + lambda * u->dx[i] * v->dy[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dx[i] * v->dy[i] + mu * u->dy[i] * v->dx[i]);
    result += - wt[i] * ((-2.0*mu/3.0) * u->dx[i] * v->dy[i] + mu * u->dy[i] * v->dx[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomJacobianElast10::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * (mu * u->dy[i] * v->dx[i] + lambda * u->dx[i] * v->dy[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dx[i] * v->dy[i] + mu * u->dy[i] * v->dx[i]);
	result += - wt[i] * ((-2.0*mu/3.0) * u->dx[i] * v->dy[i] + mu * u->dy[i] * v->dx[i]);
  }
  return result;
}
MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast10::clone() const
{
        return new CustomJacobianElast10(*this);
}

// Jacobian lin elast 1-1
double CustomWeakFormLinearElasticity::CustomJacobianElast11::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * ((lambda + 2 * mu) * u->dy[i] * v->dy[i] + mu * u->dx[i] * v->dx[i]);
	//result += wt[i] * ((kappa + 4.0 * mu / 3.0) * u->dy[i] * v->dy[i] + mu * u->dx[i] * v->dx[i]);
	result += - wt[i] * (mu * u->dx[i] * v->dx[i] + 4.0*mu/3.0 * u->dy[i] * v->dy[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomJacobianElast11::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * ((lambda + 2 * mu) * u->dy[i] * v->dy[i] + mu * u->dx[i] * v->dx[i]);
	//result += wt[i] * ((kappa + 4.0 * mu / 3.0) * u->dy[i] * v->dy[i] + mu * u->dx[i] * v->dx[i]);
	result += - wt[i] * (mu * u->dx[i] * v->dx[i] + 4.0*mu/3.0 * u->dy[i] * v->dy[i]);
  }
  return result;
}
MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast11::clone() const
{
        return new CustomJacobianElast11(*this);
}

// Jacobian lin elast 1-2
double CustomWeakFormLinearElasticity::CustomJacobianElast12::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += wt[i] * (u->val[i] * v->dy[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomJacobianElast12::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += wt[i] * (u->val[i] * v->dy[i]);
  }
  return result;
}
MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast12::clone() const
{
        return new CustomJacobianElast12(*this);
}

// Jacobian lin elast 2-0
double CustomWeakFormLinearElasticity::CustomJacobianElast20::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += wt[i] * (u->dx[i] * v->val[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomJacobianElast20::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += wt[i] * (u->dx[i] * v->val[i]);
  }
  return result;
}
MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast20::clone() const
{
        return new CustomJacobianElast20(*this);
}

// Jacobian lin elast 2-1
double CustomWeakFormLinearElasticity::CustomJacobianElast21::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += wt[i] * (u->dy[i] * v->val[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomJacobianElast21::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += wt[i] * (u->dy[i] * v->val[i]);
  }
  return result;
}
MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast21::clone() const
{
        return new CustomJacobianElast21(*this);
}

// Jacobian lin elast 2-2
double CustomWeakFormLinearElasticity::CustomJacobianElast22::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *u, Func<double> *v,
                                                   Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += wt[i] * (1.0/kappa * u->val[i] * v->val[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomJacobianElast22::ord(int n, double *wt, Func<Ord> *u_ext[],
                                              Func<Ord> *u, Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * (mu * u->dx[i] * v->dy[i] + lambda * u->dy[i] * v->dx[i]);
    //result += wt[i] * ((kappa - 2.0 * mu / 3.0) * u->dy[i] * v->dx[i] + mu * u->dx[i] * v->dy[i]);
	result += wt[i] * (1.0/kappa * u->val[i] * v->val[i]);
  }
  return result;
}
MatrixFormVol<double>* CustomWeakFormLinearElasticity::CustomJacobianElast22::clone() const
{
        return new CustomJacobianElast22(*this);
}

//residuum lin elast 0
double CustomWeakFormLinearElasticity::CustomVectorRes0::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * ((2.0 * mu * (u_ext[0]->dx[i] - Eo11) + lambda * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11 - Eo22 - Eo33)) * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12) * v->dy[i]);
	//result += wt[i] * ((2.0 * mu * (u_ext[0]->dx[i] - Eo11) + (kappa - 2.0 * mu / 3.0) * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11 - Eo22 - Eo33)) * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12) * v->dy[i]);
	//result += wt[i] * ((2.0 * mu * (u_ext[0]->dx[i] - Eo11d->value(e->x[i],e->y[i])) + (kappa - 2.0 * mu / 3.0) * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11d->value(e->x[i],e->y[i]) - Eo22d->value(e->x[i],e->y[i]) - Eo33d->value(e->x[i],e->y[i]))) * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12d->value(e->x[i],e->y[i])) * v->dy[i]);
	result += - wt[i] * ((2.0 * mu * (u_ext[0]->dx[i] - Eo11d->value(e->x[i],e->y[i])) - (2.0*mu/3.0 * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11d->value(e->x[i],e->y[i]) - Eo22d->value(e->x[i],e->y[i]) - Eo33d->value(e->x[i],e->y[i])) + u_ext[2]->val[i])) * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12d->value(e->x[i],e->y[i])) * v->dy[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomVectorRes0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * ((2.0 * mu * (u_ext[0]->dx[i] - Eo11) + lambda * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11 - Eo22 - Eo33)) * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12) * v->dy[i]);
	//result += wt[i] * ((2.0 * mu * (u_ext[0]->dx[i] - Eo11) + (kappa - 2.0 * mu / 3.0) * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11 - Eo22 - Eo33)) * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12) * v->dy[i]);
	//result += wt[i] * ((2.0 * mu * (u_ext[0]->dx[i] - Eo11d->value(e->x[i],e->y[i])) + (kappa - 2.0 * mu / 3.0) * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11d->value(e->x[i],e->y[i]) - Eo22d->value(e->x[i],e->y[i]) - Eo33d->value(e->x[i],e->y[i]))) * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12d->value(e->x[i],e->y[i])) * v->dy[i]);
	result += - wt[i] * ((2.0 * mu * (u_ext[0]->dx[i] - Eo11d->value(e->x[i],e->y[i])) - (2.0*mu/3.0 * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11d->value(e->x[i],e->y[i]) - Eo22d->value(e->x[i],e->y[i]) - Eo33d->value(e->x[i],e->y[i])) + u_ext[2]->val[i])) * v->dx[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12d->value(e->x[i],e->y[i])) * v->dy[i]);
  }
  return result;
}
VectorFormVol<double>* CustomWeakFormLinearElasticity::CustomVectorRes0::clone() const
{
        return new CustomVectorRes0(*this);
}

//residuum lin elast 1
double CustomWeakFormLinearElasticity::CustomVectorRes1::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
    //result += wt[i] * ((2.0 * mu * (u_ext[1]->dy[i] - Eo22) + lambda * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11 - Eo22 - Eo33)) * v->dy[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12) * v->dx[i]);
	//result += wt[i] * ((2.0 * mu * (u_ext[1]->dy[i] - Eo22) + (kappa - 2.0 * mu / 3.0) * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11 - Eo22 - Eo33)) * v->dy[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12) * v->dx[i]);
	//result += wt[i] * ((2.0 * mu * (u_ext[1]->dy[i] - Eo22d->value(e->x[i],e->y[i])) + (kappa - 2.0 * mu / 3.0) * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11d->value(e->x[i],e->y[i]) - Eo22d->value(e->x[i],e->y[i]) - Eo33d->value(e->x[i],e->y[i]))) * v->dy[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12d->value(e->x[i],e->y[i])) * v->dx[i]);
	result += - wt[i] * ((2.0 * mu * (u_ext[1]->dy[i] - Eo22d->value(e->x[i],e->y[i])) - (2.0*mu/3.0 * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11d->value(e->x[i],e->y[i]) - Eo22d->value(e->x[i],e->y[i]) - Eo33d->value(e->x[i],e->y[i])) + u_ext[2]->val[i])) * v->dy[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12d->value(e->x[i],e->y[i])) * v->dx[i]);
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomVectorRes1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
	//result += wt[i] * ((2.0 * mu * (u_ext[1]->dy[i] - Eo22) + lambda * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11 - Eo22 - Eo33)) * v->dy[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12) * v->dx[i]);
    //result += wt[i] * ((2.0 * mu * (u_ext[1]->dy[i] - Eo22) + (kappa - 2.0 * mu / 3.0) * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11 - Eo22 - Eo33)) * v->dy[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12) * v->dx[i]);
	//result += wt[i] * ((2.0 * mu * (u_ext[1]->dy[i] - Eo22d->value(e->x[i],e->y[i])) + (kappa - 2.0 * mu / 3.0) * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11d->value(e->x[i],e->y[i]) - Eo22d->value(e->x[i],e->y[i]) - Eo33d->value(e->x[i],e->y[i]))) * v->dy[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12d->value(e->x[i],e->y[i])) * v->dx[i]);
	result += - wt[i] * ((2.0 * mu * (u_ext[1]->dy[i] - Eo22d->value(e->x[i],e->y[i])) - (2.0*mu/3.0 * (u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11d->value(e->x[i],e->y[i]) - Eo22d->value(e->x[i],e->y[i]) - Eo33d->value(e->x[i],e->y[i])) + u_ext[2]->val[i])) * v->dy[i] + mu * (u_ext[0]->dy[i] + u_ext[1]->dx[i] - 2.0 * Eo12d->value(e->x[i],e->y[i])) * v->dx[i]);
  }
  return result;
}
VectorFormVol<double>* CustomWeakFormLinearElasticity::CustomVectorRes1::clone() const
{
        return new CustomVectorRes1(*this);
}

//residuum lin elast 2
double CustomWeakFormLinearElasticity::CustomVectorRes2::value(int n, double *wt, Func<double> *u_ext[],
                                                   Func<double> *v, Geom<double> *e, Func<double>* *ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++)
  {
	result += wt[i] * (u_ext[2]->val[i]/kappa + u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11d->value(e->x[i],e->y[i]) - Eo22d->value(e->x[i],e->y[i]) - Eo33d->value(e->x[i],e->y[i])) * v->val[i];
  }
  return result;
}
Ord CustomWeakFormLinearElasticity::CustomVectorRes2::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                              Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * (u_ext[2]->val[i]/kappa + u_ext[0]->dx[i] + u_ext[1]->dy[i] - Eo11d->value(e->x[i],e->y[i]) - Eo22d->value(e->x[i],e->y[i]) - Eo33d->value(e->x[i],e->y[i])) * v->val[i];
  }
  return result;
}
VectorFormVol<double>* CustomWeakFormLinearElasticity::CustomVectorRes2::clone() const
{
        return new CustomVectorRes2(*this);
}

// Custom exact function Eo11
double CustomExactFunctionEo11::val(double x, double y)
{
  double stepBot = dist * (1.0 / 3.0 + 0.0 * Hermes::sin(8.0 * M_PI * x) / 50 < y);
  double stepTop = dist * (2.0 / 3.0 + 0.0 * Hermes::sin(2.0 * M_PI * x) / 25 < y);
  return stepBot - stepTop;
}
double CustomExactFunctionEo11::dx(double x, double y)
{
  return 0.0; // this does not correspond to the value
}
double CustomExactFunctionEo11::ddxx(double x, double y)
{
  return 0.0; // this also probably
}

// Custom exact function Eo12
double CustomExactFunctionEo12::val(double x, double y)
{
  double stepBot = dist * (1.0 / 3.0 + 0.0 * Hermes::sin(8.0 * M_PI * x) / 50 < y);
  double stepTop = dist * (2.0 / 3.0 + 0.0 * Hermes::sin(2.0 * M_PI * x) / 25 < y);
  return stepBot - stepTop;
}
double CustomExactFunctionEo12::dx(double x, double y)
{
  return 0.0; // this does not correspond to the value
}

// Custom exact function Eo22
double CustomExactFunctionEo22::val(double x, double y)
{
  double stepBot = dist * (1.0 / 3.0 + 0.0 * Hermes::sin(8.0 * M_PI * x) / 50 < y);
  double stepTop = dist * (2.0 / 3.0 + 0.0 * Hermes::sin(2.0 * M_PI * x) / 25 < y);
  return stepBot - stepTop;
}
double CustomExactFunctionEo22::dx(double x, double y)
{
  return 0.0; // this does not correspond to the value
}

// Custom exact function Eo33
double CustomExactFunctionEo33::val(double x, double y)
{
  double stepBot = dist * (1.0 / 3.0 + 0.0 * Hermes::sin(8.0 * M_PI * x) / 50 < y);
  double stepTop = dist * (2.0 / 3.0 + 0.0 * Hermes::sin(2.0 * M_PI * x) / 25 < y);
  return stepBot - stepTop;
}
double CustomExactFunctionEo33::dx(double x, double y)
{
  return 0.0; // this does not correspond to the value
}

// Eo11 Hermes2DFunction
CustomEo11::CustomEo11(double dist)
  : Hermes2DFunction<double>(), dist(dist)
{
  cefEo11 = new CustomExactFunctionEo11(dist);
}
double CustomEo11::value(double x, double y) const
{
  return cefEo11->val(x,y);
}
Ord CustomEo11::value(Ord x, Ord y) const
{
  return Ord(10);
}
CustomEo11::~CustomEo11()
{
  delete cefEo11;
}

// Eo12 Hermes2DFunction
CustomEo12::CustomEo12(double dist)
  : Hermes2DFunction<double>(), dist(dist)
{
  cefEo12 = new CustomExactFunctionEo12(dist);
}
double CustomEo12::value(double x, double y) const
{
  return cefEo12->val(x,y);
}
Ord CustomEo12::value(Ord x, Ord y) const
{
  return Ord(10);
}
CustomEo12::~CustomEo12()
{
  delete cefEo12;
}

// Eo22 Hermes2DFunction
CustomEo22::CustomEo22(double dist)
  : Hermes2DFunction<double>(), dist(dist)
{
  cefEo22 = new CustomExactFunctionEo22(dist);
}
double CustomEo22::value(double x, double y) const
{
  return cefEo22->val(x,y);
}
Ord CustomEo22::value(Ord x, Ord y) const
{
  return Ord(10);
}
CustomEo22::~CustomEo22()
{
  delete cefEo22;
}

// Eo33 Hermes2DFunction
CustomEo33::CustomEo33(double dist)
  : Hermes2DFunction<double>(), dist(dist)
{
  cefEo33 = new CustomExactFunctionEo33(dist);
}
double CustomEo33::value(double x, double y) const
{
  return cefEo33->val(x,y);
}
Ord CustomEo33::value(Ord x, Ord y) const
{
  return Ord(10);
}
CustomEo33::~CustomEo33()
{
  delete cefEo33;
}

// Exact Function Solution Eo11
ExactSolutionEo11::ExactSolutionEo11(MeshSharedPtr mesh, double dist)
     : ExactSolutionScalar<double>(mesh), dist(dist)
{
  cefEo11 = new CustomExactFunctionEo11(dist);
}
double ExactSolutionEo11::value(double x, double y) const
{
  return cefEo11->val(x,y);
}
void ExactSolutionEo11::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = cefEo11->dx(x,y);
  dy = -cefEo11->dx(x,y);
}
Ord ExactSolutionEo11::ord(double x, double y) const
{
  return Ord(10);
}
ExactSolutionEo11::~ExactSolutionEo11()
{
  delete cefEo11;
}
MeshFunction<double>* ExactSolutionEo11::clone() const
{
  return new ExactSolutionEo11(this->mesh,this->dist);
}

// Exact Function Solution Eo12
ExactSolutionEo12::ExactSolutionEo12(MeshSharedPtr mesh, double dist)
     : ExactSolutionScalar<double>(mesh), dist(dist)
{
  cefEo12 = new CustomExactFunctionEo12(dist);
}
double ExactSolutionEo12::value(double x, double y) const
{
  return cefEo12->val(x,y);
}
void ExactSolutionEo12::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = cefEo12->dx(x,y);
  dy = -cefEo12->dx(x,y);
}
Ord ExactSolutionEo12::ord(double x, double y) const
{
  return Ord(10);
}
ExactSolutionEo12::~ExactSolutionEo12()
{
  delete cefEo12;
}
MeshFunction<double>* ExactSolutionEo12::clone() const
{
  return new ExactSolutionEo12(this->mesh,this->dist);
}

// Exact Function Solution Eo22
ExactSolutionEo22::ExactSolutionEo22(MeshSharedPtr mesh, double dist)
     : ExactSolutionScalar<double>(mesh), dist(dist)
{
  cefEo22 = new CustomExactFunctionEo22(dist);
}
double ExactSolutionEo22::value(double x, double y) const
{
  return cefEo22->val(x,y);
}
void ExactSolutionEo22::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = cefEo22->dx(x,y);
  dy = -cefEo22->dx(x,y);
}
Ord ExactSolutionEo22::ord(double x, double y) const
{
  return Ord(10);
}
ExactSolutionEo22::~ExactSolutionEo22()
{
  delete cefEo22;
}
MeshFunction<double>* ExactSolutionEo22::clone() const
{
  return new ExactSolutionEo22(this->mesh,this->dist);
}

// Exact Function Solution Eo33
ExactSolutionEo33::ExactSolutionEo33(MeshSharedPtr mesh, double dist)
     : ExactSolutionScalar<double>(mesh), dist(dist)
{
  cefEo33 = new CustomExactFunctionEo33(dist);
}
double ExactSolutionEo33::value(double x, double y) const
{
  return cefEo33->val(x,y);
}
void ExactSolutionEo33::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = cefEo33->dx(x,y);
  dy = -cefEo33->dx(x,y);
}
Ord ExactSolutionEo33::ord(double x, double y) const
{
  return Ord(10);
}
ExactSolutionEo33::~ExactSolutionEo33()
{
  delete cefEo33;
}
MeshFunction<double>* ExactSolutionEo33::clone() const
{
  return new ExactSolutionEo33(this->mesh,this->dist);
}

// Exact Function Solution tr(Eo)
ExactSolutionTrEo::ExactSolutionTrEo(MeshSharedPtr mesh, double Eo11, double Eo22, double Eo33)
     : ExactSolutionScalar<double>(mesh), Eo11(Eo11), Eo22(Eo22), Eo33(Eo33)
{
  cefEo11 = new CustomExactFunctionEo11(Eo11);
  cefEo22 = new CustomExactFunctionEo22(Eo22);
  cefEo33 = new CustomExactFunctionEo33(Eo33);
}
double ExactSolutionTrEo::value(double x, double y) const
{
  return cefEo11->val(x,y)+cefEo22->val(x,y)+cefEo33->val(x,y);
}
void ExactSolutionTrEo::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = cefEo11->dx(x,y);
  dy = -cefEo11->dx(x,y);
}
Ord ExactSolutionTrEo::ord(double x, double y) const
{
  return Ord(10);
}
ExactSolutionTrEo::~ExactSolutionTrEo()
{
  delete cefEo11, cefEo22, cefEo33;
}
MeshFunction<double>* ExactSolutionTrEo::clone() const
{
  return new ExactSolutionTrEo(this->mesh,this->Eo11,this->Eo22,this->Eo33);
}

// Custom filter S11
void CustomFilterS11::filter_fn(int n, double* x, double* y, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                      double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = 2.0 * mu * (dx.at(0)[i] - this->Eo11d->value(x[i],y[i])) - 2.0/3.0*mu * (dx.at(0)[i] + dy.at(1)[i]-this->Eo11d->value(x[i],y[i])-this->Eo22d->value(x[i],y[i])-this->Eo33d->value(x[i],y[i])) - values.at(2)[i];
    outdx[i] = 0.0;
    outdy[i] = 0.0;
  }
}

// Custom filter S12
void CustomFilterS12::filter_fn(int n, double* x, double* y, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                      double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = mu * ((dy.at(0)[i] + dx.at(1)[i]) - 2.0 * this->Eo12d->value(x[i],y[i]));
    outdx[i] = 0.0;
    outdy[i] = 0.0;
  }
}

// Custom filter S22
void CustomFilterS22::filter_fn(int n, double* x, double* y, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                      double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = 2.0 * mu * (dy.at(1)[i] - this->Eo22d->value(x[i],y[i])) - 2.0/3.0*mu * (dx.at(0)[i] + dy.at(1)[i]-this->Eo11d->value(x[i],y[i])-this->Eo22d->value(x[i],y[i])-this->Eo33d->value(x[i],y[i])) - values.at(2)[i];
    outdx[i] = 0.0;
    outdy[i] = 0.0;
  }
}

// Custom filter S33
void CustomFilterS33::filter_fn(int n, double* x, double* y, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                      double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = - 2.0 * mu * this->Eo33d->value(x[i],y[i]) - 2.0/3.0*mu * (dx.at(0)[i] + dy.at(1)[i]-this->Eo11d->value(x[i],y[i])-this->Eo22d->value(x[i],y[i])-this->Eo33d->value(x[i],y[i])) - values.at(2)[i];
    outdx[i] = 0.0;
    outdy[i] = 0.0;
  }
}

// Custom filter von Mises (with distortions)
void CustomFilter_vM::filter_fn(int n, double* x, double* y, Hermes::vector<double*> values, Hermes::vector<double*> dx, Hermes::vector<double*> dy,
                      double* out, double* outdx, double* outdy)
{
  for (int i = 0; i < n; i++)
  {
    double S11S22 = 2.0 * mu * (dx.at(0)[i] - dy.at(1)[i] - (this->Eo11d->value(x[i],y[i]) - this->Eo22d->value(x[i],y[i])));
    double S11S33 = 2.0 * mu * (dx.at(0)[i] - (this->Eo11d->value(x[i],y[i]) - this->Eo33d->value(x[i],y[i])));
    double S22S33 = 2.0 * mu * (dy.at(1)[i] - (this->Eo22d->value(x[i],y[i]) - this->Eo33d->value(x[i],y[i])));
    double S12 = mu * (dy.at(0)[i] + dx.at(1)[i] - this->Eo12d->value(x[i],y[i]));
    out[i] = 1.0/Hermes::pow(2.0, 0.5) * Hermes::pow((Hermes::pow(S11S22, 2.0) + Hermes::pow(S11S33, 2.0) + Hermes::pow(S22S33, 2.0) + 6 * Hermes::pow(S12, 2.0)), 0.5);
    outdx[i] = 0.0;
    outdy[i] = 0.0;
  }
}
