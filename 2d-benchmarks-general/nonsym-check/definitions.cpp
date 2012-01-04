#include "definitions.h"

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = Hermes::cos(x);
  dy = 0;
}

double CustomExactSolution::value(double x, double y) const 
{
  return sin(x);
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(20);
}

MeshFunction<double>* CustomExactSolution::clone()
{
  return new CustomExactSolution(*this);
}


double CustomJacobian::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                             Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * u->dx[i] * v->val[i];
  return result;
}

Ord CustomJacobian::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(20);
}

MatrixFormVol<double>* CustomJacobian::clone()
{
  return new CustomJacobian(*this);
}

double CustomResidual::value(int n, double *wt, Func<double> *u_ext[],
                             Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  double result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (u_ext[0]->dx[i] - Hermes::sin(e->x[i]) - Hermes::cos(e->x[i])) * v->val[i];
  return result;
}

Ord CustomResidual::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return Ord(20);
}

VectorFormVol<double>* CustomResidual::clone()
{
  return new CustomResidual(*this);
}

CustomWeakForm::CustomWeakForm(std::string marker_bdy_right) : WeakForm<double>(1) 
{
  // Jacobian.
  add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0));
  add_matrix_form(new CustomJacobian(0, 0));

  // Residual.
  add_vector_form(new DefaultResidualDiffusion<double>(0));
  add_vector_form(new CustomResidual(0));
  add_vector_form_surf(new DefaultVectorFormSurf<double>(0, marker_bdy_right, 
                                                         new Hermes::Hermes2DFunction<double>(1.0)));
}
