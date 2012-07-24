#include "definitions.h"

double CustomInitialConditionWave::value (double x, double y) const 
{
  return exp(-x*x - y*y);
}

void CustomInitialConditionWave::derivatives (double x, double y, double& dx, double& dy) const 
{
  dx = exp(-x*x - y*y) * (-2*x);
  dy = exp(-x*x - y*y) * (-2*y);
}

Ord CustomInitialConditionWave::ord(Ord x, Ord y) const 
{
  return Ord(10);
}
  
MeshFunction<double>* CustomInitialConditionWave::clone()
{
  return new CustomInitialConditionWave(this->mesh);
}

CustomWeakFormWave::CustomWeakFormWave(double tau, double c_squared, Solution<double>* u_prev_sln,
                                       Solution<double>* v_prev_sln) : WeakForm<double>(2) 
{
  // Volumetric matrix forms.
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(0, 1, HERMES_ANY, new Hermes2DFunction<double>(1.0), HERMES_NONSYM));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(1, 0, HERMES_ANY, new Hermes1DFunction<double>(-c_squared), HERMES_NONSYM));

  // Volumetric surface forms.
  add_vector_form(new VectorFormVolWave_0());
  add_vector_form(new VectorFormVolWave_1(c_squared));
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::VectorFormVolWave_0::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                            Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = Scalar(0);

  for (int i = 0; i < n; i++)
    result += wt[i] * u_ext[1]->val[i] * v->val[i];

  return result;
}

double CustomWeakFormWave::VectorFormVolWave_0::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                      Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormWave::VectorFormVolWave_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
                                                 ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

VectorFormVol<double>* CustomWeakFormWave::VectorFormVolWave_0::clone() 
{
  return new VectorFormVolWave_0(*this);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormWave::VectorFormVolWave_1::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                                            Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = Scalar(0);

  for (int i = 0; i < n; i++)
    result += wt[i] * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]);

  return - c_squared * result;
}

double CustomWeakFormWave::VectorFormVolWave_1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                                                      Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormWave::VectorFormVolWave_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
                                                 ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

VectorFormVol<double>* CustomWeakFormWave::VectorFormVolWave_1::clone() 
{
  return new VectorFormVolWave_1(*this);
}

