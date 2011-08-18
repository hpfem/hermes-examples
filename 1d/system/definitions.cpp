#include "definitions.h"

double CustomExactFunction1::val(double x) 
{
  return cos(M_PI*x/2);
}
  
double CustomExactFunction1::dx(double x) 
{
  return -sin(M_PI*x/2)*(M_PI/2.);
}
  
double CustomExactFunction1::ddxx(double x) 
{
  return -cos(M_PI*x/2)*(M_PI/2.)*(M_PI/2.);
}


double CustomExactFunction2::val(double x) 
{
  return 1. - (exp(K*x) + exp(-K*x))/(exp(K) + exp(-K));
}
  
double CustomExactFunction2::dx(double x) 
{
  return -K*(exp(K*x) - exp(-K*x))/(exp(K) + exp(-K));
}
  
double CustomExactFunction2::ddxx(double x) 
{
  return -K*K*(exp(K*x) + exp(-K*x))/(exp(K) + exp(-K));
}


CustomRightHandSide1::CustomRightHandSide1(double K, double d_u, double sigma)
  : Hermes::Hermes2DFunction<double>(), d_u(d_u), sigma(sigma) 
{
  cef1 = new CustomExactFunction1();
  cef2 = new CustomExactFunction2(K);
}

double CustomRightHandSide1::value(double x, double y) const 
{
  double Laplace_u = cef1->ddxx(x);
  double u = cef1->val(x);
  double v = cef2->val(x);
  return -d_u * d_u * Laplace_u - u + sigma * v;
}

Hermes::Ord CustomRightHandSide1::ord(Hermes::Ord x, Hermes::Ord y) const 
{
  return Hermes::Ord(10);
}

CustomRightHandSide1::~CustomRightHandSide1() 
{ 
  delete cef1; 
  delete cef2;
}

CustomRightHandSide2::CustomRightHandSide2(double K, double d_v)
      : Hermes::Hermes2DFunction<double>(), d_v(d_v) 
{
  cef1 = new CustomExactFunction1();
  cef2 = new CustomExactFunction2(K);
}

double CustomRightHandSide2::value(double x, double y) const 
{
  double Laplace_v = cef2->ddxx(x);
  double u = cef1->val(x);
  double v = cef2->val(x);
  return -d_v*d_v * Laplace_v - u + v;
}

Hermes::Ord CustomRightHandSide2::ord(Hermes::Ord x, Hermes::Ord y) const 
{
  return Hermes::Ord(10);
}

CustomRightHandSide2::~CustomRightHandSide2() 
{ 
  delete cef1; 
  delete cef2;
}


ExactSolutionFitzHughNagumo1::ExactSolutionFitzHughNagumo1(Mesh* mesh)
     : ExactSolutionScalar<double>(mesh) 
{
  cef1 = new CustomExactFunction1();
}

double ExactSolutionFitzHughNagumo1::value (double x, double y) const 
{
  return cef1->val(x);
}

void ExactSolutionFitzHughNagumo1::derivatives (double x, double y, double& dx, double& dy) const 
{
  dx = cef1->dx(x);
  dy = 0;
}

Hermes::Ord ExactSolutionFitzHughNagumo1::ord(Hermes::Ord x, Hermes::Ord y) const 
{
  return Hermes::Ord(10);
}

ExactSolutionFitzHughNagumo1::~ExactSolutionFitzHughNagumo1() 
{
  delete cef1;
}

ExactSolutionFitzHughNagumo2::ExactSolutionFitzHughNagumo2(Mesh* mesh, double K)
     : ExactSolutionScalar<double>(mesh) 
{
  cef2 = new CustomExactFunction2(K);
}

double ExactSolutionFitzHughNagumo2::value (double x, double y) const 
{
  return cef2->val(x);
}

void ExactSolutionFitzHughNagumo2::derivatives (double x, double y, double& dx, double& dy) const 
{
  dx = cef2->dx(x);
  dy = 0;
}

Hermes::Ord ExactSolutionFitzHughNagumo2::ord(Hermes::Ord x, Hermes::Ord y) const 
{
  return Hermes::Ord(10);
}

ExactSolutionFitzHughNagumo2::~ExactSolutionFitzHughNagumo2() 
{
  delete cef2;
}

double CustomResidual1::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                              Geom<double> *e, ExtData<double> *ext) const
{
   double result = 0;
   for (int i = 0; i < n; i++) 
   {
     result += wt[i] * (    d_u*d_u * u_ext[0]->dx[i]*v->dx[i] 
                          - u_ext[0]->val[i]*v->val[i] 
                          + sigma*u_ext[1]->val[i]*v->val[i]
                          - g1->value(e->x[i], 0)*v->val[i]
                       );
   }
 
   return result;
}

Hermes::Ord CustomResidual1::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
                         Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const 
{
   Hermes::Ord result = Hermes::Ord(0);
   for (int i = 0; i < n; i++) 
   {
     result += wt[i] * (    d_u*d_u * u_ext[0]->dx[i]*v->dx[i] 
                          - u_ext[0]->val[i]*v->val[i] 
                          + sigma*u_ext[1]->val[i]*v->val[i]
                          - g1->ord(e->x[i], Hermes::Ord(0))*v->val[i]
                        );
   }

   return result;
}

VectorFormVol<double>* CustomResidual1::clone() 
{
  return new CustomResidual1(*this);
}

double CustomResidual2::value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                              Geom<double> *e, ExtData<double> *ext) const
{
   double result = 0;
   for (int i = 0; i < n; i++) 
   {
     result += wt[i] * (    d_v*d_v * u_ext[1]->dx[i]*v->dx[i] 
                          - u_ext[0]->val[i]*v->val[i] 
                          + u_ext[1]->val[i]*v->val[i]
                          - g2->value(e->x[i], 0)*v->val[i]
                       );
   }
 
   return result;
  }

Hermes::Ord CustomResidual2::ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
                                 Geom<Hermes::Ord> *e, ExtData<Hermes::Ord> *ext) const 
{
   Hermes::Ord result = Hermes::Ord(0);
   for (int i = 0; i < n; i++) 
   {
     result += wt[i] * (    d_v*d_v * u_ext[1]->dx[i]*v->dx[i] 
                          - u_ext[0]->val[i]*v->val[i] 
                          + u_ext[1]->val[i]*v->val[i]
                          - g2->ord(e->x[i], Hermes::Ord(0))*v->val[i]
                        );
   }

  return result;
} 

VectorFormVol<double>* CustomResidual2::clone() 
{
  return new CustomResidual2(*this);
}

CustomWeakForm::CustomWeakForm(CustomRightHandSide1* g1, CustomRightHandSide2* g2) : WeakForm<double>(2) 
{
  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, Hermes::HERMES_ANY, new Hermes::Hermes1DFunction<double>(g1->d_u * g1->d_u)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(0, 0, Hermes::HERMES_ANY, new Hermes::Hermes2DFunction<double>(-1.0)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(0, 1, Hermes::HERMES_ANY, new Hermes::Hermes2DFunction<double>(g1->sigma), HERMES_NONSYM));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(1, 0, Hermes::HERMES_ANY, new Hermes::Hermes2DFunction<double>(-1.0), HERMES_NONSYM));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(1, 1, Hermes::HERMES_ANY, new Hermes::Hermes1DFunction<double>(g2->d_v * g2->d_v)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(1, 1, Hermes::HERMES_ANY, new Hermes::Hermes2DFunction<double>(1.0)));

  // Residual.
  add_vector_form(new CustomResidual1(g1->d_u, g1->sigma, g1));
  add_vector_form(new CustomResidual2(g2->d_v, g2));
}
