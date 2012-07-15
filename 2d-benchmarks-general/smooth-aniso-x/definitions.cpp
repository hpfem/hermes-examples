#include "definitions.h"

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = Hermes::cos(x);
  dy = 0;
}

double CustomExactSolution::value(double x, double y) const 
{
  return Hermes::sin(x);
}

Ord CustomExactSolution::ord(Ord x, Ord y) const 
{
  return Ord(7);
}


double CustomFunction::value(double x, double y) const 
{
  return -Hermes::sin(x);
}

Ord CustomFunction::value(Ord x, Ord y) const 
{
  return Ord(7);
}


CustomWeakFormPoisson::CustomWeakFormPoisson(std::string bdy_marker_right) : WeakForm<double>(1) 
{
  // Jacobian.
  add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0));
   
  // Residual.
  add_vector_form(new DefaultResidualDiffusion<double>(0));
  add_vector_form(new DefaultVectorFormVol<double>(0, new CustomFunction));
  add_vector_form_surf(new DefaultVectorFormSurf<double>(0, new Hermes::Hermes2DFunction<double>(1.0), bdy_marker_right));
}
