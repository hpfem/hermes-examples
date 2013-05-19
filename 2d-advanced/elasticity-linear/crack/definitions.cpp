#include "definitions.h"

CustomWeakFormLinearElasticity::CustomWeakFormLinearElasticity(double E, double nu, double rho_g,
                                 std::string surface_force_bdy, double f0, double f1) : WeakForm<double>(2)
{
  double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
  double mu = E / (2*(1 + nu));

  // SINGLE-COMPONENT FORMS. USEFUL FOR MULTIMESH, DO NOT REMOVE.
  // Jacobian.
  add_matrix_form(new DefaultJacobianElasticity_0_0<double>(0, 0, lambda, mu));
  add_matrix_form(new DefaultJacobianElasticity_0_1<double>(0, 1, lambda, mu));
  add_matrix_form(new DefaultJacobianElasticity_1_1<double>(1, 1, lambda, mu));

  // Residual - first equation.
  add_vector_form(new DefaultResidualElasticity_0_0<double>(0, HERMES_ANY, lambda, mu));
  add_vector_form(new DefaultResidualElasticity_0_1<double>(0, HERMES_ANY, lambda, mu));
  // Surface force (first component).
  add_vector_form_surf(new DefaultVectorFormSurf<double>(0, surface_force_bdy, new Hermes2DFunction<double>(-f0))); 

  // Residual - second equation.
  add_vector_form(new DefaultResidualElasticity_1_0<double>(1, HERMES_ANY, lambda, mu));
  add_vector_form(new DefaultResidualElasticity_1_1<double>(1, HERMES_ANY, lambda, mu));
  // Gravity loading in the second vector component.
  add_vector_form(new DefaultVectorFormVol<double>(1, HERMES_ANY, new Hermes2DFunction<double>(-rho_g)));
  // Surface force (second component).
  add_vector_form_surf(new DefaultVectorFormSurf<double>(1, surface_force_bdy, new Hermes2DFunction<double>(-f1))); 
}
