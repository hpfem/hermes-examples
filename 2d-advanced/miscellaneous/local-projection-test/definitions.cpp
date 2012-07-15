#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string mat_al, Hermes1DFunction<double>* lambda_al,
                                             std::string mat_cu, Hermes1DFunction<double>* lambda_cu,
                                             Hermes2DFunction<double>* src_term) : WeakForm<double>(1)
{
  // Jacobian forms.
  add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, lambda_al, mat_al));
  add_matrix_form(new DefaultJacobianDiffusion<double>(0, 0, lambda_cu, mat_cu));

  // Residual forms.
  add_vector_form(new DefaultResidualDiffusion<double>(0, lambda_al, mat_al));
  add_vector_form(new DefaultResidualDiffusion<double>(0, lambda_cu, mat_cu));
  add_vector_form(new DefaultVectorFormVol<double>(0, src_term));
};
