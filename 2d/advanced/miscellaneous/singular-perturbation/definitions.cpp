#include "definitions.h"

CustomWeakForm::CustomWeakForm(double K_squared) : WeakForm(1) 
{
  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(0, 0, HERMES_ANY, new Hermes2DFunction<double>(K_squared)));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0));
  add_vector_form(new WeakFormsH1::DefaultResidualVol<double>(0, HERMES_ANY, new Hermes2DFunction<double>(K_squared)));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(0, HERMES_ANY, new Hermes2DFunction<double>(-K_squared)));
}


