#include "definitions.h"

CustomWeakFormHeatMoistureRK::CustomWeakFormHeatMoistureRK(double c_TT, double c_ww, double d_TT, double d_Tw, double d_wT, 
    double d_ww, double k_TT, double k_ww, double T_ext, double w_ext, const std::string bdy_ext) : WeakForm<double>(2)
{
  // Jacobian - volumetric.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, HERMES_ANY, new Hermes1DFunction<double>(d_TT/c_TT), HERMES_AXISYM_Y));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 1, HERMES_ANY, new Hermes1DFunction<double>(d_Tw/c_TT), HERMES_AXISYM_Y));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(1, 0, HERMES_ANY, new Hermes1DFunction<double>(d_wT/c_ww), HERMES_AXISYM_Y));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(1, 1, HERMES_ANY, new Hermes1DFunction<double>(d_ww/c_ww), HERMES_AXISYM_Y));

  // Jacobian - surface.
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf<double>(0, 0, bdy_ext, new Hermes2DFunction<double>(-k_TT/c_TT), HERMES_AXISYM_Y));
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf<double>(1, 1, bdy_ext, new Hermes2DFunction<double>(-k_ww/c_ww), HERMES_AXISYM_Y));

  // Residual - volumetric
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, 0, HERMES_ANY, new Hermes1DFunction<double>(d_TT/c_TT), HERMES_AXISYM_Y));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, 1, HERMES_ANY, new Hermes1DFunction<double>(d_Tw/c_TT), HERMES_AXISYM_Y));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(1, 0, HERMES_ANY, new Hermes1DFunction<double>(d_wT/c_ww), HERMES_AXISYM_Y));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(1, 1, HERMES_ANY, new Hermes1DFunction<double>(d_ww/c_ww), HERMES_AXISYM_Y));

  // Residual - surface.
  add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf<double>(0, bdy_ext, 
      new Hermes2DFunction<double>(k_TT/c_TT * T_ext), HERMES_AXISYM_Y));
  add_vector_form_surf(new WeakFormsH1::DefaultResidualFormSurf<double>(0, bdy_ext, 
      new Hermes2DFunction<double>(-k_TT/c_TT), HERMES_AXISYM_Y));
  add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf<double>(1, bdy_ext, 
      new Hermes2DFunction<double>(k_ww/c_ww * w_ext), HERMES_AXISYM_Y));
  add_vector_form_surf(new WeakFormsH1::DefaultResidualFormSurf<double>(1, bdy_ext, 
      new Hermes2DFunction<double>(-k_ww/c_ww), HERMES_AXISYM_Y));
}
