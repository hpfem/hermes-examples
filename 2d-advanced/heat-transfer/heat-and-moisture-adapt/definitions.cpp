#include "definitions.h"

CustomWeakFormHeatMoistureRK::CustomWeakFormHeatMoistureRK(double c_TT, double c_ww, double d_TT, double d_Tw, double d_wT, 
    double d_ww, double k_TT, double k_ww, double T_ext, double w_ext, const std::string bdy_ext) : WeakForm<double>(2)
{
  // Jacobian - volumetric.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, HERMES_ANY, new Hermes1DFunction<double>(-d_TT/c_TT), HERMES_SYM, HERMES_AXISYM_Y));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 1, HERMES_ANY, new Hermes1DFunction<double>(-d_Tw/c_TT), HERMES_NONSYM, HERMES_AXISYM_Y));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(1, 0, HERMES_ANY, new Hermes1DFunction<double>(-d_wT/c_ww), HERMES_NONSYM, HERMES_AXISYM_Y));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(1, 1, HERMES_ANY, new Hermes1DFunction<double>(-d_ww/c_ww), HERMES_SYM, HERMES_AXISYM_Y));

  // Jacobian - surface.
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf<double>(0, 0, bdy_ext, new Hermes2DFunction<double>(-k_TT/c_TT), HERMES_AXISYM_Y));
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf<double>(1, 1, bdy_ext, new Hermes2DFunction<double>(-k_ww/c_ww), HERMES_AXISYM_Y));

  // Residual - volumetric
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, HERMES_ANY, new Hermes1DFunction<double>(-d_TT/c_TT), HERMES_AXISYM_Y));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, HERMES_ANY, new Hermes1DFunction<double>(-d_Tw/c_TT), HERMES_AXISYM_Y));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(1, HERMES_ANY, new Hermes1DFunction<double>(-d_wT/c_ww), HERMES_AXISYM_Y));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(1, HERMES_ANY, new Hermes1DFunction<double>(-d_ww/c_ww), HERMES_AXISYM_Y));

  // Residual - surface.
  add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf<double>(0, bdy_ext, new Hermes2DFunction<double>(k_TT/c_TT * T_ext), HERMES_AXISYM_Y));
  add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf<double>(0, bdy_ext, new Hermes2DFunction<double>(-k_TT/c_TT), HERMES_AXISYM_Y));
  add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf<double>(1, bdy_ext, new Hermes2DFunction<double>(k_ww/c_ww * w_ext), HERMES_AXISYM_Y));
  add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf<double>(1, bdy_ext, new Hermes2DFunction<double>(-k_ww/c_ww), HERMES_AXISYM_Y));
}

EssentialBCNonConst::EssentialBCNonConst(std::string marker, double reactor_start_time, double temp_initial, double temp_reactor_max) 
    : EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), reactor_start_time(reactor_start_time), 
    temp_initial(temp_initial), temp_reactor_max(temp_reactor_max)  
{
  markers.push_back(marker);
}

EssentialBoundaryCondition<double>::EssentialBCValueType EssentialBCNonConst::get_value_type() const 
{ 
  return BC_FUNCTION; 
}

double EssentialBCNonConst::value(double x, double y, double n_x, double n_y, double t_x, double t_y) const 
{
  double current_reactor_temperature = temp_reactor_max;
  if (current_time < reactor_start_time) 
  {
    current_reactor_temperature = temp_initial + (current_time/reactor_start_time) * (temp_reactor_max - temp_initial);
  }
  return current_reactor_temperature;
}

