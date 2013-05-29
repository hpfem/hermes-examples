#include "definitions.h"

CustomWeakFormAcoustics::CustomWeakFormAcoustics(std::string bdy_newton, double rho,
                                                 double sound_speed, double omega) : WeakForm<complex>(1) 
{
  std::complex<double>  ii =  std::complex<double>(0.0, 1.0);

  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<complex>(0, 0, HERMES_ANY, new Hermes1DFunction<complex>(1.0/rho), HERMES_SYM));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<complex>(0, 0, HERMES_ANY, new Hermes2DFunction<complex>(-sqr(omega) / rho / sqr(sound_speed)), HERMES_SYM));
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf<complex>(0, 0, bdy_newton, new Hermes2DFunction<complex>(-ii * omega / rho / sound_speed)));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<complex>(0, HERMES_ANY, new Hermes1DFunction<complex>(1.0/rho)));
  add_vector_form(new WeakFormsH1::DefaultResidualVol<complex>(0, HERMES_ANY, new Hermes2DFunction<complex>(-sqr(omega) / rho / sqr(sound_speed))));
  add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf<complex>(0, bdy_newton, new Hermes2DFunction<complex>(-ii * omega / rho / sound_speed)));
}
