#include "definitions.h"

CustomWeakFormAcoustics::CustomWeakFormAcoustics(std::string bdy_newton, double rho, double sound_speed, double omega) : WeakForm<std::complex<double> >(1) 
{
  std::complex<double>  ii =  std::complex<double>(0.0, 1.0);

  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<std::complex<double> >(0, 0, HERMES_ANY, new Hermes1DFunction<std::complex<double> >(1.0/rho), HERMES_SYM));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<std::complex<double> >(0, 0, HERMES_ANY, new Hermes2DFunction<std::complex<double> >(-sqr(omega) / rho / sqr(sound_speed)), HERMES_SYM));
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf<std::complex<double> >(0, 0, bdy_newton, new Hermes2DFunction<std::complex<double> >(-ii * omega / rho / sound_speed)));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<std::complex<double> >(0, HERMES_ANY, new Hermes1DFunction<std::complex<double> >(1.0/rho)));
  add_vector_form(new WeakFormsH1::DefaultResidualVol<std::complex<double> >(0, HERMES_ANY, new Hermes2DFunction<std::complex<double> >(-sqr(omega) / rho / sqr(sound_speed))));
  add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf<std::complex<double> >(0, bdy_newton, new Hermes2DFunction<std::complex<double> >(-ii * omega / rho / sound_speed)));
}

