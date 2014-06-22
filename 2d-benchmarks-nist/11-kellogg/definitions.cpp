#include "definitions.h"

double CustomExactSolution::value(double x, double y) const 
{
  double theta = Hermes::atan2(y,x);
  if (theta < 0) theta = theta + 2.*M_PI;
  double r = Hermes::sqrt(x*x + y*y);
  double mu;

  if (theta <= M_PI/2.)
    mu = Hermes::cos((M_PI/2. - sigma)*tau) * Hermes::cos((theta - M_PI/2. + rho)*tau);
  else if (theta <= M_PI)
    mu = Hermes::cos(rho*tau) * Hermes::cos((theta - M_PI + sigma)*tau);
  else if (theta <= 3.*M_PI/2.)
    mu = Hermes::cos(sigma*tau) * Hermes::cos((theta - M_PI - rho)*tau);
  else
    mu = Hermes::cos((M_PI/2. - rho)*tau) * Hermes::cos((theta - 3.*M_PI/2. - sigma)*tau);

  return Hermes::pow(r, tau) * mu;
}

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  double theta = Hermes::atan2(y,x);
  if (theta < 0) theta = theta + 2*M_PI;
  double r = Hermes::sqrt(x*x + y*y);
    
  // x-derivative
  if (theta <= M_PI/2.) 
    dx = tau*x*Hermes::pow(r, (2.*(-1 + tau/2.))) * Hermes::cos((M_PI/2. - sigma)*tau) * Hermes::cos(tau*(-M_PI/2. + rho + theta)) 
        + (tau*y*Hermes::pow(r, tau)*Hermes::cos((M_PI/2. - sigma)*tau) * Hermes::sin(tau*(-M_PI/2. + rho + theta))/(r*r));
  else if (theta <= M_PI)
    dx = tau*x * Hermes::pow(r, (2.*(-1 + tau/2.))) * Hermes::cos(rho*tau) * Hermes::cos(tau*(-M_PI + sigma + theta)) 
        + (tau*y * Hermes::pow(r, tau) * Hermes::cos(rho*tau) * Hermes::sin(tau*(-M_PI + sigma + theta))/(r*r));
  else if (theta <= 3.*M_PI/2.)
    dx = tau*x * Hermes::pow(r, (2.*(-1 + tau/2.))) * Hermes::cos(sigma*tau) * Hermes::cos(tau*(-M_PI - rho + theta)) 
        + (tau*y * Hermes::pow(r, tau) * Hermes::cos(sigma*tau) * Hermes::sin(tau*(-M_PI - rho + theta))/(r*r));
  else
    dx = tau*x* Hermes::pow(r, (2*(-1 + tau/2.))) * Hermes::cos((M_PI/2. - rho)*tau) * Hermes::cos(tau*(-3.*M_PI/2. - sigma + theta)) 
        + (tau*y*Hermes::pow(r, tau) * Hermes::cos((M_PI/2. - rho)*tau) * Hermes::sin(tau*(-3.*M_PI/2. - sigma + theta))/(r*r));
    
  // y-derivative
  if (theta <= M_PI/2.)
    dy = tau*y * Hermes::pow(r, (2*(-1 + tau/2.))) * Hermes::cos((M_PI/2. - sigma)*tau) * Hermes::cos(tau*(-M_PI/2. + rho + theta)) 
        - (tau * Hermes::pow(r, tau) * Hermes::cos((M_PI/2. - sigma)*tau) *Hermes::sin(tau*(-M_PI/2. + rho + theta))*x/(r*r));
  else if (theta <= M_PI)
    dy = tau*y* Hermes::pow(r, (2*(-1 + tau/2.))) * Hermes::cos(rho*tau) * Hermes::cos(tau*(-M_PI + sigma + theta)) 
        - (tau * Hermes::pow(r, tau) * Hermes::cos(rho*tau) * Hermes::sin(tau*(-M_PI + sigma + theta))*x/(r*r));
  else if (theta <= 3.*M_PI/2.)
    dy = tau*y * Hermes::pow(r, (2*(-1 + tau/2.))) * Hermes::cos(sigma*tau) * Hermes::cos(tau*(-M_PI - rho + theta)) 
        - (tau * Hermes::pow(r, tau) * Hermes::cos(sigma*tau) * Hermes::sin(tau*(-M_PI - rho + theta))*x/(r*r));
  else 
    dy = tau*y * Hermes::pow(r, (2*(-1 + tau/2.))) * Hermes::cos((M_PI/2. - rho)*tau) * Hermes::cos(tau*(-3.*M_PI/2. - sigma + theta)) 
        - (tau * Hermes::pow(r, tau) * Hermes::cos((M_PI/2. - rho)*tau) * Hermes::sin(tau*((-3.*M_PI)/2. - sigma + theta))*x/(r*r));
}

Ord CustomExactSolution::ord (double x, double y) const
{
  return Ord(6);
}

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string area_1, double r, std::string area_2) : WeakForm<double>(1)
{
  // Jacobian.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, area_1, new Hermes1DFunction<double>(r)));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(0, 0, area_2, nullptr));

  // Residual.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, area_1, new Hermes1DFunction<double>(r)));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(0, area_2, nullptr));
}
