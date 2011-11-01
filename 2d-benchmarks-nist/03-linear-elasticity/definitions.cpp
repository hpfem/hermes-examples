#include "definitions.h"

double E = 1.0;
double nu = 0.3;
// lambda for mode-1 solution.
double lambda = 0.5444837367825;
// mu for mode-1 solution (mu is the same as G)
double mu = E / (2.0 * (1.0 + nu));
// Q for mode-1 solution.
double Q = 0.5430755788367;
double k = 3.0 - 4.0 * nu;
double A = -E * (1 - nu * nu)/(1 - 2 * nu);
double B = -E * (1 - nu * nu)/(2 - 2 * nu);
double C = -E * (1 - nu * nu)/((1 - 2 * nu) * (2 - 2 * nu));
double D = 1.0 / (2 * mu);
double u_F = (k - Q * (lambda + 1));
double v_F = (k + Q * (lambda + 1));

double CustomExactSolutionU::get_angle(double y, double x) const 
{
  double theta = Hermes::atan2(y, x);
  if (theta < 0)
    theta += 2 * M_PI;
  return theta;
}

double CustomExactSolutionU::d_theta_dx(double x, double y) const 
{
  return -y/(x*x + y*y);
}

double CustomExactSolutionU::d_theta_dxd_theta_dx(double x, double y) const 
{
  return 2*x*y/((x*x + y*y)*(x*x + y*y));
}

double CustomExactSolutionU::d_theta_dy(double x, double y) const 
{
  return x/(x*x + y*y) ;
}

double CustomExactSolutionU::d_theta_dyd_theta_dy(double x, double y) const 
{
  return -2*x*y/((x*x + y*y)*(x*x + y*y));
}

double CustomExactSolutionU::d_theta_dxd_theta_dy(double x, double y) const 
{
  return (y*y - x*x)/((x*x + y*y)*(x*x + y*y));
}

double CustomExactSolutionU::r(double x, double y) const 
{
  return Hermes::pow((x*x + y*y), (lambda/2.0));  // r^labbda
}

double CustomExactSolutionU::drdx(double x, double y) const 
{
  return lambda * x * Hermes::pow((x*x + y*y), (lambda/2.0 - 1.0));
}

double CustomExactSolutionU::drdxdrdx(double x, double y) const 
{
  return lambda * (Hermes::pow((x*x + y*y), (lambda/2.0 - 1.0)) + (lambda - 2.0) 
         * x * x * Hermes::pow((x*x + y*y), (lambda/2.0 - 2.0)));
}

double CustomExactSolutionU::drdy(double x, double y) const 
{
  return lambda * y * Hermes::pow((x*x + y*y), (lambda/2.0 - 1.0));
}

double CustomExactSolutionU::drdydrdy(double x, double y) const 
{
  return lambda * (Hermes::pow((x*x + y*y), (lambda/2.0 - 1.0)) + (lambda - 2.0) 
         * y * y * Hermes::pow((x*x + y*y), (lambda/2.0 - 2.0)));
}

double CustomExactSolutionU::drdxdrdy(double x, double y) const 
{
  return lambda * 2.0 * x * y * (lambda/2.0 - 1) * Hermes::pow((x*x + y*y), 
         (lambda/2.0 - 2.0));
}

double CustomExactSolutionU::u_r(double x, double y) const 
{
  return (u_F * Hermes::cos(lambda * get_angle(y, x)) - lambda * Hermes::cos((lambda - 2) 
         * get_angle(y, x)));
}

double CustomExactSolutionU::du_rdx(double x, double y) const 
{
  return (u_F * (-1) * lambda * Hermes::sin(lambda * get_angle(y, x)) * d_theta_dx(x, y)) 
         - (lambda * (-1) * (lambda - 2) * Hermes::sin((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));
}

double CustomExactSolutionU::du_rdxdu_rdx(double x, double y) const 
{
  return (u_F * (-1) * lambda * (Hermes::sin(lambda * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
         + lambda * d_theta_dx(x, y) * d_theta_dx(x, y) * Hermes::cos(lambda * get_angle(y, x)))) 
         - (lambda * (-1) * (lambda - 2) * (Hermes::sin((lambda - 2) * get_angle(y, x)) 
         * d_theta_dxd_theta_dx(x, y) + (lambda - 2) * d_theta_dx(x, y) * d_theta_dx(x, y) 
         * Hermes::cos((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionU::du_rdy(double x, double y) const 
{
  return (u_F * (-1) * lambda * Hermes::sin(lambda * get_angle(y, x)) * d_theta_dy(x, y)) 
         - (lambda * (-1) * (lambda - 2) * Hermes::sin((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
}

double CustomExactSolutionU::du_rdydu_rdy(double x, double y) const 
{
  return (u_F * (-1) * lambda * (Hermes::sin(lambda * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
         + lambda * d_theta_dy(x, y) * d_theta_dy(x, y) * Hermes::cos(lambda * get_angle(y, x)))) 
         - (lambda * (-1) * (lambda - 2) * (Hermes::sin((lambda - 2) * get_angle(y, x)) 
         * d_theta_dyd_theta_dy(x, y) + (lambda - 2) * d_theta_dy(x, y) * d_theta_dy(x, y) 
         * Hermes::cos((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionU::du_rdxdu_rdy(double x, double y) const 
{
  return (u_F * (-1) * lambda * (Hermes::sin(lambda * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
         + lambda * d_theta_dx(x, y) * d_theta_dy(x, y) * Hermes::cos(lambda * get_angle(y, x)))) 
         - (lambda * (-1) * (lambda - 2) * (Hermes::sin((lambda - 2) * get_angle(y, x)) 
         * d_theta_dxd_theta_dy(x, y) + (lambda - 2) * d_theta_dx(x, y) * d_theta_dy(x, y) 
         * Hermes::cos((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionU::v_r(double x, double y) const 
{
  return (v_F * Hermes::sin(lambda * get_angle(y, x)) + lambda * Hermes::sin((lambda - 2) * get_angle(y, x)));
}

double CustomExactSolutionU::dv_rdx(double x, double y) const 
{
  return (v_F * lambda * Hermes::cos(lambda * get_angle(y, x)) * d_theta_dx(x, y)) 
         + (lambda * (lambda - 2) * Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));
}

double CustomExactSolutionU::dv_rdxdv_rdx(double x, double y) const 
{
  return (v_F * lambda * (Hermes::cos(lambda * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
         + lambda * d_theta_dx(x, y) * d_theta_dx(x, y) * (-1) * Hermes::sin(lambda * get_angle(y, x)))) 
         + (lambda * (lambda - 2) * (Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
         + (lambda - 2) * d_theta_dx(x, y) * d_theta_dx(x, y) * (-1) * Hermes::sin((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionU::dv_rdy(double x, double y) const 
{
  return (v_F * lambda * Hermes::cos(lambda * get_angle(y, x)) * d_theta_dy(x, y)) + (lambda * (lambda - 2) 
         * Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
}

double CustomExactSolutionU::dv_rdydv_rdy(double x, double y) const 
{
  return (v_F * lambda * (Hermes::cos(lambda * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
         + lambda * d_theta_dy(x, y) * d_theta_dy(x, y) * (-1) * Hermes::sin(lambda * get_angle(y, x)))) 
         + (lambda * (lambda - 2) * (Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
         + (lambda - 2) * d_theta_dy(x, y) * d_theta_dy(x, y) * (-1) * Hermes::sin((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionU::dv_rdxdv_rdy(double x, double y) const 
{
  return (v_F * lambda * (Hermes::cos(lambda * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
         + lambda * (-1) * d_theta_dx(x, y) * d_theta_dy(x, y) * Hermes::sin(lambda * get_angle(y, x)))) 
         + (lambda * (lambda - 2) * (Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
         + (lambda - 2) * (-1) * d_theta_dx(x, y) * d_theta_dy(x, y) * Hermes::sin((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionU::dudxdudx(double x, double y) const 
{
  return D * (drdxdrdx(x, y) * u_r(x, y) + 2 * drdx(x, y) * du_rdx(x, y) + r(x, y) * du_rdxdu_rdx(x, y) );
}

double CustomExactSolutionU::dudydudy(double x, double y) const 
{
  return D * (drdydrdy(x, y) * u_r(x, y) + 2 * drdy(x, y) * du_rdy(x, y) + r(x, y) * du_rdydu_rdy(x, y) );
}

double CustomExactSolutionU::dudxdudy(double x, double y) const 
{
  return D * (drdxdrdy(x, y) * u_r(x, y) + drdx(x, y) * du_rdy(x, y) + drdy(x, y) * du_rdx(x, y) + r(x, y) 
         * du_rdxdu_rdy(x, y) );
}

double CustomExactSolutionU::dvdxdvdx(double x, double y) const 
{
  return D * (drdxdrdx(x, y) * v_r(x, y) + 2 * drdx(x, y) * dv_rdx(x, y) + r(x, y) * dv_rdxdv_rdx(x, y) );
}

double CustomExactSolutionU::dvdydvdy(double x, double y) const 
{
  return D * (drdydrdy(x, y) * v_r(x, y) + 2 * drdy(x, y) * dv_rdy(x, y) + r(x, y) * dv_rdydv_rdy(x, y) );
}

double CustomExactSolutionU::dvdxdvdy(double x, double y) const 
{
  return D * (drdxdrdy(x, y) * v_r(x, y) + drdx(x, y) * dv_rdy(x, y) + drdy(x, y) * dv_rdx(x, y) + r(x, y) 
         * dv_rdxdv_rdy(x, y) );
}

double CustomExactSolutionU::value (double x, double y) const 
{
  return D * r(x, y) * u_r(x, y);
}

void CustomExactSolutionU::derivatives (double x, double y, double& dx, double& dy) const 
{
  dx = D * drdx(x, y) * (u_F * Hermes::cos(lambda * get_angle(y, x)) - lambda * Hermes::cos((lambda - 2) * get_angle(y, x))) +
       D * r(x, y) * (u_F * (-1) * lambda * Hermes::sin(lambda * get_angle(y, x)) * d_theta_dx(x, y)) - 
       D * r(x, y) * (lambda * (-1) * (lambda - 2) * Hermes::sin((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));

  dy = D * drdy(x, y) * (u_F * Hermes::cos(lambda * get_angle(y, x)) - lambda * Hermes::cos((lambda - 2) * get_angle(y, x))) +
       D * r(x, y) * (u_F * (-1) * lambda * Hermes::sin(lambda * get_angle(y, x)) * d_theta_dy(x, y)) - 
       D * r(x, y) * (lambda * (-1) * (lambda - 2) * Hermes::sin((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
}

Ord CustomExactSolutionU::ord (Ord x, Ord y) const 
{
  return Ord(4.0);
}

double CustomExactSolutionV::get_angle(double y, double x) const 
{
  double theta = Hermes::atan2(y, x);
  if (theta < 0)
    theta += 2 * M_PI;
  return theta;
}

double CustomExactSolutionV::d_theta_dx(double x, double y) const 
{
  return -y/(x*x + y*y);
}

double CustomExactSolutionV::d_theta_dxd_theta_dx(double x, double y) const 
{
  return 2*x*y/((x*x + y*y)*(x*x + y*y));
}

double CustomExactSolutionV::d_theta_dy(double x, double y) const 
{
  return x/(x*x + y*y) ;
}

double CustomExactSolutionV::d_theta_dyd_theta_dy(double x, double y) const 
{
  return -2*x*y/((x*x + y*y)*(x*x + y*y));
}

double CustomExactSolutionV::d_theta_dxd_theta_dy(double x, double y) const 
{
  return (y*y - x*x)/((x*x + y*y)*(x*x + y*y));
}

double CustomExactSolutionV::r(double x, double y) const 
{
  return Hermes::pow((x*x + y*y), (lambda/2.0));  // r^labbda
}

double CustomExactSolutionV::drdx(double x, double y) const 
{
  return lambda * x * Hermes::pow((x*x + y*y), (lambda/2.0 - 1.0));
}

double CustomExactSolutionV::drdxdrdx(double x, double y) const 
{
  return lambda * (Hermes::pow((x*x + y*y), (lambda/2.0 - 1.0)) + (lambda - 2.0) 
         * x * x * Hermes::pow((x*x + y*y), (lambda/2.0 - 2.0)));
}

double CustomExactSolutionV::drdy(double x, double y) const 
{
  return lambda * y * Hermes::pow((x*x + y*y), (lambda/2.0 - 1.0));
}

double CustomExactSolutionV::drdydrdy(double x, double y) const 
{
  return lambda * (Hermes::pow((x*x + y*y), (lambda/2.0 - 1.0)) + (lambda - 2.0) 
         * y * y * Hermes::pow((x*x + y*y), (lambda/2.0 - 2.0)));
}

double CustomExactSolutionV::drdxdrdy(double x, double y) const 
{
  return lambda * 2.0 * x * y * (lambda/2.0 - 1) * Hermes::pow((x*x + y*y), 
         (lambda/2.0 - 2.0));
}

double CustomExactSolutionV::u_r(double x, double y) const 
{
  return (u_F * Hermes::cos(lambda * get_angle(y, x)) - lambda * Hermes::cos((lambda - 2) 
         * get_angle(y, x)));
}

double CustomExactSolutionV::du_rdx(double x, double y) const 
{
  return (u_F * (-1) * lambda * Hermes::sin(lambda * get_angle(y, x)) * d_theta_dx(x, y)) 
         - (lambda * (-1) * (lambda - 2) * Hermes::sin((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));
}

double CustomExactSolutionV::du_rdxdu_rdx(double x, double y) const 
{
  return (u_F * (-1) * lambda * (Hermes::sin(lambda * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
         + lambda * d_theta_dx(x, y) * d_theta_dx(x, y) * Hermes::cos(lambda * get_angle(y, x)))) 
         - (lambda * (-1) * (lambda - 2) * (Hermes::sin((lambda - 2) * get_angle(y, x)) 
         * d_theta_dxd_theta_dx(x, y) + (lambda - 2) * d_theta_dx(x, y) * d_theta_dx(x, y) 
         * Hermes::cos((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionV::du_rdy(double x, double y) const 
{
  return (u_F * (-1) * lambda * Hermes::sin(lambda * get_angle(y, x)) * d_theta_dy(x, y)) 
         - (lambda * (-1) * (lambda - 2) * Hermes::sin((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
}

double CustomExactSolutionV::du_rdydu_rdy(double x, double y) const 
{
  return (u_F * (-1) * lambda * (Hermes::sin(lambda * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
         + lambda * d_theta_dy(x, y) * d_theta_dy(x, y) * Hermes::cos(lambda * get_angle(y, x)))) 
         - (lambda * (-1) * (lambda - 2) * (Hermes::sin((lambda - 2) * get_angle(y, x)) 
         * d_theta_dyd_theta_dy(x, y) + (lambda - 2) * d_theta_dy(x, y) * d_theta_dy(x, y) 
         * Hermes::cos((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionV::du_rdxdu_rdy(double x, double y) const 
{
  return (u_F * (-1) * lambda * (Hermes::sin(lambda * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
         + lambda * d_theta_dx(x, y) * d_theta_dy(x, y) * Hermes::cos(lambda * get_angle(y, x)))) 
         - (lambda * (-1) * (lambda - 2) * (Hermes::sin((lambda - 2) * get_angle(y, x)) 
         * d_theta_dxd_theta_dy(x, y) + (lambda - 2) * d_theta_dx(x, y) * d_theta_dy(x, y) 
         * Hermes::cos((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionV::v_r(double x, double y) const 
{
  return (v_F * Hermes::sin(lambda * get_angle(y, x)) + lambda * Hermes::sin((lambda - 2) * get_angle(y, x)));
}

double CustomExactSolutionV::dv_rdx(double x, double y) const 
{
  return (v_F * lambda * Hermes::cos(lambda * get_angle(y, x)) * d_theta_dx(x, y)) 
         + (lambda * (lambda - 2) * Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));
}

double CustomExactSolutionV::dv_rdxdv_rdx(double x, double y) const 
{
  return (v_F * lambda * (Hermes::cos(lambda * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
         + lambda * d_theta_dx(x, y) * d_theta_dx(x, y) * (-1) * Hermes::sin(lambda * get_angle(y, x)))) 
         + (lambda * (lambda - 2) * (Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
         + (lambda - 2) * d_theta_dx(x, y) * d_theta_dx(x, y) * (-1) * Hermes::sin((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionV::dv_rdy(double x, double y) const 
{
  return (v_F * lambda * Hermes::cos(lambda * get_angle(y, x)) * d_theta_dy(x, y)) + (lambda * (lambda - 2) 
         * Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
}

double CustomExactSolutionV::dv_rdydv_rdy(double x, double y) const 
{
  return (v_F * lambda * (Hermes::cos(lambda * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
         + lambda * d_theta_dy(x, y) * d_theta_dy(x, y) * (-1) * Hermes::sin(lambda * get_angle(y, x)))) 
         + (lambda * (lambda - 2) * (Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
         + (lambda - 2) * d_theta_dy(x, y) * d_theta_dy(x, y) * (-1) * Hermes::sin((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionV::dv_rdxdv_rdy(double x, double y) const 
{
  return (v_F * lambda * (Hermes::cos(lambda * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
         + lambda * (-1) * d_theta_dx(x, y) * d_theta_dy(x, y) * Hermes::sin(lambda * get_angle(y, x)))) 
         + (lambda * (lambda - 2) * (Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
         + (lambda - 2) * (-1) * d_theta_dx(x, y) * d_theta_dy(x, y) * Hermes::sin((lambda-2) * get_angle(y, x)) ));
}

double CustomExactSolutionV::dudxdudx(double x, double y) const 
{
  return D * (drdxdrdx(x, y) * u_r(x, y) + 2 * drdx(x, y) * du_rdx(x, y) + r(x, y) * du_rdxdu_rdx(x, y) );
}

double CustomExactSolutionV::dudydudy(double x, double y) const 
{
  return D * (drdydrdy(x, y) * u_r(x, y) + 2 * drdy(x, y) * du_rdy(x, y) + r(x, y) * du_rdydu_rdy(x, y) );
}

double CustomExactSolutionV::dudxdudy(double x, double y) const 
{
  return D * (drdxdrdy(x, y) * u_r(x, y) + drdx(x, y) * du_rdy(x, y) + drdy(x, y) * du_rdx(x, y) + r(x, y) 
         * du_rdxdu_rdy(x, y) );
}

double CustomExactSolutionV::dvdxdvdx(double x, double y) const 
{
  return D * (drdxdrdx(x, y) * v_r(x, y) + 2 * drdx(x, y) * dv_rdx(x, y) + r(x, y) * dv_rdxdv_rdx(x, y) );
}

double CustomExactSolutionV::dvdydvdy(double x, double y) const 
{
  return D * (drdydrdy(x, y) * v_r(x, y) + 2 * drdy(x, y) * dv_rdy(x, y) + r(x, y) * dv_rdydv_rdy(x, y) );
}

double CustomExactSolutionV::dvdxdvdy(double x, double y) const 
{
  return D * (drdxdrdy(x, y) * v_r(x, y) + drdx(x, y) * dv_rdy(x, y) + drdy(x, y) * dv_rdx(x, y) + r(x, y) 
         * dv_rdxdv_rdy(x, y) );
}

double CustomExactSolutionV::value (double x, double y) const 
{
  return D * r(x, y) * v_r(x, y);
}

void CustomExactSolutionV::derivatives (double x, double y, double& dx, double& dy) const 
{
    dx = D * drdx(x, y) * (v_F * Hermes::sin(lambda * get_angle(y, x)) + lambda * Hermes::sin((lambda - 2) * get_angle(y, x))) +
         D * r(x, y) * (v_F * lambda * Hermes::cos(lambda * get_angle(y, x)) * d_theta_dx(x, y)) + 
         D * r(x, y) * (lambda * (lambda - 2) * Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));

    dy = D * drdy(x, y) * (v_F * Hermes::sin(lambda * get_angle(y, x)) + lambda * Hermes::sin((lambda - 2) * get_angle(y, x))) +
         D * r(x, y) * (v_F * lambda * Hermes::cos(lambda * get_angle(y, x)) * d_theta_dy(x, y)) + 
         D * r(x, y) * (lambda * (lambda - 2) * Hermes::cos((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
}

Ord CustomExactSolutionV::ord (Ord x, Ord y) const 
{
  return Ord(4.0);
}


CustomWeakFormElasticityNIST::CustomWeakFormElasticityNIST(double E, double nu, double mu, double lambda) : WeakForm<double>(2) 
{
  // Jacobian.
  add_matrix_form(new CustomMatrixFormVolElasticityNIST_0_0(0, 0, E, nu));
  add_matrix_form(new CustomMatrixFormVolElasticityNIST_0_1(0, 1, E, nu));
  add_matrix_form(new CustomMatrixFormVolElasticityNIST_1_1(1, 1, E, nu));
  // Residual.
  add_vector_form(new CustomVectorFormVolElasticityNIST_0(0, E, nu));
  add_vector_form(new CustomVectorFormVolElasticityNIST_1(1, E, nu));
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormElasticityNIST::CustomMatrixFormVolElasticityNIST_0_0::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],     
    Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar val = Scalar(0);
  for (int i = 0; i < n; i++)
    val += wt[i] * (A * u->dx[i] * v->dx[i] + B * u->dy[i] * v->dy[i]);
  return val;
}

double CustomWeakFormElasticityNIST::CustomMatrixFormVolElasticityNIST_0_0::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormElasticityNIST::CustomMatrixFormVolElasticityNIST_0_0::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormElasticityNIST::CustomMatrixFormVolElasticityNIST_0_1::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],     
    Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar val = Scalar(0);
  for (int i = 0; i < n; i++)
    val += wt[i] * (C * u->dx[i] * v->dy[i]);
  return val;
}

double CustomWeakFormElasticityNIST::CustomMatrixFormVolElasticityNIST_0_1::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormElasticityNIST::CustomMatrixFormVolElasticityNIST_0_1::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormElasticityNIST::CustomMatrixFormVolElasticityNIST_1_1::matrix_form(int n, double *wt, Func<Scalar> *u_ext[],     
    Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar val = Scalar(0);
  for (int i = 0; i < n; i++)
    val += wt[i] * (B * u->dx[i] * v->dx[i] + A * u->dy[i] * v->dy[i]);
  return val;
}

double CustomWeakFormElasticityNIST::CustomMatrixFormVolElasticityNIST_1_1::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
}

Ord CustomWeakFormElasticityNIST::CustomMatrixFormVolElasticityNIST_1_1::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormElasticityNIST::CustomVectorFormVolElasticityNIST_0::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar val = Scalar(0);
  for (int i = 0; i < n; i++) 
  {
    // Contribution of matrix form 0, 0.
    val += wt[i] * (A * u_ext[0]->dx[i] * v->dx[i] + B * u_ext[0]->dy[i] * v->dy[i]);
    // Contribution of matrix form 0, 1.
    val += wt[i] * (C * u_ext[1]->dx[i] * v->dy[i]);
  }

  return val;
}

double CustomWeakFormElasticityNIST::CustomVectorFormVolElasticityNIST_0::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormElasticityNIST::CustomVectorFormVolElasticityNIST_0::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}

template<typename Real, typename Scalar>
Scalar CustomWeakFormElasticityNIST::CustomVectorFormVolElasticityNIST_1::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
{
  Scalar val = Scalar(0);
  for (int i = 0; i < n; i++) 
  {
    // Contribution of matrix form 1, 0.
    val += wt[i] * (C * u_ext[0]->dy[i] * v->dx[i]);
    // Contribution of matrix form 1, 1.
    val += wt[i] * (B * u_ext[1]->dx[i] * v->dx[i] + A * u_ext[1]->dy[i] * v->dy[i]);
  }

  return val;
}

double CustomWeakFormElasticityNIST::CustomVectorFormVolElasticityNIST_1::value(int n, double *wt, Func<double> *u_ext[],
    Func<double> *v, Geom<double> *e, ExtData<double> *ext) const 
{
  return vector_form<double, double>(n, wt, u_ext, v, e, ext);
}

Ord CustomWeakFormElasticityNIST::CustomVectorFormVolElasticityNIST_1::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
}
