#include "definitions.h"

/* Custom initial condition */

CustomInitialCondition::CustomInitialCondition(Mesh *mesh, double mid_x, double mid_y, double radius, double temp_water, double temp_graphite) 
  : ExactSolutionScalar<double>(mesh), mid_x(mid_x), mid_y(mid_y), radius(radius), temp_water(temp_water), temp_graphite(temp_graphite) {}

double CustomInitialCondition::value(double x, double y) const 
{
  bool in_graphite = (std::sqrt(sqr(mid_x - x) + sqr(mid_y - y)) < radius);
  double temp = temp_water;
  if (in_graphite) temp = temp_graphite;
  return temp;
}

void CustomInitialCondition::derivatives(double x, double y, double& dx, double& dy) const 
{   
  dx = 0;
  dy = 0;
}

Ord CustomInitialCondition::ord(Ord x, Ord y) const 
{
  return Ord(1);
}

/* Weak forms */

CustomWeakFormHeatAndFlow::CustomWeakFormHeatAndFlow(bool Stokes, double Reynolds, double time_step, Solution<double>* x_vel_previous_time, 
  Solution<double>* y_vel_previous_time, Solution<double>* T_prev_time, double heat_source, double specific_heat_graphite, 
  double specific_heat_water, double rho_graphite, double rho_water, double thermal_conductivity_graphite, double thermal_conductivity_water) 
  : WeakForm<double>(4), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step), x_vel_previous_time(x_vel_previous_time), y_vel_previous_time(y_vel_previous_time)
  {
    // Jacobian - flow part.
    add_matrix_form(new BilinearFormSymVel(0, 0, Stokes, Reynolds, time_step));
    add_matrix_form(new BilinearFormSymVel(1, 1, Stokes, Reynolds, time_step));
    add_matrix_form(new BilinearFormUnSymVel_0_0(0, 0, Stokes));
    add_matrix_form(new BilinearFormUnSymVel_0_1(0, 1, Stokes));
    add_matrix_form(new BilinearFormUnSymVel_1_0(1, 0, Stokes));
    add_matrix_form(new BilinearFormUnSymVel_1_1(1, 1, Stokes));
    add_matrix_form(new BilinearFormUnSymXVelPressure(0, 2));
    add_matrix_form(new BilinearFormUnSymYVelPressure(1, 2));

    // Jacobian - temperature part. 
    // Contribution from implicit Euler.
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(3, 3, "Water", new Hermes2DFunction<double>(1.0/time_step), HERMES_NONSYM));
    add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<double>(3, 3, "Graphite", new Hermes2DFunction<double>(1.0/time_step), HERMES_NONSYM));
    // Contribution from temperature diffusion. 
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(3, 3, "Water", new Hermes1DFunction<double>(thermal_conductivity_water/(rho_water * specific_heat_water)), HERMES_NONSYM));
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(3, 3, "Graphite", new Hermes1DFunction<double>(thermal_conductivity_graphite/(rho_graphite * specific_heat_graphite)), HERMES_NONSYM));
    // Contribution from temperature advection - only in water.
    add_matrix_form(new CustomJacobianTempAdvection_3_0(3, 0, "Water"));
    add_matrix_form(new CustomJacobianTempAdvection_3_1(3, 1, "Water"));
    add_matrix_form(new CustomJacobianTempAdvection_3_3(3, 3, "Water"));

    // Residual - flow part. 
    VectorFormNS_0* F_0 = new VectorFormNS_0(0, Stokes, Reynolds, time_step);
    F_0->ext.push_back(x_vel_previous_time);
    F_0->ext.push_back(y_vel_previous_time);
    add_vector_form(F_0);
    VectorFormNS_1* F_1 = new VectorFormNS_1(1, Stokes, Reynolds, time_step);
    F_1->ext.push_back(x_vel_previous_time);
    F_1->ext.push_back(y_vel_previous_time);
    add_vector_form(F_1);
    VectorFormNS_2* F_2 = new VectorFormNS_2(2);
    add_vector_form(F_2);

    // Residual - temperature part.
    // Contribution from implicit Euler method.
    VectorFormTime *vft = new VectorFormTime(3, "Water", time_step);
    vft->ext.push_back(T_prev_time);
    add_vector_form(vft);
    vft = new VectorFormTime(3, "Graphite", time_step);
    vft->ext.push_back(T_prev_time);
    add_vector_form(vft);
    // Contribution from temperature diffusion.
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(3, "Water", new Hermes1DFunction<double>(thermal_conductivity_water/(rho_water * specific_heat_water))));
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(3, "Graphite", new Hermes1DFunction<double>(thermal_conductivity_graphite/(rho_graphite * specific_heat_graphite))));
    // Contribution from heat sources.
    add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(3, "Graphite", new Hermes::Hermes2DFunction<double>(-heat_source/(rho_graphite * specific_heat_graphite))));
    // Contribution from temperature advection.
    add_vector_form(new CustomResidualTempAdvection(3, "Water"));
  };
