#include "definitions.h"

/* Custom initial condition */

CustomInitialCondition::CustomInitialCondition(Mesh *mesh, double mid_x, double mid_y, double radius, double temp_water, double temp_graphite) 
  : mid_x(mid_x), mid_y(mid_y), radius(radius), temp_water(temp_water), temp_graphite(temp_graphite) 
{
  return ;
}

double CustomInitialCondition::value(double x, double y) const 
{
  bool in_graphite = (sqrt(sqr(mid_x - x) + sqr(mid_y - y)) < radius);
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
  return 1;
}

/* Weak forms */

CustomWeakFormHeatAndFlow::CustomWeakFormHeatAndFlow(bool Stokes, double Reynolds, double time_step, Solution<double>* x_vel_previous_time, 
  Solution<double>* y_vel_previous_time, Solution<double>* T_prev_time, double heat_source, double specific_heat_graphite, 
  double specific_heat_water, double rho_graphite, double rho_water, double thermal_conductivity_graphite, double thermal_conductivity_water) 
  : WeakForm<double>(4), Stokes(Stokes), Reynolds(Reynolds), time_step(time_step), x_vel_previous_time(x_vel_previous_time), y_vel_previous_time(y_vel_previous_time)
  {
    BilinearFormSymVel* sym_form_0 = new BilinearFormSymVel(0, 0, Stokes, Reynolds, time_step);
    add_matrix_form(sym_form_0);
    BilinearFormSymVel* sym_form_1 = new BilinearFormSymVel(1, 1, Stokes, Reynolds, time_step);
    add_matrix_form(sym_form_1);

    BilinearFormUnSymVel_0_0* unsym_vel_form_0_0 = new BilinearFormUnSymVel_0_0(0, 0, Stokes);
    add_matrix_form(unsym_vel_form_0_0);
    BilinearFormUnSymVel_0_1* unsym_vel_form_0_1 = new BilinearFormUnSymVel_0_1(0, 1, Stokes);
    add_matrix_form(unsym_vel_form_0_1);
    BilinearFormUnSymVel_1_0* unsym_vel_form_1_0 = new BilinearFormUnSymVel_1_0(1, 0, Stokes);
    add_matrix_form(unsym_vel_form_1_0);
    BilinearFormUnSymVel_1_1* unsym_vel_form_1_1 = new BilinearFormUnSymVel_1_1(1, 1, Stokes);
    add_matrix_form(unsym_vel_form_1_1);

    BilinearFormUnSymXVelPressure* unsym_velx_pressure_form = new BilinearFormUnSymXVelPressure(0, 2);
    add_matrix_form(unsym_velx_pressure_form);

    BilinearFormUnSymYVelPressure* unsym_vely_pressure_form = new BilinearFormUnSymYVelPressure(1, 2);
    add_matrix_form(unsym_vely_pressure_form);

    VectorFormNS_0* F_0 = new VectorFormNS_0(0, Stokes, Reynolds, time_step);
    F_0->ext = Hermes::vector<MeshFunction<double>*>(x_vel_previous_time, y_vel_previous_time);
    add_vector_form(F_0);
    VectorFormNS_1* F_1 = new VectorFormNS_1(1, Stokes, Reynolds, time_step);
    F_1->ext = Hermes::vector<MeshFunction<double>*>(x_vel_previous_time, y_vel_previous_time);
    add_vector_form(F_1);
    VectorFormNS_2* F_2 = new VectorFormNS_2(2);
    add_vector_form(F_2);

    // Jacobian volumetric part.
    add_matrix_form(new BilinearFormTime(3, 3, "Outside", specific_heat_water, rho_water, time_step));
    add_matrix_form(new BilinearFormTime(3, 3, "Graphite Circle", specific_heat_graphite, rho_graphite, time_step));

    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(3, 3, "Outside", new Hermes1DFunction<double>(thermal_conductivity_water)));
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(3, 3, "Graphite Circle", new Hermes1DFunction<double>(thermal_conductivity_graphite)));
    add_matrix_form(new CustomJacobianAdvection(3, 3, "Outside"));

    // Residual - volumetric.
    VectorFormTime *vft = new VectorFormTime(3, "Outside", specific_heat_water, rho_water, time_step);
    Hermes::vector<MeshFunction<double>*> mesh_function_vector_outside;
    mesh_function_vector_outside.push_back(T_prev_time);
    vft->ext = mesh_function_vector_outside;
    add_vector_form(vft);

    vft = new VectorFormTime(3, "Graphite Circle", specific_heat_graphite, rho_graphite, time_step);
    Hermes::vector<MeshFunction<double>*> mesh_function_vector_graphite;
    mesh_function_vector_graphite.push_back(T_prev_time);
    vft->ext = mesh_function_vector_graphite;
    add_vector_form(vft);

    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(3, "Outside", new Hermes1DFunction<double>(thermal_conductivity_water)));
    add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<double>(3, "Graphite Circle", new Hermes1DFunction<double>(thermal_conductivity_graphite)));
    add_vector_form(new CustomResidualAdvection(3, "Outside"));
    
    add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(3, "Graphite Circle", new Hermes::Hermes2DFunction<double>(-heat_source)));
  };
