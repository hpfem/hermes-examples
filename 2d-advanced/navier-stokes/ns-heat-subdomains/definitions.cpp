#include "definitions.h"

/* Custom initial condition for temperature*/

CustomInitialConditionTemperature::CustomInitialConditionTemperature(Mesh *mesh, double mid_x, double mid_y, double radius, double temp_water, double temp_graphite) 
  : ExactSolutionScalar<double>(mesh), mid_x(mid_x), mid_y(mid_y), radius(radius), temp_water(temp_water), temp_graphite(temp_graphite) {}

double CustomInitialConditionTemperature::value(double x, double y) const 
{
  bool in_graphite = (std::sqrt(sqr(mid_x - x) + sqr(mid_y - y)) < radius);
  double temp = temp_water;
  if (in_graphite) temp = temp_graphite;
  return temp;
}

void CustomInitialConditionTemperature::derivatives(double x, double y, double& dx, double& dy) const 
{   
  dx = 0;
  dy = 0;
}

Ord CustomInitialConditionTemperature::ord(Ord x, Ord y) const 
{
  return Ord(1);
}

/* Weak forms */

CustomWeakFormHeatAndFlow::CustomWeakFormHeatAndFlow(bool Stokes, double Reynolds, double time_step, Solution<double>* x_vel_previous_time, 
    Solution<double>* y_vel_previous_time, Solution<double>* T_prev_time, double heat_source, double specific_heat_graphite, 
    double specific_heat_water, double rho_graphite, double rho_water, double thermal_conductivity_graphite, double thermal_conductivity_water,
    bool simple_temp_advection) 
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
    add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<double>(3, 3, "Graphite", 
        new Hermes1DFunction<double>(thermal_conductivity_graphite/(rho_graphite * specific_heat_graphite)), HERMES_NONSYM));
    // Contribution from temperature advection - only in water.
    if (simple_temp_advection)     
    {
      CustomJacobianTempAdvection_3_3_simple* cjta_3_3_simple = new CustomJacobianTempAdvection_3_3_simple(3, 3, "Water");
      cjta_3_3_simple->ext.push_back(x_vel_previous_time);
      cjta_3_3_simple->ext.push_back(y_vel_previous_time);
      add_matrix_form(cjta_3_3_simple);
    }
    else 
    {
      add_matrix_form(new CustomJacobianTempAdvection_3_0(3, 0, "Water"));
      add_matrix_form(new CustomJacobianTempAdvection_3_1(3, 1, "Water"));
      add_matrix_form(new CustomJacobianTempAdvection_3_3(3, 3, "Water"));
    }

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
    if (simple_temp_advection)     
    {
      CustomResidualTempAdvection_simple* crta_simple = new CustomResidualTempAdvection_simple(3, "Water");
      crta_simple->ext.push_back(x_vel_previous_time);
      crta_simple->ext.push_back(y_vel_previous_time);
      add_vector_form(crta_simple);
    }
    else add_vector_form(new CustomResidualTempAdvection(3, "Water"));
  };


bool point_in_graphite(double x, double y)
{
  double dist_from_center = std::sqrt(sqr(x - HOLE_MID_X) + sqr(y - HOLE_MID_Y));
  if (dist_from_center < 0.5 * OBSTACLE_DIAMETER) return true;
  else return false;
}

int element_in_graphite(Element* e)
{
  // Calculate element center.
  int nvert;
  if (e->is_triangle()) nvert = 3;
  else nvert = 4;
  double elem_center_x = 0, elem_center_y = 0;
  for (int i=0; i < nvert; i++)
  {
    elem_center_x += e->vn[i]->x;
    elem_center_y += e->vn[i]->y;
  }
  elem_center_x /= nvert;
  elem_center_y /= nvert;
  // Check if center is in graphite.
  if (point_in_graphite(elem_center_x, elem_center_y)) 
  {
    return 0;  // 0... refine uniformly.
  }
  else 
  {
    return -1; //-1... do not refine.
  }
}

int element_in_water(Element* e)
{
  // Calculate element center.
  int nvert;
  if (e->is_triangle()) nvert = 3;
  else nvert = 4;
  double elem_center_x = 0, elem_center_y = 0;
  for (int i=0; i < nvert; i++)
  {
    elem_center_x += e->vn[i]->x;
    elem_center_y += e->vn[i]->y;
  }
  elem_center_x /= nvert;
  elem_center_y /= nvert;
  // Check if center is in graphite.
  if (point_in_graphite(elem_center_x, elem_center_y)) 
  {
    return -1;  //-1... do not refine.
  }
  else 
  {
    return 0;  // 0... refine uniformly.
  }
}
