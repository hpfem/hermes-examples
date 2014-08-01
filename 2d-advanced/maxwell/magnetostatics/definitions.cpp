#include "definitions.h"

CustomWeakFormMagnetostatics::CustomWeakFormMagnetostatics(std::string material_iron_1, std::string material_iron_2,
  CubicSpline* mu_inv_iron, std::string material_air,
  std::string material_copper, double mu_vacuum,
  double current_density, int order_inc) : WeakForm<double>(1)
{
  // Jacobian.
  add_matrix_form(new DefaultJacobianMagnetostatics<double>(0, 0, std::vector<std::string>({ material_air, material_copper }),
    1.0, nullptr, HERMES_NONSYM, HERMES_AXISYM_Y, order_inc));
  add_matrix_form(new DefaultJacobianMagnetostatics<double>(0, 0, std::vector<std::string>({ material_iron_1, material_iron_2 }), 1.0,
    mu_inv_iron, HERMES_NONSYM, HERMES_AXISYM_Y, order_inc));
  // Residual.
  add_vector_form(new DefaultResidualMagnetostatics<double>(0, std::vector<std::string>({ material_air, material_copper }),
    1.0, nullptr, HERMES_AXISYM_Y, order_inc));
  add_vector_form(new DefaultResidualMagnetostatics<double>(0, std::vector<std::string>({ material_iron_1, material_iron_2 }), 1.0,
    mu_inv_iron, HERMES_AXISYM_Y, order_inc));
  add_vector_form(new DefaultVectorFormVol<double>(0, material_copper, new Hermes2DFunction<double>(-current_density * mu_vacuum)));
}

void FilterVectorPotential::filter_fn(int n, const std::vector<const double*>& values, double* result, GeomVol<double> *e)
{
  for (int i = 0; i < n; i++)
  {
    result[i] = 0;
    for (unsigned int j = 0; j < values.size(); j++)
      result[i] += sqr(values[j][i]);

    result[i] = std::sqrt(result[i]);
  }
}

MeshFunction<double>* FilterVectorPotential::clone() const
{
  std::vector<MeshFunctionSharedPtr<double> > fns;
  std::vector<int> items;
  for (int i = 0; i < this->solutions.size(); i++)
  {
    fns.push_back(this->solutions[i]->clone());
    items.push_back(this->items[i]);
  }
  return new FilterVectorPotential(fns, items);
}

FilterVectorPotential::FilterVectorPotential(std::vector<MeshFunctionSharedPtr<double> > solutions, std::vector<int> items)
  : MagFilter<double>(solutions, items)
{
}

FilterFluxDensity::FilterFluxDensity(std::vector<MeshFunctionSharedPtr<double> > solutions)
  : Filter<double>(solutions)
{
}

Func<double>* FilterFluxDensity::get_pt_value(double x, double y, bool use_MeshHashGrid, Element* e)
{
  throw Hermes::Exceptions::Exception("Not implemented yet"); return NULL;
}

MeshFunction<double>* FilterFluxDensity::clone() const
{
  std::vector<MeshFunctionSharedPtr<double> > fns;
  for (int i = 0; i < this->solutions.size(); i++)
    fns.push_back(this->solutions[i]->clone());
  return new FilterFluxDensity(fns);
}

void FilterFluxDensity::precalculate(int order, int mask)
{
  Quad2D* quad = quads[cur_quad];
  int np = quad->get_num_points(order, this->get_active_element()->get_mode());

  this->solutions[0]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
  this->solutions[1]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);

  const double *dudx1 = this->solutions[0]->get_dx_values();
  const double *dudy1 = this->solutions[0]->get_dy_values();

  const double *dudx2 = this->solutions[1]->get_dx_values();
  const double *dudy2 = this->solutions[1]->get_dy_values();

  const double *uval1 = this->solutions[0]->get_fn_values();
  const double *uval2 = this->solutions[1]->get_fn_values();

  update_refmap();
  double *x = this->refmap.get_phys_x(order);

  for (int i = 0; i < np; i++)
  {
    this->values[0][0][i] = std::sqrt(sqr(dudy1[i]) + sqr(dudy2[i]) +
      sqr(dudx1[i] + ((x[i] > 1e-10) ? uval1[i] / x[i] : 0.0)) +
      sqr(dudx2[i] + ((x[i] > 1e-10) ? uval2[i] / x[i] : 0.0)));
  }
}