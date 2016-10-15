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

FilterFluxDensity::FilterFluxDensity(std::vector<MeshFunctionSharedPtr<double> > solutions)
  : DXDYFilter<double>(solutions)
{
}

Func<double>* FilterFluxDensity::get_pt_value(double x, double y, bool use_MeshHashGrid, Element* e)
{
	return new Func<double>();
}

MeshFunction<double>* FilterFluxDensity::clone() const
{
  std::vector<MeshFunctionSharedPtr<double> > fns;
  for (int i = 0; i < this->solutions.size(); i++)
    fns.push_back(this->solutions[i]->clone());
  return new FilterFluxDensity(fns);
}

void FilterFluxDensity::filter_fn(int n, double* x, double* y, const std::vector<const double *>& values, const std::vector<const double *>& dx, const std::vector<const double *>& dy, double* rslt, double* rslt_dx, double* rslt_dy)
{
	for (int i = 0; i < n; i++)
	{
		rslt[i] = std::sqrt(sqr(dy[0][i]) + sqr(dx[0][i]));
	}
}