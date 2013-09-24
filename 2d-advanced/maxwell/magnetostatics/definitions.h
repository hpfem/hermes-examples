#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::WeakFormsMaxwell;

/* Weak forms */

class CustomWeakFormMagnetostatics : public WeakForm<double>
{
public:
  CustomWeakFormMagnetostatics(std::string material_iron_1, std::string material_iron_2,
                               CubicSpline* mu_inv_iron, std::string material_air,
                               std::string material_copper, double mu_vacuum,
                               double current_density, int order_inc = 3);
};

class FilterVectorPotential : public Hermes::Hermes2D::MagFilter<double>
{
public:
  FilterVectorPotential(Hermes::vector<MeshFunctionSharedPtr<double> > solutions, Hermes::vector<int> items);

  virtual MeshFunction<double>* clone() const;
protected:
  void filter_fn(int n, Hermes::vector<double*> values, double* result, Geom<double> *e);
};

class FilterFluxDensity : public Hermes::Hermes2D::Filter<double>
{
public:
  FilterFluxDensity(Hermes::vector<MeshFunctionSharedPtr<double> > solutions);

  virtual Func<double>* get_pt_value(double x, double y, bool use_MeshHashGrid = false, Element* e = NULL);
  virtual MeshFunction<double>* clone() const;

protected:
  void precalculate(int order, int mask);
};
