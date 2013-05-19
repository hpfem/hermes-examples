#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::WeakFormsElasticity;
using namespace Hermes::Hermes2D::Views;

//#define USE_MULTICOMPONENT_FORMS

class CustomWeakFormLinearElasticity : public WeakForm<double>
{
public:
  CustomWeakFormLinearElasticity(double E, double nu, double rho_g,
      std::string surface_force_bdy, double f0, double f1);
};
