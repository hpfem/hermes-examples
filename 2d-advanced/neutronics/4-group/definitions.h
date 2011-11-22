////// Weak formulation in axisymmetric coordinate system  ////////////////////////////////////

#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;
using namespace Hermes::Hermes2D::WeakFormsNeutronics::Multigroup::CompleteWeakForms::Diffusion;
using namespace Hermes::Hermes2D::WeakFormsNeutronics::Multigroup::SupportClasses;

class CustomWeakForm : public DefaultWeakFormSourceIteration<double>
{
  public:
    CustomWeakForm(
        const Hermes::Hermes2D::WeakFormsNeutronics::Multigroup::MaterialProperties::Diffusion::MaterialPropertyMaps& matprop,
        Hermes::vector<MeshFunction<double>* >& iterates,
        double init_keff, std::string bdy_vacuum);
};

// Integral over the active core.
double integrate(MeshFunction<double>* sln, std::string area);
