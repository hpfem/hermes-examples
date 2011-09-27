#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Weak form */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(double K_squared);
};


