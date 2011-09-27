#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Weak forms */

class CustomWeakFormAcoustics : public WeakForm<std::complex<double> >
{
public:
  CustomWeakFormAcoustics(std::string bdy_newton, double rho, double sound_speed, double omega);
};
