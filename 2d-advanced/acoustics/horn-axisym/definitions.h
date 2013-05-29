#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

typedef std::complex<double> complex;

/* Weak forms */

class CustomWeakFormAcoustics : public WeakForm<complex>
{
public:
  CustomWeakFormAcoustics(std::string bdy_newton, double rho,
                          double sound_speed, double omega);
};
