#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Weak forms */

class CustomWeakFormHeatMoistureRK : public WeakForm<double>
{
public:
  CustomWeakFormHeatMoistureRK(double c_TT, double c_ww, double d_TT, double d_Tw, double d_wT, double d_ww, 
      double k_TT, double k_ww, double T_ext, double w_ext, const std::string bdy_ext);
};
