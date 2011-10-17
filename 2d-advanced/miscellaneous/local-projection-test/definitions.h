#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string mat_al, Hermes1DFunction<double>* lambda_al,
                        std::string mat_cu, Hermes1DFunction<double>* lambda_cu,
                        Hermes2DFunction<double>* src_term);
};
