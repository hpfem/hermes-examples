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

/* Time-dependent Dirichlet condition for temperature */

class EssentialBCNonConst : public EssentialBoundaryCondition<double>
{
public:
  EssentialBCNonConst(std::string marker, double reactor_start_time, double temp_initial, double temp_reactor_max);
  
  ~EssentialBCNonConst() {};

  virtual EssentialBoundaryCondition<double>::EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const;

protected:
  double reactor_start_time;
  double temp_initial;
  double temp_reactor_max;
};
