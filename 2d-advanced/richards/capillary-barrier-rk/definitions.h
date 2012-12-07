#include "hermes2d.h"

#include "../constitutive.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Custom non-constant Dirichlet condition */

class RichardsEssentialBC : public EssentialBoundaryCondition<double> {
public:

  RichardsEssentialBC(std::string marker, double h_elevation, double pulse_end_time, double h_init, double startup_time) :
  EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), h_elevation(h_elevation), pulse_end_time(pulse_end_time), h_init(h_init), startup_time(startup_time)
  {
    markers.push_back(marker);
  }

  ~RichardsEssentialBC() {}

  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<double>::BC_FUNCTION; }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
    if (current_time < startup_time)
      return h_init + current_time/startup_time*(h_elevation-h_init);
    else if (current_time > pulse_end_time)
      return h_init;
    else
      return h_elevation;
  }

  // Member.
  double h_elevation;
  double pulse_end_time;
  double h_init;
  double startup_time;
};

/* Weak forms */

class CustomWeakFormRichardsRK : public WeakForm<double>
{
public:
  CustomWeakFormRichardsRK(ConstitutiveRelations* constitutive);

private:

  class CustomJacobianFormVol : public MatrixFormVol<double>
  {
  public:
    CustomJacobianFormVol(int i, int j, ConstitutiveRelations* constitutive)
      : MatrixFormVol<double>(i, j), constitutive(constitutive)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
      Func<double> *v, Geom<double> *e, Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
      Func<Ord> *v, Geom<Ord> *e, Func<Ord>* *ext) const;

    virtual MatrixFormVol<double>* clone() const;
    ConstitutiveRelations* constitutive;
  };

  class CustomResidualFormVol : public VectorFormVol<double>
  {
  public:
    CustomResidualFormVol(int i, ConstitutiveRelations* constitutive)
      : VectorFormVol<double>(i), constitutive(constitutive)
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e,
      Func<double>* *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
      Func<Ord>* *ext) const;

    virtual VectorFormVol<double>* clone() const;
    ConstitutiveRelations* constitutive;
  };
};
